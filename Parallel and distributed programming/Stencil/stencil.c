#include "stencil.h"
#include <string.h>

#define MASTER 0

int main(int argc, char* argv[]){
    if(argc != 4){
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
    }

    /* MPI initialization */
    MPI_Status status; 
    
    char* input_name = argv[1]; 
    char* output_name = argv[2]; 
    int num_steps = atoi(argv[3]); 


    int rank, size, left, right, num_values; 
    double* pIn, *pOut, *data; 
    
    MPI_Init(NULL, NULL); 
    double leftpData[2], rightpData[2]; 

    /* Assign local variables */
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    left = (rank-1+size)%size; 
    right = (rank+1+size)%size; 

    MPI_Request lsendr, rsendr, lrecvr, rrecvr; 

    /* Master PE read data*/
    if(rank == MASTER){
		num_values = read_input(input_name, &data); 
		if(num_values < 0){
			return 2; 
		}
    }

    /* Broadcast number of values to allow allocation of memory*/
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    /* Stencil settings */
    double h = 2.0*PI/num_values;
    const int STENCIL_WIDTH = 5;
    const int EXTENT = STENCIL_WIDTH/2;
    const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

    /* Allocate memory for local data */
    pIn = (double*)malloc(num_values/size * sizeof(double)); 
    pOut = (double*)malloc(num_values/size * sizeof(double)); 
    if(pIn == NULL || pOut == NULL){
		printf("Process %d: Could not allocate memory! \n", rank); 
    }

    /* Distribute data */
    if(rank == MASTER){
		memcpy(pIn, data, num_values/size*sizeof(double)); 
	/*Master sends data*/
		for(int i=1; i<size; i++){
			MPI_Ssend(&data[i*num_values/size], num_values/size, MPI_DOUBLE, i, 111*i, MPI_COMM_WORLD); 
		}
    }
    /* PEs recive data*/
    else{
		MPI_Recv(pIn, num_values/size, MPI_DOUBLE, 0, 111*rank, MPI_COMM_WORLD, &status); 
    }
    double start, my_execution_time; 
    double* exe_time;  
	start = MPI_Wtime(); 
    for(int step=0; step<num_steps; step++){
		/*Get neigbhors data*/
		MPI_Isend(pIn, EXTENT, MPI_DOUBLE, left, 1*(rank+1), MPI_COMM_WORLD, &lsendr); 
		MPI_Isend(&pIn[num_values/size-2], EXTENT, MPI_DOUBLE, right, 1111*(rank+1), MPI_COMM_WORLD, &rsendr); 
		MPI_Irecv(leftpData, EXTENT, MPI_DOUBLE, left, 1111*(left+1), MPI_COMM_WORLD, &lrecvr); 
		MPI_Irecv(rightpData, EXTENT, MPI_DOUBLE, right, 1*(right+1), MPI_COMM_WORLD, &rrecvr); 

		MPI_Wait(&lsendr, &status); 
		MPI_Wait(&rsendr, &status); 
		MPI_Wait(&lrecvr, &status); 
		MPI_Wait(&rrecvr, &status);

		/* Apply stencil on local data */
		double result = 0; 
		for(int i=0; i<EXTENT; i++){
			result=0; 
			for(int j=0; j<STENCIL_WIDTH; j++){
				double base = 0; 
				if(i+j-EXTENT<0){ 
					base = leftpData[i+j]; 
				}
				else{
					base = pIn[i+j-EXTENT]; 
					}
				result += STENCIL[j] * base;
			}
			pOut[i] = result; 
		}

		for(int i=EXTENT; i<num_values/size-EXTENT; i++){
			result = 0; 
			for(int j=0; j<STENCIL_WIDTH; j++){
				result += STENCIL[j] * pIn[i+j-EXTENT]; 
			}
			pOut[i] = result; 
		}

		for(int i=num_values/size-EXTENT; i<num_values/size; i++){
			result = 0; 
			for(int j=0; j<STENCIL_WIDTH; j++){
				double base = 0; 
				if(i+j-EXTENT>=num_values/size){ 
					base = rightpData[i+j-EXTENT-num_values/size]; 
				}
				else{ 
					base = pIn[i+j-EXTENT]; 
				}
				result += STENCIL[j] * base; 
			}
			pOut[i] = result; 
		}
		/* Swap input and output */
		if(step < num_steps-1){
		   double* temp = pIn; 
		   pIn = pOut; 
		   pOut = temp; 
		}
    }
	my_execution_time = MPI_Wtime() - start; 

    /* Master collects data and execution time from other PEs, report highest time and writes data */ 
	if(rank == MASTER){
		double worst_time = 0.0; 
		exe_time = malloc(size*sizeof(double)); 
		exe_time[0] = my_execution_time; 
		memcpy(data, pOut, num_values/size*sizeof(double)); 
		for(int i=1; i<size; i++){
			MPI_Recv(&data[i*num_values/size], num_values/size, MPI_DOUBLE, i, 111*i, MPI_COMM_WORLD, &status);
			MPI_Recv(&exe_time[i], 1, MPI_DOUBLE, i, 1349*i, MPI_COMM_WORLD, &status); 
			}
		for(int i=0; i<size; i++){
			if(exe_time[i] > worst_time){ 
				worst_time = exe_time[i]; 
			}
		}
		printf("%f\n", worst_time); 	
		if(write_output(output_name, data, num_values) != 0) {
			return 2;
		}
		}
	else{
	/* Other PEs send their data and execution time */
		MPI_Ssend(pOut, num_values/size, MPI_DOUBLE, 0, 111*rank, MPI_COMM_WORLD); 
		MPI_Ssend(&my_execution_time, 1, MPI_DOUBLE, 0, 1349*rank, MPI_COMM_WORLD); 
	}	
    /* clean */
    free(pIn); 
    free(pOut); 
    if(rank == MASTER){
		free(data);
		free(exe_time); 
    } 
    MPI_Finalize();  	
}

int read_input(const char *file_name, double **values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "r"))) {
        perror("Couldn't open input file");
        return -1;
    }
    int num_values;
    if (EOF == fscanf(file, "%d", &num_values)) {
        perror("Couldn't read element count from input file");
        return -1;
    }
    if (NULL == (*values = malloc(num_values * sizeof(double)))) {
        perror("Couldn't allocate memory for input");
        return -1;
    }
    for (int i=0; i<num_values; i++) {
        if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
	    perror("Couldn't read elements from input file");
	    return -1;
	    }
    }
    if (0 != fclose(file)){
        perror("Warning: couldn't close input file");
    }
    return num_values;
}

int write_output(char *file_name, const double *output, int num_values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "w"))) {
	perror("Couldn't open output file");
	return -1;
    }
    for (int i = 0; i < num_values; i++) {
	if (0 > fprintf(file, "%.4f ", output[i])) {
	    perror("Couldn't write to output file");
	}
    }
    if (0 > fprintf(file, "\n")) {
    	perror("Couldn't write to output file");
    }
    if (0 != fclose(file)) {
    	perror("Warning: couldn't close output file");
    }
    return 0;
}