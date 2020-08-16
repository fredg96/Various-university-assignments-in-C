#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MASTER 0 

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Program requires 2 arguments: input_file and output_file\n");
        return -1;
    }

    char *input_file_name = argv[1];
    char *output_file_name = argv[2];

    int n, rank, size, chunk;
    double *A = NULL, *B = NULL, *loc_C = NULL, *loc_A = NULL, *loc_B = NULL, *tmp = NULL;
    double t_start, t_end; 
	double* exe_time;

    MPI_Init(&argc, &argv); 
    MPI_Status status;
    MPI_Request req;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of PEs
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get my rank                  

    /* Master reads data */
    if (rank == MASTER){
        FILE *fp;
        fp = fopen(input_file_name, "r");
        if (fp == NULL){
            printf("Error finding file\n");
            exit(1);
        }
        fscanf(fp, "%d", &n);
        if (n % size != 0){
            printf("Can't divide N by number of processors, exiting program.");
            exit(1);
        }

        A = (double *)malloc(n * n * sizeof(double));
        B = (double *)malloc(n * n * sizeof(double));

        for (int i = 0; i < n*n; i++){ //Read A row wise    
            fscanf(fp, "%lf", &A[i]);
        }
        for (int i = 0; i < n; i++){ //Read B column wise
            for (int j = 0; j < n; j++){
                fscanf(fp, "%lf", &B[j * n + i]);
            }
        }

        fclose(fp);
    }

    
    t_start = MPI_Wtime(); // Input finished start timer to meassure all communication

    MPI_Bcast(&n, 1, MPI_INT, MASTER, MPI_COMM_WORLD); // broadcast to other processes
    chunk = n / size; // Data for each PE 
    loc_A = (double *)malloc(chunk * n * sizeof(double));
    loc_B = (double *)malloc(chunk * n * sizeof(double));

    MPI_Scatter(A, chunk * n, MPI_DOUBLE, loc_A, chunk * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(B, chunk * n, MPI_DOUBLE, loc_B, chunk * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    tmp = (double *)malloc(chunk * n * sizeof(double));
    memcpy(tmp, loc_B, chunk * n * sizeof(double));

    int left = ((rank - 1) % size + size) % size;
    int right = (rank + 1) % size;
    int offset = chunk * rank;
    loc_C = (double *)malloc(chunk * n * sizeof(double));

	/* Calculate local elemet of C */
    for (int i = 0; i < size; i++){
        if (i > 0){
            offset = ((left + 1) % size) * chunk;
        }
        for (int row = 0; row < chunk; row++){
            for (int col = 0; col < chunk; col++){
                double sum = 0.0;
                for (int k = 0; k < n; k++){
                    sum += loc_A[row * n + k] * tmp[col * n + k];
                }
                loc_C[row * n + (col + offset)] = sum;
            }
        }
        if (i != size - 1){
            MPI_Isend(loc_B, chunk * n, MPI_DOUBLE, right, 666, MPI_COMM_WORLD, &req);
            MPI_Isend(&offset, 1, MPI_INT, right, 999, MPI_COMM_WORLD, &req);
            MPI_Recv(tmp, chunk * n, MPI_DOUBLE, left, 666, MPI_COMM_WORLD, &status);
            MPI_Recv(&offset, 1, MPI_INT, left, 999, MPI_COMM_WORLD, &status);
        }
        right = (right + 1) % size;
        left = ((left - 1) % size + size) % size;
    }

    double *C = NULL;
    if (rank == MASTER){
        C = (double *)malloc(n * n * sizeof(double));
	}
    MPI_Gather(loc_C, chunk * n, MPI_DOUBLE, C, chunk * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD); // Master gather results

    t_end = MPI_Wtime(); // End timer before writing output

    // Write data to file
    if (rank == MASTER){
		double worst_time = 0.0; 
		exe_time = malloc(size*sizeof(double)); 
		exe_time[0] = t_end;
		for(int i=1; i<size; i++){
			MPI_Recv(&exe_time[i], 1, MPI_DOUBLE, i, 1349*i, MPI_COMM_WORLD, &status); 
			}
		for(int i=0; i<size; i++){
			if(exe_time[i] > worst_time){ 
				worst_time = exe_time[i]; 
			}
		}
        printf("%.2f\n", worst_time);

        FILE *fp = fopen(output_file_name, "a");
        if (fp == NULL){
            printf("File Error\n");
            exit(1);
        }

        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                if (j == n - 1){
                    fprintf(fp, "%f\n", C[i * n + j]);
				}
                else{
                    fprintf(fp, "%f ", C[i * n + j]);
				}
            }
        }
        fclose(fp);
        free(C);
    }
	else{
		MPI_Ssend(&t_end, 1, MPI_DOUBLE, 0, 1349*rank, MPI_COMM_WORLD); 
	}

	free(A);
    free(B);
    free(loc_A);
    free(loc_B);
    free(tmp);
    free(loc_C);

    MPI_Finalize(); /* Shut down and clean up MPI */

    return 0;
}