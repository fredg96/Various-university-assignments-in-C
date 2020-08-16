#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define MASTER 0 
void prop(int *x, double *w);

int main(int argc, char *argv[]){
    if (argc != 3)
    {
        printf("Program requires 2 arguments: N and output_file\n");
        return -1;
    }
	int N = atoi(argv[1]);
    char *output_file_name = argv[2];

    int rank, size, chunk, R, minMax[1], absvals[1], bins, hist_idx, loc_hist[20], hist[20];
	int *x = NULL, *loc_res = NULL, *minMaxs= NULL, *hist_res = NULL;
	double t, t_start, t_end, T, interval, tau;
	double *w = NULL, *P = NULL, *exe_time = NULL;

    MPI_Init(&argc, &argv); 
    MPI_Status status;
    MPI_Request req;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of PEs
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get local rank  

	/*Check that N is divisible by rank, ensures load balance*/
	if(rank == MASTER){
		if(N % size != 0){
			printf("Number of simulations not divisible by number of PEs, exiting\n");
			return -1;
		}
	}
	
	/* Initialize data and set up variables */
	chunk = N / size; // How many experiments will each PE do 
    x = (int *)malloc(7 * sizeof(int)); //state vector
    w = (double *)malloc(15 * sizeof(double)); //vector with derived properties
	loc_res = (int *)malloc(chunk*sizeof(int)); //vector to store local results from MC experiments
	P = (double *)malloc(15 * 7 * sizeof(double)); //Initialize state update matrix

	T = 100; //end time
	bins = 20; // number of bins
	minMax[0] = 200000; // set minimum value to be large to ensure that we find a smaller value
 	minMax[1] = 0; // set maimum value small to ensure that we find a larger value
	int loc_seed = (int)time(NULL);
    srand (loc_seed+rank); //try to ensure each PE has unique seed 
	
	/*Set up state update matrix*/
	R = 15;
	for(int i=0; i<15*7; i++){
		P[i] = 0.0;
	}
	P[0*7+0] = 1.0;
	P[1*7+0] = -1.0;
	P[2*7+0] = -1.0;
	P[2*7+2] = 1.0;
	P[3*7+1] = 1.0;
	P[4*7+1] = -1.0;
	P[5*7+1] = -1.0;
	P[5*7+3] = 1.0;
	P[6*7+2] = -1.0;
	P[7*7+2] = -1.0;
	P[7*7+4] = 1.0;
	P[8*7+3] = -1.0;
	P[9*7+3] = -1.0;
	P[9*7+5] = 1.0;
	P[10*7+4] = -1.0;
	P[11*7+4] = -1.0;
	P[11*7+6] = 1.0;
	P[12*7+5] = -1.0;
	P[13*7+0] = 1.0;
	P[13*7+6] = -1.0;
	P[14*7+6] = -1.0;
	
	t_start = MPI_Wtime(); // Input finished start timer 

	/*Do experiments and collect data */
	for(int k = 0; k<chunk; k++){	
		t = 0.0;
		/*Set state vector to initial values*/
		x[0] = 900;
		x[1] = 900;
		x[2] = 30;
		x[3] = 330;
		x[4] = 50;
		x[5] = 270;
		x[6] = 20;
		while(t<T){
			prop(x,w); // Line3 in Algorithm 1

			/*Line 4 in Algorithm 1*/
			double a0 = 0.0; 
			for(int i=0; i<15; i++){
				a0 += w[i];
			}

			/*Line 5 in Algorithm 1*/
			double u1 =  rand() / ((double) RAND_MAX);
			double u2 =  rand() / ((double) RAND_MAX);
			double thresh = a0*u2;

			tau = -log(u1)/a0; // Line 6 in Algorithm 1

			/*Line 7 in Algorithm 1*/
			double curr_sum = 0.0;
			int r = 0;
			int idx = 0;
			while(curr_sum<thresh && idx < 15){
				curr_sum += w[idx];
				idx += 1;
			}
			r = idx - 1;

			/*Line 8 in Algorithm 1*/
			x[0] = x[0]+P[r*7+0];
			x[1] = x[1]+P[r*7+1];
			x[2] = x[2]+P[r*7+2];
			x[3] = x[3]+P[r*7+3];
			x[4] = x[4]+P[r*7+4];
			x[5] = x[5]+P[r*7+5];
			x[6] = x[6]+P[r*7+6];
				
			t += tau; // Line 9 in Algorithm 1
					}
		loc_res[k] = x[0]; //store result of MC experiment
		if(x[0]<minMax[0]){ // check if current results is the smallest and store it
			minMax[0] = x[0];
		}
		if(x[0]>minMax[1]){ // check if current result is the largest and store it
			minMax[1] = x[0];
		}
	}

	/*Start taking meassurements for histogram, min max value and frequency*/
	for(int i=0; i<bins; i++){ // set each bin to zero for each PE
		loc_hist[i] = 0;
	}

	minMaxs = (int *)malloc(2*size * sizeof(int)); // allocate to store everyones largest and smalest value	
	MPI_Allgather(&minMax, 2, MPI_INT, minMaxs, 2, MPI_INT, MPI_COMM_WORLD); // gather smallest/largest value from everyone
	for(int i=0; i<size;i++){ //Find largest and smallest value
		if(minMax[0]>minMaxs[2*i]){ //smalest values are on even positions, 0, 2, ...
			minMax[0]=minMaxs[2*i];
		}
		if(minMax[1]<minMaxs[2*i+1]){ //largest values are on odd positions 1,3,...
			minMax[1]=minMaxs[2*i+1];
		}	
	}

	interval = (double)(minMax[1]-minMax[0])/bins; // calculate width of bins
	hist_idx = 0; // index for histogram [0,bins]
	for(int i = 0; i<chunk; i++){
		if(loc_res[i] != minMax[1]){
			hist_idx = ((loc_res[i]-minMax[0])/interval); // calculate bin index idx for results from experiment i
		}
		else{
			hist_idx = bins-1; //Handle upper edge case and place them/it in last bin ensure that the highest value get placed in bin with index bins-1
		}
		loc_hist[hist_idx] += 1;
	}
	if(rank!=MASTER){
		MPI_Send(loc_hist, bins, MPI_INT, MASTER, 666, MPI_COMM_WORLD); //send each local histogram, i.e. the frequency for each bin
	}
	else{
		hist_res = (int *)malloc(size*bins * sizeof(double)); //allocate memory to store each local histogram
		memcpy(hist_res, loc_hist, bins * sizeof(int)); //master copies own local histogram to global
		for (int i = 1; i < size; i++){ // receive each local histogram
            		MPI_Recv(hist_res + (i*bins), bins, MPI_INT, i, 666, MPI_COMM_WORLD, &status);
		}
		for(int i= 0;i<bins;i++){ //make sure that global histogram starts with eachbin having zero samples
			hist[i] = 0;
		}
		for(int i=0;i<size;i++){ //calculate global histogram from local
			for(int j=0;j<bins;j++){
				hist[j] = hist[j]+hist_res[i*bins+j];
			}
		}
		free(hist_res);
	}
	
	t_end = MPI_Wtime(); // End timer before writing output

	if (rank != MASTER){ //since most are not master	
		MPI_Ssend(&t_end, 1, MPI_DOUBLE, 0, 1349*rank, MPI_COMM_WORLD); //non master send their end time to master

		/*Everyone except master freees own data*/
		free(x);
		free(w);
		free(P);
		free(loc_res);
		free(minMaxs);
    }
	else{
		double worst_time = 0.0; 
		double slowest_hist = 0.0;
		exe_time = (double *)malloc(size*sizeof(double)); //to store everyones end times
		exe_time[0] = t_end; // set own end time
		for(int i=1; i<size; i++){
			MPI_Recv(&exe_time[i], 1, MPI_DOUBLE, i, 1349*i, MPI_COMM_WORLD, &status); //receive end time from everyone	
		}
		
		/*find slowest PEs end time*/
		for(int i=0; i<size; i++){ 
			if(exe_time[i] > worst_time){ 
				worst_time = exe_time[i]; 
			}	
		}
        printf("%f\n", worst_time-t_start);
        FILE *fp = fopen(output_file_name, "a");
        if (fp == NULL){
            printf("File Error\n");
            exit(1);
        }
		//fprintf(fp, "%f ", worst_time-t_start);
		fprintf(fp, "%d ", minMax[0]); //write smallest, lower bound, value to file
		fprintf(fp, "%d ", minMax[1]); //write largest, upper bound, value to file
		for (int i = 0; i < bins; i++){
			fprintf(fp, "%d ", hist[i]); //write number of samples for each bin in order 0 to 19 for 20 bins in total
		}
        fclose(fp);

		/* master frees own data*/
		free(x);
		free(w);
		free(P);
		free(loc_res);	
		free(minMaxs);
	}
    MPI_Finalize(); /* Shut down and clean up MPI */
    return 0;
}
/*Function provided on studentportalen*/
/**
 * Compute propensities for the Malaria model.
 * @param x State vector Should be of length 7!
 * @param w Result vector (propensities). Should be of length 15!
 *
 */
void prop(int *x, double *w) {
	// Birth number, humans
	const double LAMBDA_H = 20;
	// Birth number, mosquitoes
	const double LAMBDA_M = 0.5;
	// Biting rate of mosquitoes
	const double B = 0.075;
	/* Probability that a bite by an infectious mosquito results in transmission
	   of disease to human*/
	const double BETA_H = 0.3;
	/* Probability that a bite results in transmission of parasite to a
	   susceptible mosquito*/
	const double BETA_M = 0.5;
	// Human mortality rate
	const double MU_H = 0.015;
	// Mosquito mortality rate
	const double MU_M = 0.02;
	// Disease induced death rate, humans
	const double DELTA_H = 0.05;
	// Disease induced death rate, mosquitoes
	const double DELTA_M = 0.15;
	// Rate of progression from exposed to infectious state, humans
	const double ALFA_H = 0.6;
	// Rate of progression from exposed to infectious state, mosquitoes
	const double ALFA_M = 0.6;
	// Recovery rate, humans
	const double R = 0.05;
	// Loss of immunity rate, humans
	const double OMEGA = 0.02;
	/* Proportion of an antibody produced by human in response to the incidence
	   of infection caused by mosquito. */
	const double NU_H = 0.5;
	/* Proportion of an antibody produced by mosquito in response to the
	   incidence of infection caused by human. */
	const double NU_M = 0.15;

	w[0] = LAMBDA_H;
	w[1] = MU_H * x[0];
	w[2] = (B * BETA_H * x[0] * x[5]) / (1 + NU_H * x[5]);
	w[3] = LAMBDA_M;
	w[4] = MU_M * x[1];
	w[5] = (B * BETA_M * x[1]*x[4]) / (1 + NU_M * x[4]);
	w[6] = MU_H * x[2];
	w[7] = ALFA_H * x[2];
	w[8] = MU_M * x[3];
	w[9] = ALFA_M * x[3];
	w[10] = (MU_H + DELTA_H) * x[4];
	w[11] = R * x[4];
	w[12] = (MU_M + DELTA_M) * x[5];
	w[13] = OMEGA * x[6];
	w[14] = MU_H * x[6];
}