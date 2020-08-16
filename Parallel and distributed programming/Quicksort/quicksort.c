#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MASTER 0 /* task ID of master task */

int cmpfunc (const void * a, const void * b){
   return ( *(int*)a - *(int*)b );
}

int main(int argc, char *argv[]){
    if (argc != 4){
        printf("Give 3 arguments: name of input, name of output, pivot strategy 1-3:\n");
        return 0;
    }

    char *input_file_name = argv[1];
    char *output_file_name = argv[2];
    int pivot_strategy = atoi(argv[3]);

    int rank, size;
    int n;
    int *data = NULL;
    int chunk;             
    int i, loc_start, loc_end;  
    double t_start, t_end; 

    MPI_Init(&argc, &argv); 
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

    /*Master reads data*/
    if (rank == MASTER){
        FILE *f;
        f = fopen(input_file_name, "r");
        if (f == NULL)
        {
            printf("File Error\n");
            exit(1);
        }
        fscanf(f, "%d", &n);
        data = (int *)malloc(n * sizeof(int));
        for (int i = 0; i < n; i++)
        {
            fscanf(f, "%d", &data[i]);
        }
        fclose(f);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != MASTER){
        data = (int *)malloc(n * sizeof(int));
    }
    MPI_Bcast(data, n, MPI_INT, 0, MPI_COMM_WORLD);

    chunk = n / size; // number per PE
    loc_start = rank * chunk;          
    loc_end = (rank + 1) * chunk - 1; 
    if (rank == size - 1){
        loc_end = n - 1; // last PE compute until end
	}
    t_start = MPI_Wtime();

    int loc_rank, loc_size;
    MPI_Comm group_comm;
    MPI_Comm last_comm = MPI_COMM_WORLD;
    loc_rank = rank;
    loc_size = size;

    /* new local array*/
    int loc_n = loc_end - loc_start + 1;
    int *loc_arr = NULL;
    loc_arr = (int *)malloc(loc_n * sizeof(int));
    memcpy(loc_arr, data + loc_start, loc_n * sizeof(int));

    qsort(loc_arr, loc_n, sizeof(int), cmpfunc);

    while (loc_size > 1){
        int pivot;
        if (pivot_strategy == 1){ // median of on PE, master         
            if (loc_rank == 0)
                pivot = loc_arr[loc_n / 2];
        }
        else if (pivot_strategy == 2){ // median of all medians in each group          
            pivot = loc_arr[loc_n / 2];
            if (loc_rank == 0){
				int *pivots = NULL; // array to store medians
                pivots = (int *)malloc(loc_size * sizeof(int));
                pivots[0] = pivot;
                for (int i = 1; i < loc_size; i++)                {
                    MPI_Recv(pivots + i, 1, MPI_INT, i, 111, last_comm, &status);
                }
                qsort(pivots, loc_size, sizeof(int), cmpfunc );
                pivot = pivots[loc_size / 2]; 
               
            }
            else{
				MPI_Send(&pivot, 1, MPI_INT, 0, 111, last_comm);
            }          
        }
        else if (pivot_strategy == 3){ // mean value of all medians in each group         
            pivot = loc_arr[loc_n / 2];
            if (loc_rank == 0){
				int *pivots = NULL; // array to store medians
                pivots = (int *)malloc(loc_size * sizeof(int));
                pivots[0] = pivot;
                for (int i = 1; i < loc_size; i++){
                    MPI_Recv(pivots + i, 1, MPI_INT, i, 111, last_comm, &status);
                }
                unsigned long int sum_median = 0;
                for (int i = 0; i < loc_size; i++){
                    sum_median += pivots[i];
                }
                pivot = sum_median / loc_size; // mean of medians
                
            }
            else{
				MPI_Send(&pivot, 1, MPI_INT, 0, 111, last_comm);
            }
        }
        
        MPI_Bcast(&pivot, 1, MPI_INT, 0, last_comm);
        
        /* find the spliting position */
        int pivot_index = loc_n / 2;
        while (pivot < loc_arr[pivot_index] && pivot_index > 0){
            pivot_index--;
		}
        while (pivot > loc_arr[pivot_index] && pivot_index < loc_n){
            pivot_index++;
		}
        int n_send, n_keep, n_recv;
        int *kept = NULL;
        int *received = NULL;

        if (loc_rank<(loc_size/2)){
            n_send = loc_n - pivot_index;
            n_keep = pivot_index;
			/*send large part, keep small part*/
            MPI_Send(&n_send, 1, MPI_INT, loc_rank + (loc_size/2), 0, last_comm);
			MPI_Send(loc_arr + pivot_index, n_send, MPI_INT, loc_rank + (loc_size/2), 1, last_comm); 
			/*receive data*/
            MPI_Recv(&n_recv, 1, MPI_INT, loc_rank + (loc_size/2), 0, last_comm, &status);
            received = (int *)malloc(n_recv * sizeof(int));
			MPI_Recv(received, n_recv, MPI_INT, loc_rank + (loc_size/2), 1, last_comm, &status);
        
            kept = (int *)malloc(n_keep * sizeof(int));
            memcpy(kept, loc_arr, n_keep * sizeof(int));
        }
		else{
            n_send = pivot_index;
            n_keep = loc_n - pivot_index;
			/*receive data*/
            MPI_Recv(&n_recv, 1, MPI_INT, loc_rank - (loc_size/2), 0, last_comm, &status);
			received = (int *)malloc(n_recv * sizeof(int));
			MPI_Recv(received, n_recv, MPI_INT, loc_rank - (loc_size/2), 1, last_comm, &status);
			/*send small part, keep large part*/
            MPI_Send(&n_send, 1, MPI_INT, loc_rank - (loc_size/2), 0, last_comm);                                    
            MPI_Send(loc_arr, n_send, MPI_INT, loc_rank - (loc_size/2), 1, last_comm);
            kept = (int *)malloc(n_keep * sizeof(int));
            memcpy(kept, loc_arr + pivot_index, n_keep * sizeof(int));
        }

        /* update local array by merging in a sorted way*/
        free(loc_arr);
        loc_n = n_keep + n_recv; // local array size
        loc_arr = NULL;
        loc_arr = (int *)malloc(loc_n * sizeof(int));
		/*Merge array in sorted order*/
		int i = 0, j = 0, k = 0;
		while (i < n_keep && j < n_recv){ //compare between arrays
			if (kept[i] < received[j]){
				loc_arr[k++] = kept[i++];
			}
			else{
				loc_arr[k++] = received[j++];
			}
		}
		while (i < n_keep){ // fill in remaining from kept part
			loc_arr[k++] = kept[i++];
		}
		while (j < n_recv){ // fill in remaining from received part
			loc_arr[k++] = received[j++];
		}

        free(received);
        free(kept);

        MPI_Comm_split(last_comm, loc_rank<(loc_size/2), loc_rank, &group_comm);
        MPI_Comm_rank(group_comm, &loc_rank);
        MPI_Comm_size(group_comm, &loc_size);
        last_comm = group_comm;
    }
    /*MASTER merges local arrays into whole array*/
	if (rank == MASTER){
		free(data);
        data = NULL;
        data = (int *)malloc(n * sizeof(int));
        memcpy(data, loc_arr, loc_n * sizeof(int));
        int *loc_arr_sizes = NULL;
        loc_arr_sizes = (int *)malloc(size * sizeof(int));

        for (int i = 1; i < size; i++){
            MPI_Recv(loc_arr_sizes + i, 1, MPI_INT, i, 999, MPI_COMM_WORLD, &status);
            MPI_Recv(data + loc_n, loc_arr_sizes[i], MPI_INT, i, 666, MPI_COMM_WORLD, &status);
            loc_n += loc_arr_sizes[i];
        }
        free(loc_arr_sizes);
        
    }
    else{
        MPI_Send(&loc_n, 1, MPI_INT, MASTER, 999, MPI_COMM_WORLD);
        MPI_Send(loc_arr, loc_n, MPI_INT, MASTER, 666, MPI_COMM_WORLD);
    }

    t_end = MPI_Wtime();

    /*Master writes data*/
    if (rank == MASTER){
        printf("%.2f\n", t_end - t_start);
        
        for (int i = 0; i < n - 1; i++){
            if (data[i] > data[i + 1]){
                printf("Wrong!\n");
                break;
            }
        }
        FILE *fp = fopen(output_file_name, "a");
        if (fp == NULL){
            printf("File Error\n");
            exit(1);
        }
        for (int i = 0; i < n; i++){
            fprintf(fp, "%d  ", data[i]);
        }
        fclose(fp);
    }
    free(loc_arr);
    free(data);
    MPI_Finalize(); /* Shut down and clean up MPI */

    return 0;
}