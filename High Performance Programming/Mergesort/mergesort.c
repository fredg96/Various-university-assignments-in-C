/* 
*Program thet generates an array with random doubles and sorts it with parallel mergesort
*Purpose: Individual assignment in the course high performance programing (Uppsala University)
*Author: Fredrik Gutafsson
*Date of last modification: 21-03-2019
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void print(int size, double array[]);
int main(int argc, char** argv);
void sort(double array[], int left, int right, int threads, double sorted[]);
void combine(double array[], int leftIndex, int midIndex, int rightIndex, double sorted[]);
void serialSort(double array[], int left, int right, double sorted[]);
void run_omp (double array[], int size, double sorted[], int nThreads);
double randomNum();

int size;
int nThreads;
/* 
*Main function to generate an array with psudorandom numbers and then sort it with mergesort
*Inputs: Size of array, number of threads, if printing
*Returns: sorted array
*/
int main(int argc, char** argv){
    if(argc!= 4){
        printf("Incorrect amount of input argument. \n Size of array, number of threads and if printing(1, printing, or 0, not printing) expected.");
        return 1; 
    }
    size = atoi(argv[1]); // size of array
    nThreads = atoi(argv[2]); //number of threads
    if(nThreads<0){
        printf("Negative amount of threads not possible \n");
        return 1;
    }
    int iPrint = atoi(argv[3]); //if printing iprint = 1
    double* array = (double *)malloc(sizeof(double)*size);
    double *sorted = (double *)malloc(sizeof(double)*size); //Array for merging
    int i;
    srand(34241243); //Set a seed for reproducibility
    for(i= 0; i<size;i++){ //Fill array with pseudorandom numbers
        array[i] = randomNum();
    }
    iPrint != 1?:print(size,array);  //if printing iprint = 1
    printf("\n Sorting \n");   
    omp_set_nested (1);
    sort(array, 0, size-1, nThreads, sorted);  //parallel mergesort  
    iPrint != 1?:print(size,array);
    free(array); //Free array
    free(sorted); //Free copying array
    return 0;
}
/* 
*Prints array of size N
*Inputs: size of array, array to print
*/
void print(int size, double array[]){
    int i;
    for(i = 0; i<size; i++){
        printf("%f ", array[i]);
    }
}

/* 
*Recursively sort array in parallel 
*Inputs: unsorted array, left index, right index, number of available threads, dummy array to sort into
*/
void sort(double array[], int left, int right, int threads, double sorted[]){
     if(left<right){ //Check if we have an array to sort
        if(threads == 1){ //If  one thread don't do parallel
            serialSort(array, left, right,sorted);
        }else {
        int mid  = (right+left)/2; //Index for middle of array
        #pragma omp parallel sections 
        {
            #pragma omp section 
            {
                sort(array, left, mid,  threads/2,sorted); //split left half
            }
            #pragma omp section 
            {
                sort(array, mid+1, right, threads/2,sorted); //split right half
            }
        } 
        combine(array, left, mid, right,sorted); //combine
        }   
    }
}
/* 
*Recursively serial sorting of array 
*Inputs: unsorted array, left index, right index, dummy array to sort into
*/
void serialSort(double array[], int left, int right, double sorted[]){ 
    if(left<right){ //Check if we have an array to sort
        int mid  = (right+left)/2; //Index for middle of array
        serialSort(array, left, mid, sorted);
        serialSort(array, mid+1, right,sorted);
        combine(array,left,mid,right,sorted);
    }
}

/* 
*Sort and merge two arrays 
*Inputs: unsorted array, left index, middle index, right index, dummy array to sort into
*Returns: sorted array
*/
void combine(double array[], int leftIndex, int midIndex, int rightIndex, double sorted[]){
    int i;
    int left = leftIndex; //create indicies to increase
    int right = midIndex+1;
    for(i=left; i<=right;){ //Do from start to end of array
        if((left<=midIndex) && (right<=rightIndex)){ //Check that we're still within left and right bounds
            if(array[left]<array[right]){ //Check which element is smaller and place it at position i
                sorted[i]=array[left];
                left++;
                i++;
            }
            else {
                sorted[i]=array[right];
                right++;
                i++;
            }
        }        
        if(left>midIndex){ //If left index is higher than middle only right to merge
            for(;right<=rightIndex;right++,i++){   
                sorted[i]=array[right];
            }   
            break;
        }
        if(right>rightIndex){ //If right index higher than max then only left to merge
            for(;left<=midIndex;left++, i++){
                sorted[i]=array[left];
            }
            break;
        }
    }
    for(i=leftIndex;i<=rightIndex;i++){ //Copy array to original
        array[i]=sorted[i];
    }
}

/* 
*Generates pseudorandom double numbers between -500 and 500
*Returns: double between -500 and 500
*/
double randomNum(){
    int max = 500;
    int min = -500;
    double quotient = RAND_MAX/(max - min); 
    return min + (rand()/quotient);
}
