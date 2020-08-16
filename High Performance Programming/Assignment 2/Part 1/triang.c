#include <stdio.h>
#include <stdlib.h>



// The main program.
int main( int argc, char *argv[] )
{
    if(argc != 2) {
        printf("Incorrect number of input arguemnts. \n"); 
        return 1; 
    } 
    int num = atof(argv[1]);
    int i, j, val;
    for(i = 0; i < num; i++) {
        for(j = 0; j <= i; j++) {
            if(j == 0) {
                val = 1; 
                printf("%d\t", val);
            } else{
                val = val * (i+1-j)/j;
                printf("%d\t", val);
            }
        }
        printf("\n");
    }

    return 0; 
}
