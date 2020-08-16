#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct vector
{
   double x;
   double y;
} vector;

void fileRead(char * file, double **stored_data_ptr, int N);
void fileWrite(char * file, double **stored_data_ptr, int N);
void step(vector **pos, double *mass, vector **vel, int N, double epsilon, double G, double delta_t);

int main(int argc, char *argv[]){
    if(argc != 6){
        printf("Wrong number of input igruments.\nExpexted N, filename, nsteps, delta_t graphics.");
        return 1; 
    }
    int N = atoi(argv[1]);
    char * input_file = argv[2];
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics = atoi(argv[5]);
    if(graphics != 0 && graphics != 1){
        printf("Invalid input graphics, terminating.\n");
        return 1;
    } 
    double epsilon = 1E-3;
    double G = (double)100/N;
    double *read_data = (double *)malloc(sizeof(double)*6*N); //buffer for inital condition
    fileRead(input_file, &read_data, N);
    double brightness[N];
    double *mass = (double *)malloc(sizeof(double)*N);
    vector *pos = (vector *)malloc(sizeof(vector)*N);
    vector *vel = (vector *)malloc(sizeof(vector)*N);
    int i, j;
    for(i=0; i<N; i++){
        j = 6*i;
        pos[i].x = read_data[j];
        pos[i].y = read_data[j+1];
        mass[i] = read_data[j+2];
        vel[i].x = read_data[j+3];
        vel[i].y = read_data[j+4];
        brightness[i] = read_data[j+5];
    }

    free(read_data);
    
    for(i=0; i<nsteps; i++){
        step(&pos, mass, &vel, N, epsilon, G, delta_t);
    }

    double *write_data = (double *)malloc(sizeof(double)*6*N); //buffer for final timestep data
    j = 0;
    for(i = 0; i<N; i++){ //set final timestep data
        j = 6*i;
        write_data[j] = pos[i].x;
        write_data[j+1] = pos[i].y;
        write_data[j+2] = mass[i];
        write_data[j+3] = vel[i].x ;
        write_data[j+4] = vel[i].y; 
        write_data[j+5] = brightness[i];
    }
    free(pos);
    free(mass);
    free(vel);
    fileWrite("./result.gal", &write_data, N);
    free(write_data);
    return 0;
}

void fileRead(char * file, double **stored_data_ptr, int N){ 
    FILE *fp = fopen(file, "r");
    int n = fread((*stored_data_ptr), sizeof(double), N*6, fp);    
    if(n != N*6){
        printf("Error in reading input.\n");
    }   
    fclose(fp); 
}

void fileWrite(char * file, double **stored_data_ptr, int N){ 
    FILE *fp = fopen(file, "w");
    fwrite((*stored_data_ptr), sizeof(double), N*6, fp);       
    fclose(fp); 
}

inline void step(vector **pos, double *mass, vector **vel, int N, double epsilon, double G, double delta_t){
    int i, j;
    double dx, dy, r, tempx, tempy; 
    for(i=0; i<N; i++){
        tempx = 0;
        tempy = 0;
        for(j=0; j<N; j++){
            if(i != j){
                dx = (*pos)[i].x-(*pos)[j].x;
                dy = (*pos)[i].y-(*pos)[j].y;
                r = sqrt(dx*dx+dy*dy);
                tempx += (-G*mass[j]*dx)/pow(r+epsilon, 3);
                tempy += (-G*mass[j]*dy)/pow(r+epsilon, 3);
            }
        }
        (*vel)[i].x = (*vel)[i].x + tempx*delta_t;
        (*vel)[i].y = (*vel)[i].y + tempy*delta_t;
    }
    for(i=N-1; i>=0; i--){
        (*pos)[i].x = (*pos)[i].x + (*vel)[i].x*delta_t;
        (*pos)[i].y = (*pos)[i].y + (*vel)[i].y*delta_t;
    }
} 