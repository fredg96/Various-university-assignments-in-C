#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

typedef struct data_collection{
    double pos_x;
    double pos_y;
    double mass;
} data_t;

typedef struct vector  //(0,0) in top left corner
{
    double x;
    double y;
} vector;

typedef struct tree_node
{
    struct vector center;  //boundaries of subsquare
    struct vector centerOfMass;
    double width; //distance from center to edge of square
    double totalMass;
    struct tree_node *topLeftTree; //children trees
    struct tree_node *topRightTree; 
    struct tree_node *botLeftTree; 
    struct tree_node *botRightTree; 
    struct data_collection *object_data; //data for each interstallar object
    int isLeaf;
} node_t;

void fileRead(char * file, double **stored_data, int N);
void fileWrite(char * file, double **stored_data, int N);
void printTree(node_t *node);
void findSubSquare(node_t *parent, data_t *data, node_t **next_search, vector *nextCenter);
void initialiseChild(node_t **child, vector center, double width);
void split(node_t **node);
void add(node_t **node, data_t *data_to_add, vector center, double width);
void getAcceleration(node_t **node_trav, data_t *curr_object_i, vector *acc);
int inBound(double pos_x, double pos_y);
void deleatTree(node_t **node);

double G;
double theta_max;
double delta_t;

int main(int argc, char *argv[]){
    if(argc != 8){
        printf("Invalid number of inputs, terminating...\n");
        return 1;
    }
    const int N = atoi(argv[1]);
    char * input_file = argv[2];
    const int nsteps = atoi(argv[3]);
    delta_t = atof(argv[4]);
    theta_max = atof(argv[5]);
    const int graphics = atoi(argv[6]);
    int n_threads = atoi(argv[7]);
    if(graphics != 0){
        printf("Graphics not implemented, terminating.\n");
        return 1;
    } 

    G = (double)100/N;

    double *read_data = (double *)malloc(sizeof(double)*6*N); //buffer for inital condition
    data_t *stored_data = (data_t *)malloc(sizeof(data_t)*N); //node data
    vector *vel = (vector *)malloc(sizeof(vector)*N); //velocities
    double brightness[N];
    int i, j, k;
    j = 0;
    fileRead(input_file, &read_data, N);
    for(i = 0; i<N; i++){
        j = 6*i;
        stored_data[i].pos_x = read_data[j];
        stored_data[i].pos_y = read_data[j+1];
        stored_data[i].mass = read_data[j+2];
        vel[i].x = read_data[j+3];
        vel[i].y = read_data[j+4];
        brightness[i] = read_data[j+5];
    }
    
    free(read_data);

    node_t *root = NULL;
    vector center = {0.5, 0.5};
    double width = 0.5;
  
    //time step
    for(k=0; k<nsteps; k++){
        for(i=0; i<N; i++){ 
            add(&root, &stored_data[i], center, width);
        }
        #pragma omp parallel num_threads(n_threads)
        {
        #pragma omp for
        for(i=0; i<N; i++){  
            vector acc;
            acc.x = 0.0;
            acc.y = 0.0;
            getAcceleration(&root, &stored_data[i], &acc);
            vel[i].x += delta_t*acc.x;
            vel[i].y += delta_t*acc.y; 
        }

        #pragma omp for
        for(i=0;i<N;i++){              
            stored_data[i].pos_x += delta_t*vel[i].x;
            stored_data[i].pos_y += delta_t*vel[i].y;
            if(!inBound(stored_data[i].pos_x, stored_data[i].pos_y)){
                printf("Object out of bound. Terminating...\n");
                i = N; //quit current thread
                k = nsteps;
            }
        }
        }
        deleatTree(&root); 
    
    }
    double *write_data = (double *)malloc(sizeof(double)*6*N); //buffer for final timestep data
    j = 0;
    for(i = 0; i<N; i++){ //set final timestep data
        j = 6*i;
        write_data[j] = stored_data[i].pos_x;
        write_data[j+1] = stored_data[i].pos_y;
        write_data[j+2] = stored_data[i].mass;
        write_data[j+3] = vel[i].x ;
        write_data[j+4] = vel[i].y; 
        write_data[j+5] = brightness[i];
    }
    free(stored_data);
    free(vel);
    fileWrite("./result.gal", &write_data, N);
    free(write_data);
    return 0;
}

void fileRead(char * file, double **stored_data_ptr, int N){ 
    FILE *fp = fopen(file, "r");
    int n = fread((*stored_data_ptr), sizeof(double), N*6, fp);       
    if(n != N*6){
        printf("Error in reading input. \n");
    }
    fclose(fp); 
}

void fileWrite(char * file, double **stored_data_ptr, int N){ 
    FILE *fp = fopen(file, "w");
    fwrite((*stored_data_ptr), sizeof(double), N*6, fp);       
    fclose(fp); 
}

void printTree(node_t *node){
    if(node == NULL){
        printf("Tree is empty!\n");
        return;
    } 
    if(node != NULL){
        printf("Box: center: (%f, %f), width: %f, center of mass:(%f, %f), total mass: %f\n", node->center.x, node->center.y, node->width, node->centerOfMass.x, node->centerOfMass.y, node->totalMass);
        if(node->object_data != NULL){
            printf("Object: pos_x: %f, pos_y: %f, mass: %f\n", node->object_data->pos_x, node->object_data->pos_y, node->object_data->mass);
        }
    }
    if(node->topLeftTree != NULL){
        printTree(node->topLeftTree);
    }
    if(node->topRightTree != NULL){
        printTree(node->topRightTree);
    }
    if(node->botLeftTree!= NULL){
        printTree(node->botLeftTree);
    }
    if(node->botRightTree != NULL){
        printTree(node->botRightTree);
    }
}

// Input current node and the data to add. Finds the subnode where the data shold be put. 
void findSubSquare(node_t *parent, data_t *data, node_t **next_search, vector *nextCenter){
    double w = parent->width;
    double org_x = parent->center.x;
    double org_y = parent->center.y;
    double data_x = data->pos_x;
    double data_y = data->pos_y;
    if(org_x-w < data_x && data_x <= org_x && org_y-w < data_y && data_y < org_y){
        (*next_search) = parent->topLeftTree;
        nextCenter->x = org_x-w/2;
        nextCenter->y = org_y-w/2;
    } else if(org_x < data_x && data_x < org_x + w && org_y-w < data_y && data_y <= org_y){
        (*next_search) = parent->topRightTree;
        nextCenter->x = org_x+w/2; 
        nextCenter->y = org_y-w/2;
    } else if(org_x-w < data_x && data_x < org_x && org_y <= data_y && data_y < org_y+w){
        (*next_search) = parent->botLeftTree;
        nextCenter->x = org_x-w/2;
        nextCenter->y = org_y+w/2;
    } else if(org_x <= data_x && data_x < org_x+w && org_y < data_y && data_y < org_y+w){
        (*next_search) = parent->botRightTree;
        nextCenter->x = org_x+w/2; 
        nextCenter->y = org_y+w/2;
    } else{
        printf("Values: %f, %f\n", data_x, data_y);
        printf("Current box: center: (%f, %f), width: %f\n", org_x, org_y, w);
        printf("Unhadeled case. Terminating...\n");
        exit(0);
    }
}

void initialiseChild(node_t **child, vector center, double width){ 
    (*child)->object_data = NULL;
    (*child)->topLeftTree = NULL;
    (*child)->topRightTree = NULL;
    (*child)->botLeftTree = NULL;
    (*child)->botRightTree = NULL;
    (*child)->center = center;
    (*child)->centerOfMass.x = 0;
    (*child)->centerOfMass.y = 0;
    (*child)->width = width;
    (*child)->totalMass = 0;
    (*child)->isLeaf = 1;
}

void split(node_t **node){ 
    vector oldCenter = (*node)->center;
    double width = ((*node)->width)/2;
    (*node)->topLeftTree = (node_t *) malloc(sizeof(node_t));
    vector newCenter = {oldCenter.x-width, oldCenter.y-width};
    initialiseChild(&((*node)->topLeftTree), newCenter, width);
    (*node)->topRightTree = (node_t *) malloc(sizeof(node_t)); 
    newCenter.x = oldCenter.x+width; 
    initialiseChild(&((*node)->topRightTree), newCenter, width);
    (*node)->botLeftTree = (node_t *) malloc(sizeof(node_t)); 
    newCenter.x = oldCenter.x-width; 
    newCenter.y = oldCenter.y+width;
    initialiseChild(&((*node)->botLeftTree), newCenter, width);
    (*node)->botRightTree = (node_t *) malloc(sizeof(node_t));
    newCenter.x = oldCenter.x+width; 
    initialiseChild(&((*node)->botRightTree), newCenter, width);
    (*node)->isLeaf = 0;
}

void add(node_t **node, data_t *data_to_add, vector center, double width){
    if((*node) == NULL){ //leaf node, add here
        (*node) = (node_t *) malloc(sizeof(node_t));
        (*node)->object_data = data_to_add; 
        (*node)->topLeftTree = NULL;
        (*node)->topRightTree = NULL;
        (*node)->botRightTree = NULL;
        (*node)->botRightTree = NULL;
        (*node)->center = center;
        (*node)->centerOfMass.x = 0; 
        (*node)->centerOfMass.x = 0; 
        (*node)->totalMass = 0; 
        (*node)->width = width;
        (*node)->isLeaf = 1;
        return;
    } else if((*node)->object_data == NULL && !((*node)->isLeaf)){ //In the middle of the tree keep stepping
        (*node)->totalMass += data_to_add->mass;
        (*node)->centerOfMass.x += data_to_add->pos_x*data_to_add->mass;
        (*node)->centerOfMass.y += data_to_add->pos_y*data_to_add->mass;
        vector nextCenter;
        node_t *next_search;
        findSubSquare(*(node), data_to_add, &(next_search), &(nextCenter));
        add(&next_search, data_to_add, nextCenter, width/2);
    } else if((*node)->object_data == NULL && (*node)->isLeaf){ //Newly created node, add here
        (*node) -> object_data = data_to_add;
    } else{ //leaf node that we need to split
        if((*node)->object_data->pos_x == data_to_add->pos_x && (*node)->object_data->pos_y == data_to_add->pos_y){
            printf("Collision, objects at same location. Terminating...\n");
            exit(0);
        }
        split(node);
        data_t *old_data = (*node)->object_data;
        (*node)->object_data = NULL;
        (*node)->totalMass = old_data->mass + data_to_add->mass;
        (*node)->centerOfMass.x = old_data->pos_x*old_data->mass + data_to_add->pos_x*data_to_add->mass;
        (*node)->centerOfMass.y = old_data->pos_y*old_data->mass + data_to_add->pos_y*data_to_add->mass;
        vector nextCenter;
        node_t *next_search;
        findSubSquare(*(node), old_data, &(next_search), &(nextCenter));
        add(&next_search, old_data, nextCenter, width/2);
        findSubSquare(*(node), data_to_add, &(next_search), &(nextCenter));
        add(&next_search, data_to_add, nextCenter, width/2);
    }
}

void getAcceleration(node_t **node_trav, data_t *curr_object_i, vector *acc){
    double epsilon = 0.001;
    if((*node_trav) == NULL){
        return;
    } else if ((*node_trav)->isLeaf && (*node_trav)->object_data != NULL){ //leaf node where we have an object
        if(!((*node_trav)->object_data->pos_x == curr_object_i->pos_x && (*node_trav)->object_data->pos_y == curr_object_i->pos_y)){
            double dx = curr_object_i->pos_x-(*node_trav)->object_data->pos_x;
            double dy = curr_object_i->pos_y-(*node_trav)->object_data->pos_y;
            double dist = sqrt(dx*dx+dy*dy);
            (*acc).x += -G*((*node_trav)->object_data->mass)/pow(epsilon+dist,3)*dx;
            (*acc).y += -G*((*node_trav)->object_data->mass)/pow(epsilon+dist,3)*dy;
        }
    } else if(!(*node_trav)->isLeaf){ //middle of tree
        double dx =  (curr_object_i)->pos_x - (*node_trav)->center.x;
        double dy =  (curr_object_i)->pos_y - (*node_trav)->center.y;
        double dist = sqrt(dx*dx+dy*dy); //dist from object_i to center of region
        double theta = (((*node_trav)->width)*2)/dist; //calc local theta
        if(theta<=theta_max && (*node_trav)->totalMass > 0){ 
            dx = curr_object_i->pos_x - ((*node_trav)->centerOfMass.x)/(*node_trav)->totalMass;
            dy = curr_object_i->pos_y - ((*node_trav)->centerOfMass.y)/(*node_trav)->totalMass;
            dist = sqrt(dx*dx+dy*dy);
            (*acc).x += -G*((*node_trav)->totalMass)/pow(epsilon+dist,3)*dx;
            (*acc).y += -G*((*node_trav)->totalMass)/pow(epsilon+dist,3)*dy;
        } else{
            getAcceleration(&((*node_trav)->topLeftTree), curr_object_i, acc);
            getAcceleration(&((*node_trav)->topRightTree), curr_object_i, acc);
            getAcceleration(&((*node_trav)->botLeftTree), curr_object_i, acc);
            getAcceleration(&((*node_trav)->botRightTree), curr_object_i, acc);
        }
    } //last case for empty leaf nodes, do nothing
}

int inBound(double pos_x, double pos_y){
    return (pos_x < 1.0) && (pos_x > 0) && (pos_y < 1.0) && (pos_y > 0); //in original square
}

void deleatTree(node_t **node){
    if((*node) != NULL){
        if(!((*node)->isLeaf)){
            deleatTree(&((*node)->topLeftTree));
            deleatTree(&((*node)->topRightTree));
            deleatTree(&((*node)->botLeftTree));
            deleatTree(&((*node)->botRightTree));
        } 
        free((*node));
        (*node) = NULL;
    }
}