#include <stdio.h>
#include <stdlib.h>

typedef struct list_node{
    int index;
    double min; 
    double max;
    struct list_node *next;
} node_t; 

void print_database(node_t *first);
void add(int index, double min, double max, node_t **first);
void free_list(node_t **first);
void deleat_index(int index, node_t **first);

int main(){
    node_t *first = NULL;
    char op = 1; 
    while(1){
        
        printf("Enter command:");
        fflush(stdin);
        int count = scanf(" %c", &op); //return number of read characters
        if(count != 1){ 
            printf("Invalid number of inputs, try again.\n");
            continue;
        }
        if (op == 'A'){ //Add entry to database
            int index;
            double min;
            double max; 
            scanf(" %d", &index);
            scanf(" %le", &min);
            scanf(" %le", &max);
            fflush(stdin);
            if(index > 31 || index < 1){
                printf("Date to add is not in the range of january, try again.\n");
                continue;
            }
            if(min > max){
                printf("Minimum temerature is higher than maximum temerature, try again.\n");
                continue;
            }
            add(index, min, max, &first);
        } else if(op == 'D'){ //Deleat index
            int index;
            scanf(" %d", &index);
            fflush(stdin);
            if(index > 31 || index < 1){
                printf("Date to deleat is not in the range of january, try again.\n");
                continue;
            }
            deleat_index(index, &first);
        } else if(op == 'P') { //Print
            print_database(first);
        } else if(op == 'Q'){ //Quit
            free_list(&first);
            break;
        } else{
            printf("Invalid operator.\n");
        }
    }
    return 0;
}

void print_database(node_t *first){
    if(first == NULL) {
        printf("Empty database.\n");
        return;
    }
    printf("day\tmin\t\tmax\n");
    while(first != NULL){
        printf("%d\t%f\t%f\n", first->index, first->min, first->max);
        first = first -> next;
    }
}

void add(int index, double min, double max, node_t **first){
    if((*first) == NULL){ //add to empty list or last in list
        (*first) = (node_t *) malloc(sizeof(node_t));
        (*first)->index = index;
        (*first)->min = min;
        (*first)->max = max; 
        (*first)->next = NULL;
    } else if((*first)->index == index){ //if index already exist
        (*first)->min = min; 
        (*first)->max = max; 
    } else if((*first) != NULL && index < (*first)->index){ //push
        node_t *new_node = (node_t *) malloc(sizeof(node_t));
        new_node->index = index;
        new_node->min = min; 
        new_node->max = max;
        new_node->next = (*first);
        (*first) = new_node;
    } else if((*first)->next != NULL && (*first)->next->index > index){ //add in the middle of the list
        node_t *new_node = (node_t *) malloc(sizeof(node_t));
        new_node->index = index;
        new_node->min = min; 
        new_node->max = max;
        new_node->next = (*first)->next;
        (*first)->next = new_node;
    } else{ //step
        add(index, min, max, &((*first)->next));
    }
}

void free_list(node_t **first){ //belong to quit command
    if((*first) != NULL){
        if((*first)->next != NULL){
            free_list(&((*first)->next));
        }
        free((*first));
        (*first) = NULL;
    }
}

void deleat_index(int index, node_t **first){
    if((*first) == NULL || (*first)->index > index){ //element not found
        printf("Day to deleat is not in the database.\n");
        return;
    } else if((*first)->index == index){ //remove first element
        node_t *temp = (*first);
        (*first) = (*first)->next;
        free(temp);
    } else if((*first)->next != NULL && (*first)->next->index == index){ //remove later
        node_t *temp = (*first)->next;
        (*first)->next = (*first)->next->next;
        free(temp);
    } else{
        deleat_index(index, &((*first)->next)); //step
    }
} 

