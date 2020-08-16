#include <stdio.h>
#include <stdlib.h>

int main(){
    FILE *fp = fopen("little_bin_file", "r"); // open file for reading
    int i; 
	double dou; 
	char ch; 
	float floa; 
	fread(&i, sizeof(int), 1, fp);
	fread(&dou, sizeof(double), 1, fp);
	fread(&ch, sizeof(char), 1, fp);
	fread(&floa, sizeof(float), 1, fp);
	fclose(fp);
	printf("%d\n%f\n%c\n%f\n", i, dou, ch, floa);
    return 0;
}

