/*
Parallel and Distributed Systems
--Pagerank Sequential Algorithm--
-Author: Mitsios Georgios
-September 2014
*/

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

void import_data(char *file, int N, int **E){

	int tmp, i, j;
	FILE *fp;

	fp = fopen(file, "r+");
	if (fp == NULL){
   		printf("Error Opening Input Data File\n");
    	exit(1);
	}


	for(i=0; i<N; i++){
		tmp = fscanf(fp, "%d %d", &E[0][i], &E[1][i]);
		if( tmp != 2 ){
			printf("Error Reading Input Data File!");
			exit(2);
		}
	}

	fclose(fp);
}

void save_data(double *xi, int nd, char *file){

	int i;

	FILE *fp;
   	fp = fopen (file, "w+");

   	for(i=0; i<nd; i++){
   		fprintf(fp, "%f\n", xi[i]);
   	}

   	printf("Data stored succesfully\n");

   	fclose(fp);
}

//This function returns the error of x[i] - last_x[i] used to terminate the loop in the power method function
double error_fun(int nd2, double *x2, double *last_x2){

	int i;
	double tmp=0, error;

	error = fabs( x2[0] - last_x2[0] );

	for(i=1;i<nd2;i++){
		tmp = fabs( x2[i]-last_x2[i] );
		if (tmp > error){ error = tmp; }
	}

	return error;
}

//This function computes the L matrix
void link_struct_fun(int nd, int ed, int *c, int **L, int **E){

	int i,j;

	printf("Generating Arrays of Links per Page...\n");

	for (i=0; i<nd; i++){

		c[i] = 0;
		for (j=0; j<ed; j++){
			if (E[1][j] == i){
				c[i]++;
				L[i] = (int *) realloc(L[i], c[i]*sizeof(int));
				L[i][c[i]-1] = E[0][j];
			}
		}
	}
}

//This is the function that computes the possibility of being directed to each page
//using the power method
void power_method_fun(int nd, int *c1, int **L1, double *x1, double *last_x1){

	int pos, i, j, count=0;
	double err;
	double delta, p=0.85;

	printf("Calculating Possibility of each Page...\n");

	delta = (1-p)/(double)nd;
	count = 0;

	for(i=0;i<nd;i++){
		x1[i] = 1/(double)nd;
		last_x1[i] = 0;
	}

	do{
		for(i=0;i<nd;i++){
			last_x1[i] = x1[i];
			x1[i] = 0;

		}

		for(i=0;i<nd;i++){
			if ( c1[i] == 0 ){
				for(j=0;j<nd;j++){
					x1[j] = x1[j] + last_x1[i]/(double)nd;
				}
			}
			else {
				for(j=0;j<c1[i];j++){
					pos = L1[i][j];
					x1[pos] = x1[pos] + last_x1[i]/(double)c1[i];
				}
			}
		}

		for(i=0;i<nd;i++){
			x1[i] = p*x1[i] + delta;
		}

		count++;
		err = error_fun(nd, x1, last_x1);
	}while(err > 0.0000001);

	printf("Finished after %d iterations.\n", count);

}


void main(int argc, char **argv){

	struct timeval first, second, lapsed;
  	struct timezone tzp;

	int nodes, edges, count, **G , *c, **L;
	int i, j;
	double *x, *last_x;
	char *input_file;
	char *output_file;

	input_file = argv[1];
	nodes = atoi(argv[2]);
	edges = atoi(argv[3]);
	output_file = "results.txt";

	printf("Dataset File: %s \n", input_file);
	printf("Number of Nodes: %d\nNumber of Edges: %d\n\n", nodes, edges );


	//Memory Allocation
	G = (int **)malloc(2*sizeof(int *));
	for (i=0; i<2; i++){
		G[i] = (int *)malloc(edges*sizeof(int));
		for (j=0; j<edges; j++){
			G[i][j] = 0;
		}
	}
	c = (int *)malloc(nodes*sizeof(int));
	L = (int **)malloc(nodes*sizeof(int *));

	import_data(input_file, edges, G);

	x = (double *)malloc(nodes*sizeof(double));
	last_x = (double *)malloc(nodes*sizeof(double));

	gettimeofday(&first, &tzp);

	//Call the main functions to determine pageranks
	link_struct_fun(nodes, edges, c, L, G);
	power_method_fun(nodes, c, L, x, last_x);

	gettimeofday(&second, &tzp);

	lapsed.tv_usec = second.tv_usec - first.tv_usec;
  	lapsed.tv_sec = second.tv_sec - first.tv_sec;
 	printf("Time elapsed: %d.%d s\n", (int)(lapsed.tv_sec), (int)(lapsed.tv_usec));

	save_data(x, nodes, output_file);

	free(G);
	free(L);
	free(c);
	free(x);
	free(last_x);

}
