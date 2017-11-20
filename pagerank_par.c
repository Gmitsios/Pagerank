/*
Parallel and Distributed Systems
--Pagerank PThread Algorithm--
-Author: Mitsios Georgios
-September 2014
*/

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>


#define checkResults(string, val) { \
    if (val){                       \
        printf("Failed with %d at %s\n", val, string);  \
        exit(1);\
    }\
}\

//General Variables accesed by all functions
int nodes, edges, bounce, **G , *c, **L;
int num_threads, chk=0, count=0;
double *x, *last_x;

//Declaration PThread attributes 
pthread_t *thread;
pthread_attr_t joinable;
pthread_barrier_t barrier;
pthread_mutex_t mutex;

void import_data(char *file){

	int tmp, i;
	FILE *fp;

	fp = fopen(file, "r+");
	if (fp == NULL){
   		printf("Error Opening Input Data File\n");
    	exit(1);
	}


	for(i=0; i<edges; i++){
		tmp = fscanf(fp, "%d %d", &G[0][i], &G[1][i]);
		if( tmp != 2 ){
			printf("Error Reading Input Data File!");
			exit(2);
		}
	}

	fclose(fp);
}

void save_data(char *file){

	int i;

	FILE *fp;
   	fp = fopen (file, "w+");

   	for(i=0; i<nodes; i++){
   		fprintf(fp, "%f\n", x[i]);
   	}

   	printf("Data stored succesfully\n");

   	fclose(fp);
}

//This function returns the error of x[i] - last_x[i] used to terminate the loop in the power method function
double error_fun(){

	int i;
	double tmp=0, error;

	error = fabs( x[0] - last_x[0] );

	for(i=1;i<nodes;i++){
		tmp = fabs( x[i]-last_x[i] );
		if (tmp > error){ error = tmp; }
	}

	return (error);
}

//This function computes the c[i] of each node
void *max_references(void *arg){

	int i, j, tid = (int)(intptr_t)arg;

	for (i=tid; i<nodes; i+=num_threads){

		c[i] = 0;
		for (j=0; j<edges; j++){
			if (G[1][j] == i){
				c[i]++;
			}
		}
	}
}

//This function computes the L matrix
void *link_struct_fun(void *arg){

	int i, j, tid = (int)(intptr_t)arg;

	if(tid == 0){
		printf("Generating Arrays of Links per Page...\n");
	}


	for (i=tid*bounce; i<(tid*bounce+bounce); i++){

		c[i] = 0;
		for (j=0; j<edges; j++){
			if (G[1][j] == i){
				L[i][c[i]] = G[0][j];
				c[i]++;

			}
		}
	}
}

//Same as the above funtion for the remaining nodes (nodes%num_threads)
void *link_struct_fun_rem(void *arg){

	int i, j, tid = (int)(intptr_t)arg;

		c[bounce*num_threads+tid] = 0;
		for (j=0; j<edges; j++){
			if (G[1][j] == bounce*num_threads+tid){
				L[i][c[bounce*num_threads+tid]] = G[0][j];
				c[bounce*num_threads+tid]++;
			}
		}
}

//This is the function that computes the possibility of being directed to each page
//using the power method
void *power_method_fun(void *arg){

	int pos, i, j;
	int tid = (int)(intptr_t)arg;
	double err;
	double delta, p=0.85;

	if(tid == 0){
		printf("Calculating Possibility of each Page...\n");
	}

	delta = (1.0-p)/(double)nodes;

	for(i=tid;i<nodes;i+=num_threads){
		x[i] = 1.0/(double)nodes;
		last_x[i] = 0.0;
	}
	pthread_barrier_wait(&barrier);

	do{
		for(i=tid;i<nodes;i+=num_threads){
			last_x[i] = x[i];
			x[i] = 0.0;

		}
		pthread_barrier_wait(&barrier);

		for(i=tid;i<nodes;i+=num_threads){

			if ( c[i] == 0 ){
				chk = pthread_mutex_lock(&mutex);
				checkResults("mutex lock",chk);
				for(j=0;j<nodes;j++){
					x[j] = x[j] + last_x[i]/(double)nodes;
				}
				chk = pthread_mutex_unlock(&mutex);
				checkResults("mutex unlock",chk);
			}
			else {
				chk = pthread_mutex_lock(&mutex);
				checkResults("mutex lock",chk);
				for(j=0;j<c[i];j++){
					pos = L[i][j];
					x[pos] = x[pos] + last_x[i]/(double)c[i];
				}
				chk = pthread_mutex_unlock(&mutex);
				checkResults("mutex unlock",chk);
			}


		}

		pthread_barrier_wait(&barrier);

		for(i=tid;i<nodes;i+=num_threads){
			x[i] = p*x[i] + delta;
		}
		pthread_barrier_wait(&barrier);


		if(tid == 0){
			count++;
		}

		err = error_fun();

	}while(err > 0.0000001);

	if(tid == 0){
		printf("Finished after %d iterations.\n", count);
	}

}


void main(int argc, char **argv){

	struct timeval first, second, lapsed;
  	struct timezone tzp;

	int i, j, cmax=1;
	double grp;
	char *input_file;
	char *output_file;

	input_file = argv[1];
	nodes = atoi(argv[2]);
	edges = atoi(argv[3]);
	output_file = "results.txt";

	num_threads = atoi(argv[4]);

	printf("NUM THREADS %d \n", num_threads);
	printf("Dataset File: %s \n", input_file);
	printf("Number of Nodes: %d\n", nodes);
	printf("Number of Edges: %d\n\n", edges);

	//Initialization of PThread attributes 
	pthread_barrier_init(&barrier, NULL, num_threads);
	pthread_mutex_init(&mutex, NULL);
  	thread = (pthread_t *)malloc(num_threads*sizeof(pthread_t));


	grp = ((double)nodes)/((double)num_threads);
	printf("Group size: %f \n", grp);
	bounce = (int)floor(grp);

	printf("Bounce step: %d \n\n", bounce);

	//Memory Allocation
	G = (int **)malloc(2*sizeof(int *));
	for (i=0; i<2; i++){
		G[i] = (int *)malloc(edges*sizeof(int));
		for (j=0; j<edges; j++){
			G[i][j] = 0;
		}
	}
	c = (int *)malloc(nodes*sizeof(int));

	x = (double *)malloc(nodes*sizeof(double));
	last_x = (double *)malloc(nodes*sizeof(double));

	//Import Data
	import_data(input_file);

	printf("Initializing Data... \n");



	//Call the "max_references" function to find number of links per page
	for(i=0;i<num_threads;i++){
		chk = pthread_create( &thread[i], NULL, max_references, (void *)(intptr_t)i );
		checkResults("finding max amount of references",chk);
	}

	for(i=0; i<num_threads; i++) {
		chk = pthread_join(thread[i],NULL);
		checkResults("pthread join no1",chk);
	}

	//Allocating memory for L since we now got its number of columns per line c[i]
	L = (int **)malloc(nodes*sizeof(int *));
	for(i=0; i<nodes; i++){
		L[i] = (int *)malloc(c[i]*sizeof(int));
		for (j=0; j<c[i]; j++){
			L[i][j] = -1;
		}
	}


	printf("\n --Start of Parallel Section-- \n");


	gettimeofday(&first, &tzp);

	//Parallel call of the "link_struct_fun" function, "num_thread" times 
	for(i=0;i<num_threads;i++){
		chk = pthread_create( &thread[i], NULL, link_struct_fun, (void *)(intptr_t)i );
		checkResults("link structure initialization",chk);
	}

	for(i=0; i<num_threads; i++) {
		chk = pthread_join(thread[i],NULL);
		checkResults("pthread join no3",chk);
	}

	//Parallel call of the "link_struct_fun" function, "nodes(mod)num_thread" times 
	for(i=0;i<(nodes%num_threads);i++){
		chk = pthread_create( &thread[i], &joinable, link_struct_fun_rem, (void *)(intptr_t)i );
		checkResults("link structure initialization",chk);
	}

	for(i=0; i<(nodes%num_threads); i++) {
		chk = pthread_join(thread[i],NULL);
		checkResults("pthread join no4",chk);
	}

	//Parallel call of the "power_method_fun" function, "num_thread" times 
	for(i=0;i<num_threads;i++){
		chk = pthread_create( &thread[i], NULL, power_method_fun, (void *)(intptr_t)i );
		checkResults("power method initialization",chk);
	}

	for(i=0; i<num_threads; i++) {
		chk = pthread_join(thread[i],NULL);
		checkResults("pthread join no5",chk);
	}

	gettimeofday(&second, &tzp);

	printf("--End of Parallel Section-- \n\n");

	lapsed.tv_usec = second.tv_usec - first.tv_usec;
  	lapsed.tv_sec = second.tv_sec - first.tv_sec;
 	printf("Time elapsed: %d.%d s\n", (int)(lapsed.tv_sec), (int)(lapsed.tv_usec));

	pthread_mutex_destroy(&mutex);
	pthread_barrier_destroy(&barrier);

	//Save Data
	save_data(output_file);


}
