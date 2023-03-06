#ifndef TSP_H
#define TSP_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*constants*/

#define INFINITE 1e15
#define DEFAULT_RAND 1

/*data structures*/

typedef struct {
	double x;
	double y;
} point;

typedef struct {
/*input data*/
	int nnodes; 	
	double* demand;   
	point* pts;
	
/*parameters*/
	double timelimit; /*total time limit*/
	int randseed;
	char verbosity;
	char fileIn[100];
	
/*results*/
	double zbest; /*best value for the obj. function*/
	double tstart;
	double tbest; /*time for the best solution*/
	int* best_sol; 
} instance;	


/*utilities*/

double ptdist(point*, point*);
double dist(int, int, instance*);


/*input elaboration*/

void parse_cmd(int, char**, instance*);
 
#endif