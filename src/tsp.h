#ifndef TSP_H
#define TSP_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

/*constants*/

#define INFINITE 1e30
#define EPSILON 0.00001
#define DEFAULT_RAND 1

/*errors*/

#define FILE_OPEN_ERR 1
#define FILE_STRUCT_ERR 2
#define FORMAT_ERR 3
#define INDEX_ERR 4

/*data structures*/

typedef struct {
	double x;
	double y;
} point;

typedef struct {
/*input data*/
	int nnodes; 	
	point* pts;
	
/*costs*/
	double* costs;
	
/*parameters*/
	double timelimit; /*total time limit*/
	int randseed;
	int n_sim;
	char verbosity;
	char fileIn[100];
	
/*results*/
	double zbest; /*best value for the obj. function*/
	double tstart;
	double tbest; /*time for the best solution*/
	int* best_sol; 
	int best_prob;
} instance;	


/*generic utilities*/

void swapInt(int*, int*);
void swapDouble(double*, double*);

/*distance utilities*/

double dist(int, int, const point*);
double minDist(int, int*, const instance*);

/*input elaboration and initialization*/

void initInst(instance*);
void parse_cmd(int, char**, instance*);
int read_fileIn(instance*);
void compute_costs(instance*);

/*managing errors and debug*/

bool checkSol(double, const int*, const instance*);

/*managing solutions*/

void updateBest(double, const int*, instance*);

/*output elaboration*/

int write_out_file(const instance*, const int*, const char*);
int write_plotting_script(const char*);

/*memory management*/

void freeInst(instance*);
 
#endif