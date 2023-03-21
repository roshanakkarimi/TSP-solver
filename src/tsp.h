#ifndef TSP_H
#define TSP_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>

/*constants*/

#define INFINITE LONG_MAX
#define EPSILON 0.00001
#define DEFAULT_RAND 1

/*errors*/

#define FILE_OPEN_ERR 1
#define FILE_STRUCT_ERR 2
#define FORMAT_ERR 3
#define INDEX_ERR 4

/*heur. modes*/

#define GREEDY 0
#define GRASP 1

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
	int prob; /*probability*/
	int mode;
	char verbosity;
	char fileIn[100];
	bool two_opt;
	
/*results*/
	double zbest; /*best value for the obj. function*/
	int tbest; /*time for the best solution*/
	int* best_sol; 
	int best_prob;
} instance;	

/*generic functions*/

typedef double (*cost)(int ind1, int ind2, const point* pts);
typedef double (*node_picker)(int i, int nearest_prob, int* sol, const instance* inst);

/*generic utilities*/

void swapInt(int*, int*);
void swapDouble(double*, double*);

/*cost utilities*/

double dist(int, int, const point*);
double sq_dist(int, int, const point*);
double minCost(int, int*, const instance*);

/*input elaboration and initialization*/

void initInst(instance*);
bool parse_cmd(int, char**, instance*);
int read_fileIn(instance*);
void compute_costs(instance*, cost);

/*managing errors and debug*/

bool checkSol(double, const int*, const instance*);

/*managing solutions*/

void updateBest(double, const int*, instance*);

/*output elaboration*/

int write_out_file(const instance*, const int*, const char*, const char*);
int write_plotting_script(const char*, int, int);

/*memory management*/

void freeInst(instance*);
 
#endif