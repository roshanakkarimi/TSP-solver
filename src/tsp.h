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
#include <float.h>

/*constants*/

#define INFINITE_DBL DBL_MAX
#define EPSILON 0.00001
#define DEFAULT_RAND 1
#define DEFAULT_TT 1
#define MIN_NODES 75
#define MAX_NODES 150
#define MAX_COORD 7500

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
	int mode;
	char verbosity;
	char fileIn[100];
	
/*options*/
	bool two_opt;
	bool tabu;
	
/*specific pars.*/
	int n_sim; /*for test mode*/
	int prob; /*for grasp policy functions*/
	int tabu_tenure; /*for tabu alg.*/
	
/*running pars.*/
	double t_start;
	
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
double secMinCost(int, int*, const instance*);
double getSolCost(const int*, const instance*);

/*input elaboration and initialization*/

void initInst(instance*);
bool parse_cmd(int, char**, instance*);
void alloc_inst(instance*);
int read_fileIn(instance*);
void rand_points(instance*);
void compute_costs(instance*, cost);

/*managing errors and debug*/

int myError(const char*, int);
bool checkSol(double, const int*, const instance*);

/*solving*/

void updateBest(double, const int*, instance*);
double greedy_picker(int, int, int*, const instance*); /*just a wrapper to use node_picker standard*/
double grasp_picker(int, int, int*, const instance*);

/*output elaboration*/

int write_out_file(const instance*, const int*, const char*, const char*);
int write_plotting_script(const char*, int, int);

/*memory management*/

void freeInst(instance*);

/*global*/

double (*pickers[2])(int, int, int*, const instance*);
const char mods[2][50];
 
#endif