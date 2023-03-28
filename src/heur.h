#ifndef HEUR_H
#define HEUR_H

#include "tsp.h"

/*global*/

double (*kickers[1])(const instance*, int*, int*, int);
const char mods[2][50];

/*heuristic node pickers*/

typedef double (*node_picker)(int i, int nearest_prob, int* sol, const instance* inst);
double (*pickers[2])(int, int, int*, const instance*);
double greedy_picker(int, int, int*, const instance*); /*just a wrapper to use node_picker standard*/
double grasp_picker(int, int, int*, const instance*);

/*solver*/

void gr_solve(instance*, int*, int, node_picker);

/*refinement utilities*/

void reverse_sub(int*, int, int, int);
double impr(double*, int*, int, int, int);

/*kickers from local min.*/

double tabu_kicker(const instance*, int*, int*, int);

/*tabu utilities*/

int choose_tabu(int*, int, int, int);

/*refinement*/

double two_opt_move(const instance*, int*, int*, int, bool);
void refine(instance*, int*);

/*output files and scripts management*/

int open_init_stats(FILE**, const instance*);
int open_out_files(instance*, FILE**, FILE**);



#endif