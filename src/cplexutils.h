#include "tsp.h"
#include <cplex.h>

void initCplex(instance *inst, CPXENVptr env);
int xpos(int i, int j, instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void print_error(const char *err);
