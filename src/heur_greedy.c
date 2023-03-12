#include "tsp.h"
#include <stdbool.h>

void swapInt(int* i1, int* i2) {
	int temp;
	temp = *i2;
	*i2 = *i1;
	*i1 = temp;
} /*swapInt*/

double minDist(int i, int* sol, const instance* inst) {
	int j, indMin;
	double dj;
	double minDist = inst->costs[i * inst->nnodes + i + 1];
	indMin = i + 1;
	for(j = i + 2; j < inst->nnodes; j++)
		if((dj = inst->costs[i * inst->nnodes + j]) < minDist){
			indMin = j;
			minDist = dj;
		} /*if*/
	swapInt(sol + i, sol + indMin);
	return minDist;
} /*minDist*/

bool checkSol(double z, const int* sol, const instance* inst){
	int i;
	double z_check = 0;
	int* count = calloc(inst->nnodes, sizeof(int));
	assert(count != NULL);
	for(i = 0; i < inst->nnodes; i++)
		count[sol[i]]++;
	for(i = 0; i < inst->nnodes; i++)
		if(count[i] != 1){
			free(count);
			return false;
		} /*if*/
	free(count);
	for(i = 0; i < inst->nnodes; i++)
		z_check += inst->costs[i * inst->nnodes + sol[i]];
	if(abs(z_check - z) > EPSILON)
		return false;
	return true;
} /*checkSol*/

void updateBest(int z, const int* sol, instance* inst){
	int i;
	if(checkSol(z, sol, inst) && (z < inst->zbest || inst->zbest == -1)) { /*lazy eval.*/
		inst->zbest = z;
		for(i = 0; i < inst->nnodes; i++)
			inst->best_sol[i] = sol[i];
	} /*if*/
} /*updateBest*/

int main(int argc, char **argv)
{
	int* sol;
	int i;
	double z = 0;
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	initInst(inst);
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	for(i = 0; i < inst->nnodes; i++) 
		sol[i] = i;
	for(i = 0; i < inst->nnodes - 1; i++)
        z += minDist(i, sol, inst);
	updateBest(z, sol, inst);
    /*write_plotting_script(inst);
    system("gnuplot gnuplot_out.p");*/
	free(sol);
	free(inst);
	return 0;
}