#include "tsp.h"

double minDist(int i, int* sol, const instance* inst) {
	int j, indMin;
	double dj;
	double minDist = inst->costs[sol[i] * inst->nnodes + sol[i + 1]];
	indMin = i + 1;
	for(j = i + 2; j < inst->nnodes; j++)
		if((dj = inst->costs[sol[i] * inst->nnodes + sol[j]]) < minDist){
			indMin = j;
			minDist = dj;
		} /*if*/
	swapInt(sol + i + 1, sol + indMin);
	return minDist;
} /*minDist*/


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
	compute_costs(inst);
	
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	for(i = 0; i < inst->nnodes; i++)
		sol[i] = i;
	for(i = 0; i < inst->nnodes - 1; i++)
        z += minDist(i, sol, inst);
	updateBest(z, sol, inst);
	
	write_out_file(inst, sol, "h_greedy");
    write_plotting_script("h_greedy");
    system("gnuplot gnuplot_out.p");
	
	free(sol);
	freeInst(inst);
	return 0;
}