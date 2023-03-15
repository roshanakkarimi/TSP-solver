#include "tsp.h"

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
	z += inst->costs[sol[i]];
	updateBest(z, sol, inst);
	
	write_out_file(inst, inst->best_sol, "h_greedy");
    write_plotting_script("h_greedy");
    system("gnuplot gnuplot_out.p");
	/*printf("z = %f\n", z);
	for(i = 0; i < inst->nnodes; i++)
		printf("%d, %s", sol[i], (((i+1)%4==0) ? "\n" : ""));*/
	
	free(sol);
	freeInst(inst);
	return 0;
}