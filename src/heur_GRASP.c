#include "tsp.h"

double secMinDist(int i, int* sol, const instance* inst){
	int j, iMin, iSecMin;
	double dj;
	int last = sol[i];
	double minDist = inst->costs[last * inst->nnodes + sol[i + 1]];
	double secMinDist = inst->costs[last * inst->nnodes + sol[i + 2]];
	iMin = i + 1;
	iSecMin = i + 2;
	if(secMinDist < minDist){
		swapInt(&iMin, &iSecMin);
		swapDouble(&minDist, &secMinDist);
	} /*if*/
	for(j = i + 3; j < inst->nnodes; j++)
		if((dj = inst->costs[last * inst->nnodes + sol[j]]) < minDist){
			iSecMin = iMin;
			secMinDist = minDist;
			iMin = j;
			minDist = dj;
		} /*if*/
		else if(dj < secMinDist){
			iSecMin = j;
			secMinDist = dj;
		} /*if*/
	swapInt(sol + i + 1, sol + iSecMin);
	return secMinDist;
} /*secMinDist*/

double chooseNode(int i, int nearest_prob, int* sol, const instance* inst){
	if(rand() % 100 < nearest_prob)
		return minDist(i, sol, inst);
	else 
		return secMinDist(i, sol, inst);
} /*chooseNode*/

double testN(instance* inst, int* sol, int prob) {
	int i, j;
	double z_avg = 0;
	double z;
	for(i = 0; i < inst->n_sim; i++){
		srand(i);
		z = 0;
		for(j = 0; j < inst->nnodes; j++)
			sol[j] = j;
		for(j = 0; j < inst->nnodes - 2; j++)
			z += chooseNode(j, prob, sol, inst);
		z += inst->costs[sol[j] * inst->nnodes + sol[j + 1]] + inst->costs[sol[j + 1]]; /*add last two edges costs, which don't need computation*/
		updateBest(z, sol, inst);
		z_avg += z;
	} /*for*/
	return z_avg / inst->n_sim;
} /*testN*/

int main(int argc, char **argv)
{
	int* sol;
	int prob;
	double z_avg_best, z_avg_curr;
	FILE* stats;
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	
	initInst(inst);
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
	compute_costs(inst);
	
	stats = fopen("../out/stats_GRASP.dat", "w");
	if (stats == NULL) {
		printf("Couldn't create the stats file!");
		return FILE_OPEN_ERR;
	} /*if*/
	fprintf(stats, "#PROBABILITY and AVG COST\n");
	
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	inst->best_prob = 0;
	z_avg_best = testN(inst, sol, 0);
	fprintf(stats, "%d  %.2f\n", 0, z_avg_best);
	for(prob = 1; prob <= 100; prob++) {
		if((z_avg_curr = testN(inst, sol, prob)) < z_avg_best){
			z_avg_best = z_avg_curr;
			inst->best_prob = prob;
		} /*if*/
		fprintf(stats, "%d  %.2f%s", prob, z_avg_curr, ((prob == 100) ? "" : "\n"));
	} /*for*/
	
	printf("Average best cost %f obtained with prob. of choosing nn equals to %d\n", z_avg_best, inst->best_prob);
	write_out_file(inst, inst->best_sol, "h_GRASP");
    write_plotting_script("h_GRASP");
    system("gnuplot gnuplot_out.p");
	system("gnuplot gnuplot_stats.p");
	
	free(sol);
	freeInst(inst);
	return 0;
}