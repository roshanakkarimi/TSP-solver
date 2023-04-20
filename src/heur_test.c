#include "tsp.h"

double secMinCost(int i, int* sol, const instance* inst){
	int j, iMin, iSecMin;
	double dj;
	int last = sol[i];
	double minCost = inst->costs[last * inst->nnodes + sol[i + 1]];
	double secMinCost = inst->costs[last * inst->nnodes + sol[i + 2]];
	iMin = i + 1;
	iSecMin = i + 2;
	if(secMinCost < minCost){
		swapInt(&iMin, &iSecMin);
		swapDouble(&minCost, &secMinCost);
	} /*if*/
	for(j = i + 3; j < inst->nnodes; j++)
		if((dj = inst->costs[last * inst->nnodes + sol[j]]) < minCost){
			iSecMin = iMin;
			secMinCost = minCost;
			iMin = j;
			minCost = dj;
		} /*if*/
		else if(dj < secMinCost){
			iSecMin = j;
			secMinCost = dj;
		} /*if*/
	swapInt(sol + i + 1, sol + iSecMin);
	return secMinCost;
} /*secMinCost*/

double greedy_picker(int i, int nearest_prob, int* sol, const instance* inst){ /*just a wrapper to use node_picker standard*/
	return minCost(i, sol, inst);
} /*greedy_picker*/

double grasp_picker(int i, int nearest_prob, int* sol, const instance* inst){
	if(rand() % 100 < nearest_prob)
		return minCost(i, sol, inst);
	else 
		return secMinCost(i, sol, inst);
} /*grasp_picker*/

double solve(instance* inst, int* sol, int prob, int rs, node_picker chooseNode){
	double z = 0;
	int j;
	srand(rs);
	for(j = 0; j < inst->nnodes; j++)
		sol[j] = j;
	for(j = 0; j < inst->nnodes - 2; j++)
		z += chooseNode(j, prob, sol, inst);
	z += inst->costs[sol[j] * inst->nnodes + sol[j + 1]] + inst->costs[sol[j + 1]]; /*add last two edges costs, which don't need computation*/
	updateBest(z, sol, inst);
	return z;
} /*solve*/

double testN(instance* inst, int* sol, int prob) {	
	double (*pickers[2])(int i, int nearest_prob, int* sol, const instance* inst) = 
	{
	greedy_picker,
	grasp_picker
    };
	int i;
	double z_avg = 0;
	for(i = 0; i < inst->n_sim; i++){
		z_avg += solve(inst, sol, prob, i, (node_picker)pickers[1]);
	} /*for*/
	return z_avg / inst->n_sim;
} /*testN*/

int init_stats(FILE** pstats){
	*pstats = fopen("../out/stats_GRASP.dat", "w");
	if (*pstats == NULL) {
		printf("Couldn't create the stats file!");
		return FILE_OPEN_ERR;
	} /*if*/
	fprintf(*pstats, "#PROBABILITY and AVG COST\n");
	return 0;
} /*init_stats*/

double compute_best_avg(instance* inst, int* sol, FILE* stats){
	int prob;
	double best_prob = 0;
	double z_avg_curr;
	double z_avg_best = testN(inst, sol, 0);
	fprintf(stats, "%d  %.2f\n", 0, z_avg_best);
	for(prob = 1; prob <= 100; prob++) {
		if((z_avg_curr = testN(inst, sol, prob)) < z_avg_best){
			z_avg_best = z_avg_curr;
			best_prob = prob;
		} /*if*/
		fprintf(stats, "%d  %.2f\n", prob, z_avg_curr);
	} /*for*/
	inst->best_prob = best_prob;
	return z_avg_best;
} /*compute_best_avg*/

int main(int argc, char **argv)
{
	int* sol;
	FILE* stats;
	
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	
	initInst(inst);
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
	compute_costs(inst, (cost)sq_dist);
	
	init_stats(&stats);
	printf("Average best cost %f ", compute_best_avg(inst, sol, stats));
	printf("obtained with prob. of choosing nn equals to %d\n", inst->best_prob);
	printf("Best sol.: %f\n", inst->zbest);
	system("gnuplot gnuplot_stats.p");
	
	write_out_file(inst, inst->best_sol, "h_GRASP");
    write_plotting_script("h_GRASP", inst->nnodes);
    system("gnuplot gnuplot_out.p");
		
	free(sol);
	freeInst(inst);
	return 0;
}