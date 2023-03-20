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

void reverse_sub(int* sol, int i1, int i2, int n){
	while(i1 != i2 && i1 != i2 + 1){
		swapInt(sol + i1, sol + i2);
		i1 = (i1 + 1) % n;
		i2 = (i2 - 1) % n;
    } /*while*/
} /*reverse_sub*/

double two_opt(instance* inst, int* sol){
	int i, j, e1, e2;
	int n = inst->nnodes;
	double curr_impr, max_impr = 0;
	int first = sol[0];
	int last = sol[n - 1];
	double last_cost = inst->costs[first * n + last];
	for(i = 1; i < n - 2; i++)
		for(j = i + 2; j < n; j++)
			if((curr_impr = inst->costs[sol[i - 1] * n + sol[i]] + inst->costs[sol[j - 1] * n + sol[j]] 
			             - inst->costs[sol[i - 1] * n + sol[j - 1]] - inst->costs[sol[i] * n + sol[j]]) > max_impr) {
				max_impr = curr_impr;
				e1 = i; 
				e2 = j;
			} /*if*/
	for(i = 2; i < n - 1; i++)
		if((curr_impr = inst->costs[sol[i - 1] * n + sol[i]] + last_cost 
			            - inst->costs[sol[i - 1] * n + last] - inst->costs[sol[i] * n + first]) > max_impr) {
				max_impr = curr_impr;
				e1 = i; 
				e2 = 0;
			} /*if*/
	if(max_impr > 0) {
		printf("Switched edges ending in %d and %d\n", sol[e1], sol[e2]);
		reverse_sub(sol, e1, (e2 - 1) % n, n);
	} /*if*/
	return max_impr;
} /*two_opt*/

int main(int argc, char **argv)
{
	double (*pickers[2])(int i, int nearest_prob, int* sol, const instance* inst) = 
	{
	greedy_picker,
	grasp_picker
	
	
    };
	
	double impr;
	int* sol;
	int ord = 0;
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	
	initInst(inst);
    parse_cmd(argc, argv, inst);
	read_fileIn(inst);
	compute_costs(inst, (cost)sq_dist);
		
	solve(inst, sol, inst->prob, inst->randseed, (node_picker)pickers[inst->mode]);
	write_plotting_script("h_GRASP", inst->nnodes, ord++);
	write_out_file(inst, sol, "h_GRASP", "w");
	
	if(inst->two_opt) while((impr = two_opt(inst, sol))) {
		updateBest(inst->zbest - impr, sol, inst);
		write_out_file(inst, sol, "h_GRASP", "a");
		write_plotting_script("h_GRASP", inst->nnodes, ord++);
	} /*while*/
	
	printf("Best sol.: %f\n", inst->zbest);
	system("gnuplot gnuplot_out.p");
	
	free(sol);
	freeInst(inst);
	return 0;
} /*main*/