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

/*ROSHAA*/

double min(double* list, int size){
	double minValue=INFINITE;
	for (int i=0; i < size; i++){
		if(list[i] < minValue){
			minValue = list[i];
		}
	}

	return minValue;
}

double max(double* list, int size){
	double maxValue=0;
	for (int i=0; i < size; i++){
		if(list[i] > maxValue){
			maxValue = list[i];
		}
	}

	return maxValue;
}
void plot_sol(instance* inst, int arr[], int size) {

	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	fprintf(gnuplotPipe, "set xlabel 'X'\n");
	fprintf(gnuplotPipe, "set ylabel 'Y'\n");

	fprintf(gnuplotPipe, "plot '-' using 2:3:1 with linespoints\n");
	for (int i = 0; i <= size; i++) {
		fprintf(gnuplotPipe, "%d %lf %lf\n", i, inst->pts[arr[i]].x, inst->pts[arr[i]].y);
	}
	fprintf(gnuplotPipe, "e\n");

	// Close the Gnuplot pipe
	pclose(gnuplotPipe);
}

void furthest_initialization(instance* inst, int* arr){ /*make it more functional*/
	double maxDistance=0;
	int* extrem = malloc(sizeof(int)*inst->nnodes);
	double* xArr = malloc(sizeof(int)*inst->nnodes);
	double* yArr = malloc(sizeof(int)*inst->nnodes);

	for (int i=0; i < inst->nnodes; i++){
		xArr[i] = inst->pts[i].x;
	} 

	for (int i=0; i < inst->nnodes; i++){
		yArr[i] = inst->pts[i].y;
	}
	

	double minX = min(xArr, inst->nnodes);
	double maxX = max(xArr, inst->nnodes);
	double minY = min(yArr, inst->nnodes);
	double maxY = max(yArr, inst->nnodes); 

	for (int i = 0; i < inst->nnodes; i++)
	{	
		if (inst->pts[i].x == minX || 
			inst->pts[i].x == maxX ||
			inst->pts[i].y == minY ||
			inst->pts[i].y == maxY) extrem[i] = 1;
	} /* finding extremum points */

	for (int i = 0; i <inst->nnodes; i++){
		for (int j = 0; j <inst->nnodes; j++){
			if (extrem[i] * extrem[j] == 1){ /* both be extremum */
				double currentDistance = dist(i, j, inst->pts); /* extrem = 0 or 1*/
				if (currentDistance > maxDistance){
					arr[0] = i;
					arr[1] = j;
					maxDistance = currentDistance;
				}
			}
		}
	} /*computing maximum distance*/

	free(xArr);
	free(yArr);
	free(extrem);
}

bool unvisited_nodes(int* prev, int size){
	for (int i=0; i<size; i++){
		if (prev[i] == -1){
			return true;
		}
	}
	return false;
}

double extra_mileage(instance* inst){

	int firstNodeIndex =0, secondNodeIndex=0, visitedNode, newNode;
	//int two_first[2];

	int* prev = malloc(sizeof(int)*inst->nnodes);

	double minCost, currentCost;

	/*nodes structure*/
	for (int i=0; i<inst->nnodes; i++){
		prev[i] = -1;
	}

	/*finding furthest points*/
	//furthest_initialization(inst, two_first);
	firstNodeIndex = 6;
	secondNodeIndex = 16; 

	prev[secondNodeIndex] = firstNodeIndex;
	prev[firstNodeIndex] = secondNodeIndex;

	/*finding closest next option*/
	while(unvisited_nodes(prev, inst->nnodes)){
		minCost = INFINITE;
		for (int i=0; i<inst->nnodes; i++){
			if (prev[i] == -1){ /*unvisited*/
				for (int j=0; j<inst->nnodes; j++){
					if (prev[j] != -1){ /*visited*/
						currentCost = inst->costs[prev[j] * inst->nnodes + i] + inst->costs[i * inst->nnodes + j];
						if(currentCost < minCost){
						minCost = currentCost;
						newNode = i;
						visitedNode = j;
						}
					}
				}
			}
		}
		prev[newNode] = prev[visitedNode];
		prev[visitedNode] = newNode;
	}

	int* solution_sequence = malloc(sizeof(int) * inst->nnodes);
	int curr_index = 6;
	for (int i=0; i < inst->nnodes; i++){
		solution_sequence[i] = prev[curr_index];
		curr_index = prev[curr_index];
	} 	/*creating solution sequence*/
	plot_sol(inst, solution_sequence, inst->nnodes);

	free(prev);
	free(solution_sequence);
	return 0.0;

}


double testN(instance* inst, int* sol, int prob) {
	int i;
	double z_avg = 0;
	for(i = 0; i < inst->n_sim; i++){
		z_avg += solve(inst, sol, prob, i, grasp_picker);
	} /*for*/
	return z_avg / inst->n_sim;
} /*testN*/

int main(int argc, char **argv)
{
	int* sol;
	bool test;
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	
	initInst(inst);
	test = parse_cmd(argc, argv, inst);
	read_fileIn(inst);
	compute_costs(inst, (cost)sq_dist);
	
	if(test){
		int prob;
		double z_avg_best, z_avg_curr;
		FILE* stats = fopen("../out/stats_GRASP.dat", "w");
		if (stats == NULL) {
			printf("Couldn't create the stats file!");
			return FILE_OPEN_ERR;
		} /*if*/
		fprintf(stats, "#PROBABILITY and AVG COST\n");
	
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
		system("gnuplot gnuplot_stats.p");
	} /*if*/
	
	else
		//solve(inst, sol, inst->prob, inst->randseed, (node_picker)((inst->mode == GREEDY) ? greedy_picker : grasp_picker));
		extra_mileage(inst);
	printf("Best sol.: %f\n", inst->zbest);
	write_out_file(inst, inst->best_sol, "h_GRASP");
    write_plotting_script("h_GRASP", inst->nnodes);
    system("gnuplot gnuplot_out.p");
	
	free(sol);
	freeInst(inst);
	return 0;
}
