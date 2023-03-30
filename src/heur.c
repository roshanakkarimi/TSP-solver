#include "tsp.h"
#include "heur.h"

#define SOLV_PERC 0.2


/*greedy heuristics node pickers*/

double greedy_picker(int i, int nearest_prob, int* sol, const instance* inst){ /*just a wrapper to use node_picker standard*/
	return minCost(i, sol, inst);
} /*greedy_picker*/

double grasp_picker(int i, int nearest_prob, int* sol, const instance* inst){
	if(rand() % 100 < nearest_prob)
		return minCost(i, sol, inst);
	else 
		return secMinCost(i, sol, inst);
} /*grasp_picker*/

double (*pickers[2])(int i, int nearest_prob, int* sol, const instance* inst) = 
{
	greedy_picker,
	grasp_picker
};


/*greedy solver*/

void gr_solve(instance* inst, int* sol, int prob, node_picker chooseNode){
	double z = 0;
	int j, n = inst->nnodes;
	while(time(NULL) - inst->t_start < inst->timelimit * SOLV_PERC){
		z = 0;
		for(j = 0; j < n; j++)
			sol[j] = j;
		swapInt(sol, sol + rand() % n);
		for(j = 0; j < n - 2; j++)
			z += chooseNode(j, prob, sol, inst);
		z += inst->costs[sol[j] * n + sol[j + 1]] + inst->costs[sol[j + 1] * n + sol[0]]; /*add last two edges costs, which don't need computation*/
		updateBest(z, sol, inst);
	}
} /*gr_solve*/


/*refinement utilities*/

void reverse_sub(int* sol, int i1, int i2, int n){
	while(i1 != i2 && i1 != (i2 + 1) % n){
		swapInt(sol + i1, sol + i2);
		i1 = (i1 + 1) % n;
		i2 = (i2 + n - 1) % n;
    } /*while*/
} /*reverse_sub*/

double impr(double* costs, int* sol, int i, int j, int n){
	double old_costs = costs[sol[(i + n - 1) % n] * n + sol[i]] + costs[sol[(j + n - 1) % n] * n + sol[j]];
	double new_costs = costs[sol[(i + n - 1) % n] * n + sol[(j + n - 1) % n]] + costs[sol[i] * n + sol[j]];
	return old_costs - new_costs;
} /*impr*/


/*metaheuristics*/

double (*meta_moves[1])(const instance*, int*, int*, int) = 
{
	tabu_move
};

double tabu_move(const instance* inst, int* sol, int* tabu_list, int iter){
	return two_opt_move(inst, sol, tabu_list, iter, true);
} /*tabu_move*/


/*tabu utilities*/

int choose_tabu(int* sol, int e1, int e2, int n){
	int r = rand() % 4;
	if(r == 0) return sol[(e1 + n - 1) % n];
	else if(r == 1) return sol[e1];
	else if(r == 2) return sol[(e2 + n - 1) % n];
	else if(r == 3) return sol[e2];
	else return -1; /*just for the compiler*/
} /*choose_tabu*/


/*refinement*/

double two_opt_move(const instance* inst, int* sol, int* tabu_list, int iter, bool second_best){
	int i, j, e1, e2;
	int n = inst->nnodes;
	double curr_impr, max_impr = -1 * INFINITE_DBL;
	for(i = 0; i < n - 2; i++)
		if(iter > tabu_list[sol[i]]) /*if i is not tabu*/
			for(j = i + 2; j < n - 1 * (i == 0); j++)
				if(iter > tabu_list[sol[j]] && 
			       (curr_impr = impr(inst->costs, sol, i, j, inst->nnodes)) > max_impr) { /*if j is not tabu and better improvement*/
					max_impr = curr_impr;
					e1 = i; 
					e2 = j;
				} /*if*/
	if(max_impr > 0 || second_best) /*only after real improvements, or tabu moves*/
		reverse_sub(sol, e1, (e2 + n - 1) % n, n);
	if(inst->ref_mode == TWO_TABU && second_best){
		while(tabu_list[i = choose_tabu(sol, e1, e2, n)] > iter){} /*set tabu only if not already tabu*/
		tabu_list[i] = iter + inst->tabu_tenure;
	} /*if*/
	if(second_best)
		return max_impr; /*could be negative*/
	else if(max_impr > 0)
		return max_impr;
	return 0.0;
} /*two_opt_move*/ 

void refine(instance* inst, int* sol){ 
	int iter = 1; 
	int n = inst->nnodes;
	int* tabu_list;
	bool stop;
	double impr_curr, z_curr, end_time;
	FILE* stats = NULL;
    z_curr = getSolCost(sol, inst);
	end_time = inst->t_start + inst->timelimit * (1.0 - SOLV_PERC);
	
	tabu_list = calloc(n, sizeof(int));
	assert(tabu_list != NULL);
	
	open_init_stats(&stats, inst);
	fprintf(stats, "# ITER. and CORRESP. COST\n");
	fprintf(stats, "0 %f\n", z_curr);
	
	do{
		impr_curr = two_opt_move(inst, sol, tabu_list, iter, false);
		if(impr_curr > 0){
			z_curr -= impr_curr;
			updateBest(z_curr, sol, inst);
		} /*if*/
		else{
			if(inst->ref_mode == TWO) stop = true;
			else z_curr -= meta_moves[inst->ref_mode - 1](inst, sol, tabu_list, iter); /*ref. with metaheur. correspond to constants >= 1*/
		} /*else*/
		/*printf("Iter. %d ended with z_curr = %f\n", iter, z_curr);*/
		if(iter < 5000) fprintf(stats, "%d %f\n", iter, z_curr);
		iter++;
	}while(!stop && (time(NULL) < end_time));
	
	fclose(stats);
	free(tabu_list);
} /*refine*/


/*output files and scripts management*/

int open_init_stats(FILE** stats, const instance* inst){
	FILE* stats_script;
	char out_name[50];
	
	sprintf(out_name, "../out/%s_stats.dat", inst->fileIn);
	*stats = fopen(out_name, "w");
	if(*stats == NULL) return myError("Couldn't create the stats file!", FILE_OPEN_ERR);
	
	stats_script = fopen("stats_script.p", "w");
	if(stats_script == NULL) return myError("Couldn't create the stats script!", FILE_OPEN_ERR);
	fprintf(stats_script, "set terminal qt persist size 1000,500\n");
	fprintf(stats_script, "set xlabel \"Iteration of refinement phase\"\nshow xlabel\n");
	fprintf(stats_script, "set ylabel \"Cost of solution\"\nshow ylabel\n");
	fprintf(stats_script, "plot \"%s\" with linespoints pointtype 0\n", out_name);
	/*fprintf(stats_script, "plot \"%s\" with points pointtype 0\n", out_name);*/
	
	fclose(stats_script);
	return 0;
} /*open_init_stats*/

int open_out_files(instance* inst, FILE** out, FILE** script){
	int i;
	char out_name[50];
	sprintf(out_name, "../out/%s.dat", inst->fileIn);
	
	*out = fopen(out_name, "w");
	if(*out == NULL) return myError("Couldn't create the output file!", FILE_OPEN_ERR);
	fprintf(*out, "#MODE: %s", heur_modes[inst->heur_mode]);
	if(inst->ref_mode != NOTHING)
		fprintf(*out, " with %s", ref_modes[inst->ref_mode]);
	fprintf(*out, "\n#NODES\n");
	for(i = 0; i < inst->nnodes; i++)
		fprintf(*out, "%d %f %f\n", i, inst->pts[i].x, inst->pts[i].y);
	
	*script = fopen("gnuplot_out.p", "w");
	if (*script == NULL) return myError("Couldn't create the script!", FILE_OPEN_ERR);
	fprintf(*script, "set terminal qt persist size 500,500\n");
	fprintf(*script, "set key off\n");
	fprintf(*script, "plot \"%s\" index 0 using 2:3:1 with labels\n", out_name);
	
	return 0;
} /*open_out_files*/

int main(int argc, char **argv)
{
	int *sol;
	int i;
	FILE *out, *script;

	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	
	initInst(inst);
    parse_cmd(argc, argv, inst);
	srand(inst->randseed);
	if(strcmp(inst->fileIn, "\0") == 0) {
		rand_points(inst);
		strcpy(inst->fileIn, "random");
	} /*if*/
    else read_fileIn(inst);
	compute_costs(inst, (cost)sq_dist);
	open_out_files(inst, &out, &script);
	
	inst->t_start = time(NULL);
	gr_solve(inst, sol, inst->prob, (node_picker)pickers[inst->heur_mode]);
	update_plotting_script(script, inst->fileIn, 0);
	update_out_file(inst, inst->best_sol, out);
	
	printf("\nSol. without refinement: %f\n", inst->zbest);
	
	if(inst->ref_mode != NOTHING){
		for(i = 0; i < inst->nnodes; i++) sol[i] = inst->best_sol[i];
		inst->t_start = time(NULL);
		printf("Wait for refinement...\n");
		refine(inst, sol);
	} /*if*/
		
    update_out_file(inst, inst->best_sol, out);
	update_plotting_script(script, inst->fileIn, 1);
	
	printf("Best sol. found: %f\n", inst->zbest);
	
	fclose(out);
	fclose(script);
	
	system("gnuplot gnuplot_out.p");
	if(inst->ref_mode != NOTHING) system("gnuplot stats_script.p");
	
	free(sol);
	freeInst(inst);
	return 0;
} /*main*/