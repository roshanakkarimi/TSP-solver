#include "tsp.h"

#define SOLV_PERC 0.2

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
} /*solve*/

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

int choose_tabu(int* sol, int e1, int e2, int n){
	int r = rand() % 4;
	if(r == 0) return sol[(e1 + n - 1) % n];
	else if(r == 1) return sol[e1];
	else if(r == 2) return sol[(e2 + n - 1) % n];
	else if(r == 3) return sol[e2];
	else return -1; /*just for the compiler*/
} /*choose_tabu*/

double two_opt_move(instance* inst, int* sol, int* tabu_list, int iter, bool not_impr){
	int i, j, e1, e2;
	int n = inst->nnodes;
	double curr_impr, max_impr = -1 * INFINITE_DBL;
	for(i = 0; i < n - 2; i++)
		if(iter > tabu_list[sol[i]])
			for(j = i + 2; j < n - 1 * (i == 0); j++)
				if(iter > tabu_list[sol[j]] &&
			       (curr_impr = impr(inst->costs, sol, i, j, inst->nnodes)) > max_impr) {
					max_impr = curr_impr;
					e1 = i; 
					e2 = j;
				} /*if*/
	if(max_impr > 0 || not_impr)
		reverse_sub(sol, e1, (e2 + n - 1) % n, n);
	if(inst->tabu && not_impr)
		tabu_list[(i = choose_tabu(sol, e1, e2, n))] = iter + inst->tabu_tenure;
	return max_impr * (max_impr > 0 || not_impr);
} /*two_opt*/

void refine(instance* inst, int* sol){ /*remove tabu*/
	int iter = 1; 
	int n = inst->nnodes;
	int* tabu_list;
	bool t = inst->tabu;
	bool not_impr = false;
	double impr_curr, z_curr;
    z_curr = getSolCost(sol, inst);
	
	tabu_list = calloc(n, sizeof(int));
	assert(tabu_list != NULL);
	
	do{
		impr_curr = two_opt_move(inst, sol, tabu_list, iter, not_impr);
		z_curr -= impr_curr;
		if(!not_impr) updateBest(z_curr, sol, inst);
		not_impr = (impr_curr <= 0);
		iter++;
	}while((t && (time(NULL) - inst->t_start < inst->timelimit * (1.0 - SOLV_PERC))) || (!t && !not_impr));
	free(tabu_list);
} /*refine*/

/*int open_out_files(instance* inst, FILE* out, FILE* script){
	char out_name[50];
	sprintf(out_name, "%s_%s%s", inst->fileIn, mods[inst->mode], (inst->two_opt ? "_2opt" : ""));
	script = fopen("gnuplot_out.p", "w");
	if (script == NULL) return myError("Couldn't create the script!", FILE_OPEN_ERR);
	out = fopen(fname, mode);
	if(out == NULL) return myError("Couldn't create the output file!", FILE_OPEN_ERR);
	return 0;
} *//*open_out_files*/

int main(int argc, char **argv)
{
	int *sol;
	int i;
	char out_name[50];
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	sol = malloc(inst->nnodes * sizeof(int));
	assert(sol != NULL);
	
	initInst(inst);
    parse_cmd(argc, argv, inst);
	srand(inst->randseed);
	if(strcmp(inst->fileIn, "rand") == 0) rand_points(inst); /*set as input*/
    else read_fileIn(inst);
	compute_costs(inst, (cost)sq_dist);
	
	sprintf(out_name, "%s_%s%s", inst->fileIn, mods[inst->mode], (inst->two_opt ? "_2opt" : ""));
	
	inst->t_start = time(NULL);
	gr_solve(inst, sol, inst->prob, (node_picker)pickers[inst->mode]);
	write_plotting_script(out_name, inst->nnodes, 0);
	write_out_file(inst, sol, out_name, "w");
	
	printf("\nSol. without 2opt: %f\n", inst->zbest);
	
	for(i = 0; i < inst->nnodes; i++) sol[i] = inst->best_sol[i];
	inst->t_start = time(NULL);
	if(inst->two_opt)
		refine(inst, sol);
		
    write_out_file(inst, sol, out_name, "a");
	write_plotting_script(out_name, inst->nnodes, 1);
	
	printf("Best sol. found: %f\n", inst->zbest);
	system("gnuplot gnuplot_out.p");
	
	free(sol);
	freeInst(inst);
	return 0;
} /*main*/