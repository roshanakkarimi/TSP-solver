#include "heur.h"

#define N_COMP 4


int main(int argc, char **argv){
	int i, j;
	int* sol;
	FILE* stats;
	instance* inst;
	double zbests[N_COMP], t;
	char comp_modes[N_COMP][100] =
	{
		"determ. greedy",
		"GRASP greedy",
		"GRASP greedy with 2opt ref.",
		"GRASP greedy with 2opt ref. and tabu"
		
	};
	
	if(argc < 2) 
		return myError("Usage: performance {number of random instances to be gen.}", COMM_ERR);
	
	inst = malloc(sizeof(instance));
	assert(inst != NULL);
	sol = malloc(N_DEF_NODES * sizeof(int));
	assert(sol != NULL);
	
	/*stats file initialization*/
	stats = fopen("../pp/heur_perf.csv", "w");
	if(stats == NULL) 
		return myError("Couldn't create the stats file!", FILE_OPEN_ERR);
	fprintf(stats, "%d", N_COMP);
	for(i = 0; i < N_COMP; i++)
		fprintf(stats, ", %s", comp_modes[i]);
	
	initInst(inst);
	alloc_inst(inst);
	
	inst->timelimit = 5; /*input?*/
	
	srand(inst->randseed);
	
	for(i = 0; i < atoi(argv[1]); i++){		
		rand_points(inst);
		compute_costs(inst, (cost)sq_dist);
		printf("---------------------------\nStarting elaboration on instance %d out of %d\n", i+1, atoi(argv[1]));
		t = time(NULL);
		
		inst->zbest = INFINITE_DBL; /*greedy initializ.*/
		inst->t_start = time(NULL);
		gr_solve(inst, sol, (node_picker)greedy_picker);
		zbests[0] = inst->zbest;
		
		inst->zbest = INFINITE_DBL; /*grasp greedy initializ.*/
		inst->t_start = time(NULL);
		gr_solve(inst, sol, (node_picker)grasp_picker);
		zbests[1] = inst->zbest;
		
		inst->ref_mode = TWO;
		refine(inst, sol, false);
		zbests[2] = inst->zbest;
		
		inst->t_start = time(NULL);
		inst->ref_mode = TWO_TABU;
		refine(inst, sol, false);
		zbests[3] = inst->zbest;
		
		printf("Elaboration of instance %d out of %d ended in %f secs\n", i+1, atoi(argv[1]), time(NULL)-t);
		fprintf(stats, "\nrand_%d", i);
		for(j = 0; j < N_COMP; j++){
			fprintf(stats, ", %f", zbests[j]);
			printf("%s: %.2f\n", comp_modes[j], zbests[j]);
		} /*for*/
	}/*for*/
	
	fclose(stats);
	free(sol);
	freeInst(inst);
	return 0;
} /*main*/