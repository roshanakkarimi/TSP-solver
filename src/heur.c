#include "heur.h"


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
		alloc_inst(inst);
		rand_points(inst);
		strcpy(inst->fileIn, "random");
	} /*if*/
    else read_fileIn(inst);
	compute_costs(inst, (cost)sq_dist);
	open_out_files(inst, &out, &script);
	
	inst->t_start = time(NULL);
	gr_solve(inst, sol, (node_picker)pickers[inst->heur_mode]);
	update_plotting_script(script, inst->fileIn, 0);
	update_out_file(inst, inst->best_sol, out);
	
	printf("\nSol. without refinement: %f\n", inst->zbest);
	
	if(inst->ref_mode != NOTHING){
		for(i = 0; i < inst->nnodes; i++) sol[i] = inst->best_sol[i];
		inst->t_start = time(NULL);
		printf("Wait for refinement...\n");
		refine(inst, sol, true);
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