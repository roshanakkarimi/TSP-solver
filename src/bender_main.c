#include "tsp_cplex.h"
#define DEBUG if(inst->verbosity >=1)
#define DEBUGPLUS if(inst->verbosity >=10)

void plot_sol(const instance* inst, const int* arr, int size, bool lines) {
	FILE* gnuplotPipe = _popen("gnuplot -persist", "w");
	fprintf(gnuplotPipe, "set xlabel 'X'\n");
	fprintf(gnuplotPipe, "set ylabel 'Y'\n");
	
	if (lines) {
		fprintf(gnuplotPipe, "plot '-' with lines\n");
		for (int i = 0; i < size; i++) {
			fprintf(gnuplotPipe, "%lf %lf\n", inst->pts[i].x, inst->pts[i].y);
			fprintf(gnuplotPipe, "%lf %lf\n\n", inst->pts[arr[i]].x, inst->pts[arr[i]].y);
		}
	}
	else {
		fprintf(gnuplotPipe, "plot '-' using 2:3:1 with labels\n");
		for (int i = 0; i < size; i++) {
			fprintf(gnuplotPipe, "%d %lf %lf\n", arr[i], inst->pts[arr[i]].x, inst->pts[arr[i]].y);
		}
	}
	fprintf(gnuplotPipe, "e\n");

	// Close the Gnuplot pipe
	_pclose(gnuplotPipe);
}

void addSECs(const instance* inst, CPXENVptr env, CPXLPptr lp, const int* comp, int ncomp, int iter) {
	int** comp_nodes;
	int *izeros, *indexes, *nncomp;
	double *values, *rhs;
	int i, j, h, nnz = 0;
	char* senses;
	char** cnames;
	
	comp_nodes = malloc(ncomp * sizeof(int*));
	assert(comp_nodes != NULL);
	rhs = malloc(ncomp * sizeof(double));
	assert(rhs != NULL);
	izeros = calloc(ncomp, sizeof(int));
	assert(izeros != NULL);
	nncomp = calloc(ncomp, sizeof(int));
	assert(nncomp != NULL);
	senses = malloc(ncomp * sizeof(char));
	assert(senses != NULL);
	cnames = malloc(ncomp * sizeof(char*));
	assert(cnames != NULL);

	for (i = 0; i < ncomp; i++) {
		comp_nodes[i] = malloc(inst->nnodes * sizeof(int));
		assert(comp_nodes[i] != NULL);
		senses[i] = 'L';
		cnames[i] = malloc(50 * sizeof(char));
		assert(cnames[i] != NULL);
		sprintf(cnames[i], "SEC_comp_%d_iter_%d", i, iter);
	} /*for*/

	DEBUG printf("-----\n@addSECs: allocation successful!\n");

	for (i = 0; i < inst->nnodes; i++) {
		comp_nodes[comp[i] - 1][nncomp[comp[i] - 1]] = i; /*store node index in the row corresp. to its comp., in the first available pos.*/
		nncomp[comp[i] - 1]++;
	}


	DEBUG printf("@addSECs: components listed!\n");

	for (i = 0; i < ncomp; i++) {
		nnz += nncomp[i] * (nncomp[i] - 1) / 2; /*possible edges connecting nodes in comp. i*/
		rhs[i] = nncomp[i] - 1;
		if (i > 0) izeros[i] = izeros[i - 1] + nncomp[i - 1] * (nncomp[i - 1] - 1) / 2;
	} /*for*/

	DEBUG printf("@addSECs: nnz, izeros, rhs computed! nnz = %d\n", nnz);

	values = malloc(nnz * sizeof(double));
	assert(values != NULL);
	indexes = malloc(nnz * sizeof(int));
	assert(indexes != NULL);

	nnz = 0;
	for (i = 0; i < ncomp; i++)
		for (j = 0; j < nncomp[i]; j++)
			for (h = j + 1; h < nncomp[i]; h++) {
				indexes[nnz] = xpos(comp_nodes[i][j], comp_nodes[i][h], inst);
				values[nnz++] = 1.0;
			} /*for*/

	DEBUG printf("@addSECs: indexes written!\n");

	if(CPXaddrows(env, lp, 0, ncomp, nnz, rhs, senses, izeros, indexes, values, NULL, cnames)) print_error("CPXaddrows(): error 1");
	DEBUG CPXwriteprob(env, lp, "model.lp", NULL); /*remove from here!*/

	for (i = 0; i < ncomp; i++) {
		free(comp_nodes[i]);
		free(cnames[i]);
	} /*for*/
	free(comp_nodes);
	free(indexes);
	free(cnames);
	free(senses);
	free(values);
	free(rhs);
	free(izeros);
	free(nncomp);
} /*addSECs*/

int main(int argc, char** argv) {
	int *succ, *comp;
	double *xstar;
	int ncomp, error, i = 0;
	CPXENVptr env;
	CPXLPptr lp;
	instance* inst;
	int ncols;
	
	inst = malloc(sizeof(instance));
	assert(inst != NULL);
	
	env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if (error) print_error("CPXcreateprob() error");
			
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

	build_model(inst, env, lp);
	DEBUG printf("@main: CPLEX model built successfully!\n");
	/* Cplex's parameter setting */
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	DEBUGPLUS CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); /* Cplex output on screen */
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->randseed);
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);

	succ = (int*)calloc(inst->nnodes, sizeof(int));
	assert(succ != NULL);
	comp = (int*)calloc(inst->nnodes, sizeof(int));
	assert(comp != NULL);
	ncols = CPXgetnumcols(env, lp);
	xstar = (double*)calloc(ncols, sizeof(double));
	assert(xstar != NULL);
	
	inst->t_start = time(NULL);
	while (true) {
		if(time(NULL) - inst->t_start >= inst->timelimit) break;
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit - (time(NULL) - inst->t_start));
		error = CPXmipopt(env,lp);
		DEBUG printf("-----\n@main @bender_loop: EXECUTED CPXmipopt!\n");
		if (error) 
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}
		if (CPXgetx(env, lp, xstar, 0, ncols-1)) print_error("CPXgetx() error");	
		build_sol(xstar, inst, succ, comp, &ncomp);
		DEBUG printf("@main @bender_loop: SOLUTION BUILT SUCCESSFULLY! %d COMPS\n", ncomp);
		if(ncomp == 1) break;
		addSECs(inst, env, lp, comp, ncomp, i); /*add SECs for each comp.*/
		DEBUG printf("@main @bender_loop: ADDED SECs FOR %d COMPS AT ITER %d!\n", ncomp, i);
		i++;
	} /*while*/

	CPXwriteprob(env, lp, "model.lp", NULL);
	DEBUGPLUS plot_sol(inst, succ, inst->nnodes, false);
	plot_sol(inst, succ, inst->nnodes, true);

	
	/* free and close cplex model  */ 
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 
	
	free(xstar);
	freeInst(inst);
	return 0;
} /*main*/
