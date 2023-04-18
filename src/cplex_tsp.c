#include "tsp.h"
#include "cplexutils.h"
#include <cplex.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#define FILE_OPEN_ERR 1
#define FILE_STRUCT_ERR 2
#define FORMAT_ERR 3
#define INDEX_ERR 4
#define DEBUG    // comment out to avoid debugging 
#define EPS 1e-5

void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
{   
	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;
		int done = 0;
		while ( !done )  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if ( i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;  // last arc to close the cycle
		
		// go to the next component...
	}
}

int main(int argc, char **argv)
{
	// initializing the tsp instance
	instance* inst = malloc(sizeof(instance));
	int* succ = malloc(sizeof(int)*inst->nnodes);
	int* comp = malloc(sizeof(int)*inst->nnodes);
	int ncomp;

//first try
	initInst(inst);
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
  
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error); // creating environment of problem (pointer), error = 0 if no problems (enough memory and ,...)
	if ( error ) print_error("CPXopenCPLEX() error"); //
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); //empty linear program model (pointer)
	if ( error ) print_error("CPXcreateprob() error"); 

	build_model(inst, env, lp, false, succ, comp, &ncomp, "model1.lp");

	// setting up the cplex parameters
	initCplex(inst,env);
	
	CPXsolver(inst, env, lp);
	
	// use the optimal solution found by CPLEX
	int ncols = CPXgetnumcols(env, lp);
	double *xstar = (double *) calloc(ncols, sizeof(double));
	if ( CPXgetx(env, lp, xstar, 0, ncols-1) /*from 0 to ncols-1 of xstar will be filled */ ) print_error("CPXgetx() error");	
	
	build_sol(xstar, inst, succ, comp, &ncomp);

	// free and close cplex model   
	CPXfreeprob(env, &lp); //based on the sequence (reversed)
	CPXcloseCPLEX(&env); 

//after component try
	if ( ncomp > 1){
	// open CPLEX model
		int error;
		CPXENVptr env = CPXopenCPLEX(&error); // creating environment of problem (pointer), error = 0 if no problems (enough memory and ,...)
		if ( error ) print_error("CPXopenCPLEX() error"); //
		CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); //empty linear program model (pointer)
		if ( error ) print_error("CPXcreateprob() error"); 

		build_model(inst, env, lp, false, succ, comp, &ncomp, "model2.lp");

	// setting up the cplex parameters
		initCplex(inst,env);
		CPXsolver(inst, env, lp);
	
	// use the optimal solution found by CPLEX
		int ncols = CPXgetnumcols(env, lp);
		double *xstar = (double *) calloc(ncols, sizeof(double));
		if ( CPXgetx(env, lp, xstar, 0, ncols-1) /*from 0 to ncols-1 of xstar will be filled */ ) print_error("CPXgetx() error");	
	

	// free and close cplex model   
		CPXfreeprob(env, &lp); //based on the sequence (reversed)
		CPXcloseCPLEX(&env); 
	}
	
	free(succ);
	free(comp);
	free(xstar);
	freeInst(inst);
}

/*
for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			if ( xstar[xpos(i,j,inst)] > 0.5 ) printf("  ... x(%3d,%3d) = 1\n", i+1,j+1); // xpos(i,j,inst) the position in cplex of xij, x* = ?
            // > 0.5 because approximate arithmetic is not exactly 0 or 1 but in a specfific threshold > 0.5 equals to 1 and otherwise is 0
		}
	}
	*/
/*
**** LAZY CONTRAINTS IN THE INPUT MODEL ****

Ex: MZT formulation with directed-arc variables x_ij and x_ji --> xpos_compact(i,j,inst)

...

	int izero = 0;
	int index[3]; 
	double value[3];

	// add lazy constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0	
	double big_M = inst->nnodes - 1.0;
	double rhs = big_M -1.0;
	char sense = 'L';
	int nnz = 3;
	for ( int i = 1; i < inst->nnodes; i++ ) // excluding node 0
	{
		for ( int j = 1; j < inst->nnodes; j++ ) // excluding node 0 
		{
			if ( i == j ) continue;
			sprintf(cname[0], "u-consistency for arc (%d,%d)", i+1, j+1);
			index[0] = upos(i,inst);	
			value[0] = 1.0;	
			index[1] = upos(j,inst);
			value[1] = -1.0;
			index[2] = xpos_compact(i,j,inst);
			value[2] = big_M;
			if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints() for u-consistency");
		}
	}
	
	// add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1, for each arc (i,j) with i < j
	rhs = 1.0; 
	char sense = 'L';
	nnz = 2;
	for ( int i = 0; i < inst->nnodes; i++ ) 
	{
		for ( int j = i+1; j < inst->nnodes; j++ ) 
		{
			sprintf(cname[0], "SEC on node pair (%d,%d)", i+1, j+1);
			index[0] = xpos_compact(i,j,inst);
			value[0] = 1.0;
			index[1] = xpos_compact(j,i,inst);
			value[1] = 1.0;
			if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) print_error("wrong CPXlazyconstraints on 2-node SECs");
		}
	}

...

*** SOME MAIN CPLEX'S PARAMETERS ***


	// increased precision for big-M models
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);		// very important if big-M is present
	CPXsetdblparam(env, CPX_PARAM_EPRHS, 1e-9);   						

	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
	if ( VERBOSE >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
	
	CPXsetdblparam(env, CPX_PARAM_TILIM, CPX_INFBOUND+0.0); 
	
	CPXsetintparam(env, CPX_PARAM_NODELIM, 0); 		// abort Cplex after the root node
	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);	// abort Cplex after the first incumbent update

	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-4);  	// abort Cplex when gap below 0.01%	 
	

	
*** instance TESTBED for exact codes:

		all TSPlib instances with n <= 500	
*/
