#include "tsp.h"
#include "cplexutils.h"
#include <stdbool.h>


void initCplex(instance *inst, CPXENVptr env)
{
	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); //CPX_OFF => nothing to show of Cplex outputs
	if ( inst->verbosity >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->randseed);	//int
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit); //double
	// ...
}     

int xpos(int i, int j, instance *inst)      // to be verified                                           
{ 
	if ( i == j ) printf(" %d == %d in xpos", i, j);
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2; 
	return pos;
}

void build_model(instance *inst, CPXENVptr env, CPXLPptr lp, bool loop, int *succ, int *comp, int *ncomp, const char *mlname)
{    
	//double zero = 0.0;  
	char binary = 'B';

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

// add binary var.s x(i,j) for i < j  

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);  		// ... x(1,2), x(1,3) ....
			double obj = dist(i,j,inst->pts); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error(" wrong CPXnewcols on x var.s"); //&binary shows if the variable is binary or not we use & becuase cplex need a pointer so we use this trick that is also why cname is an array of characters
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) ) print_error(" wrong position for x var.s"); //automatic check for xpos trick
		
		}
	} 

// add the degree constraints 

	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));

	for ( int h = 0; h < inst->nnodes; h++ )  		// add the degree constraint on node h
	{
		double rhs = 2.0;
		char sense = 'E';                            // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h+1);   
		int nnz = 0;
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i,h, inst);
			value[nnz] = 1.0;
			nnz++;
		}
		int izero = 0;
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1"); // 0 means not having only 1 row at every iteration and no extra variables , nnz = non zeros, 
	} 

// add the subtour elimination constraints
	if (*ncomp>1){
		int *index2 = (int *) calloc(inst->nnodes, sizeof(int));
		double *value2 = (double *) calloc(inst->nnodes, sizeof(double));
		int *count_comp = (int *) calloc(*ncomp, sizeof(int));
		for (int i = 0; i < *ncomp ; i++)
		{
			for (int h = 0; h < inst->nnodes; h++)
			{
				if(comp[h] == i)
				{
					count_comp[i] += 1;
				}
			}
			char sense = 'L';
			double rhs = count_comp[i];
			int nnz = 0;
			for (int j = 0; j < inst->nnodes; j++)
			{
				if(comp[j] == i)
				{
					index2[nnz] = xpos(j, succ[j], inst);
					value2[nnz] = 1;
					nnz++;
					int izero = 0;
					if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index2, value2, NULL, &cname[0]) ) print_error("CPXaddrows(): error 1");
				}
			}
		}
		free(count_comp);
		free(value2);
		free(index2);
	}
	
	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	/*if ( VERBOSE >= 100 )*/
	CPXwriteprob(env, lp, mlname, NULL);   //write the model into a file .lp 

}

void print_error(const char *err) 
{ 
	printf("\n\n ERROR: %s \n\n", err); 
	fflush(NULL); 
	exit(1); 
}     

void CPXsolver(instance *inst, CPXENVptr env, CPXLPptr lp){
	int error;
	error = CPXmipopt(env,lp); //optimize the mixed integer program env, lp
	if ( error ) 
	{
		printf("CPX error code %d\n", error);
		print_error("CPXmipopt() error"); 
	} // infeasibility is not an error
	
}

