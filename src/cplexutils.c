#include "tsp.h"
#include "cplexutils.h"


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
/***************************************************************************************************************************/
{ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2; //can be wrong :3 
	return pos;
}

void build_model(instance *inst, CPXENVptr env, CPXLPptr lp)
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

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	/*if ( VERBOSE >= 100 )*/
	CPXwriteprob(env, lp, "model.lp", NULL);   //write the model into a file .lp 

}

void print_error(const char *err) 
{ 
	printf("\n\n ERROR: %s \n\n", err); 
	fflush(NULL); 
	exit(1); 
}     
