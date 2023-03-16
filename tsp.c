#include "tsp.h"
#include <stdbool.h>
#include <stdlib.h>

/*utilities*/

double ptdist(point* pt1, point* pt2) {
	return sqrt(pow(pt1->x - pt2->x, 2) + pow(pt1->y - pt2->y, 2)); /* to be changed using C builtin ecludian function */
} /*ptdist*/

double dist(int i, int j, instance* inst) {
	return ptdist(&inst->pts[i], &inst->pts[j]);
} /*dist*/


/*managing errors*/

int myError(const char* err, int errType){
	printf("\nFATAL ERROR:\n%s", err);
	return errType;
} /*myError*/


/* input elaboration*/

void cmdHelp() {
	printf("\n--Available options:--\n -f to set input file \n -tl to set overall time limit\n -rs to set random seed\n -help to see options\n");
} /*help*/

void dispPars(const instance* inst) {
	printf("\n--Current parameters:--\n");
	printf(" time limit: %f\n", inst->timelimit);
	printf(" randseed: %d\n", inst->randseed);
	printf(" verbosity level: %d\n", inst->verbosity);
	printf(" input file: %s\n", inst->fileIn);
} /*dispPars*/

void parse_cmd(int argc, char** argv, instance *inst) {
	int i = 0;
	bool invalid_opt;
	strcpy(inst->fileIn, "NOT_SET");
	inst->verbosity = 1;
	inst->randseed = DEFAULT_RAND;
	inst->timelimit = INFINITE;
	while(i < argc - 1) {
		i++;
		invalid_opt = true;
		if(strcmp(argv[i],"-f") == 0){
			strcpy(inst->fileIn, argv[++i]);
			invalid_opt = false;
			continue;
		} /*set input file*/
		if(strcmp(argv[i],"-tl") == 0 ){ 
			inst->timelimit = atof(argv[++i]);
			invalid_opt = false;
			continue;
		} /*set time limit*/
		if(strcmp(argv[i],"-rs") == 0){
			inst->randseed = abs(atoi(argv[++i])); 
			invalid_opt = false;
			continue;
		} /*set random seed*/
		if(strcmp(argv[i],"-v") == 0){
			inst->verbosity = atoi(argv[++i]);
			invalid_opt = false;
			continue;
		} /*set verbosity param.*/
		if(strcmp(argv[i],"-disp") == 0) {
			dispPars(inst);
			invalid_opt = false;
			continue;
		} /*print the values of the parameters*/
		if(strcmp(argv[i],"-help") == 0) invalid_opt = false;
		else if(invalid_opt) printf("\nINVALID OPTION:");
		cmdHelp();
	} /*while*/
} /*parse_cmd*/

int read_fileIn(instance* inst) {
	
	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	bool node_sect = false;
	
	FILE* fileIn = fopen(inst->fileIn, "r");
	if(fileIn == NULL) return myError("Couldn't open the file!", FILE_OPEN_ERR);
	
	inst->nnodes = -1;

	while (fgets(line, sizeof(line), fileIn) != NULL) 
	{
		if (strlen(line) <= 1 ) continue; /*skip empty lines*/
	    par_name = strtok(line, " :");
		
		if (strncmp(par_name, "NAME", 4) == 0 || strncmp(par_name, "COMMENT", 7) == 0) continue;
		if (strncmp(par_name, "TYPE", 4) == 0) 
		{
			token1 = strtok(NULL, " :");  
			if(strncmp(token1, "TSP", 3) != 0) return myError("Format error: only TYPE == CSP allowed!", FORMAT_ERR); 
			continue;
		} /*if*/
		if (strncmp(par_name, "DIMENSION", 9) == 0) 
		{
			if(inst->nnodes >= 0) return myError("Repeated DIMENSION section in input file!", FILE_STRUCT_ERR);
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1); 
			inst->pts = malloc(inst->nnodes * sizeof(point));      
			continue;
		} /*if*/
		if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0) 
		{
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "ATT", 3) != 0) return myError("Format error: only EDGE_WEIGHT_TYPE == ATT allowed!", FORMAT_ERR); 
			continue;
		} /*if*/            
		if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0) 
		{
			if (inst->nnodes <= 0) return myError("DIMENSION section should appear before NODE_COORD_SECTION section!", FILE_STRUCT_ERR);
			node_sect = true;   
			continue;
		} /*if*/
		if (strncmp(par_name, "EOF", 3) == 0) 
		{
			node_sect = false;
			break;
		} /*if*/

		if (node_sect)
		{
			int i = atoi(par_name) - 1; 
			if (i < 0 || i >= inst->nnodes) return myError("Unknown node in NODE_COORD_SECTION section!", INDEX_ERR);     
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->pts[i].x = atof(token1);
			inst->pts[i].y = atof(token2); 
			continue;
		} /*if*/
		
	} /*while*/
	
	/*initializing distance matrix*/
	calc_distance_matrix(inst);

	fclose(fileIn);
	return 0;
} /*read_fileIn*/

/*output procedures*/

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

double solution_value(instance* inst, int arr[]) {
	double total_cost = 0;
	for ( int i=0; i<inst->nnodes; i++){
		total_cost += inst->distance_matrix[arr[i]][arr[i+1]];
	}
	return total_cost;
}

/*tsp huristics*/
void calc_distance_matrix(instance* inst) {

	double** distance_matrix = (double **)malloc(sizeof(double) * inst->nnodes);
	for ( int i =0; i< inst->nnodes; i++){
		distance_matrix[i] = (double *)malloc(sizeof(double) * inst->nnodes);
	}

	for(int i=0; i< inst->nnodes; i++){
		for(int j=0; j< inst->nnodes; j++){
			distance_matrix[i][j] = ptdist(&inst->pts[i], &inst->pts[j]);
		} 
	}

	inst->distance_matrix = distance_matrix;

	/*for (int i = 0; i < inst->nnodes; i++) {
		for ( int j=0; j < inst->nnodes; j++){
			free(distance_matrix[i]);
		}
    } */
}

void greedy(instance* inst) {
	
	int tmp;
	double current_dist;
	int* solution_sequence = malloc(sizeof(int) * inst->nnodes + 1); 

	for (int i=0; i < inst->nnodes; i++){
		solution_sequence[i] = i;
	} /*initialization of the solution sequence*/
		
	for (int i=0; i < inst->nnodes -1; i++){

		int best_index;
		double min_dist = INFINITE;

		for (int j=i+1; j < inst->nnodes; j++){

			current_dist = inst->distance_matrix[solution_sequence[i]][solution_sequence[j]];
			if (current_dist < min_dist ){
				best_index = j;
				min_dist = current_dist;
			}

		}
		tmp = solution_sequence[i+1];
		solution_sequence[i+1] = solution_sequence[best_index];
		solution_sequence[best_index] = tmp;
	}

	solution_sequence[inst->nnodes] = solution_sequence[0];
	plot_sol(inst, solution_sequence, inst->nnodes);
	solution_value(inst, solution_sequence);

	/*print outputs*/
	printf("%f,", solution_value(inst, solution_sequence));
	for (int i=0; i<= inst->nnodes; i++){
		printf("%d,", solution_sequence[i]);
	}
	/*end*/
	
	free(solution_sequence);
}

/*freeing the memory*/

void freeInst(instance* inst) {
	free(inst->pts);
	free(inst);
} /*freeInst*/
