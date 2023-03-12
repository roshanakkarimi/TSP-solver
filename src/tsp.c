#include "tsp.h"
#include <stdbool.h>

/*utilities*/

double dist(int i, int j, const point* pts) {
	return sqrt(pow(pts[i].x - pts[j].x, 2) + pow(pts[i].y - pts[j].y, 2));
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

void initInst(instance *inst) {
	strcpy(inst->fileIn, "../data/");
	inst->verbosity = 1;
	inst->randseed = DEFAULT_RAND;
	inst->timelimit = INFINITE;
	inst->zbest = -1;
} /*initInst*/

void parse_cmd(int argc, char** argv, instance *inst) {
	int i = 0;
	bool invalid_opt;
	while(i < argc - 1) {
		i++;
		invalid_opt = true;
		if(strcmp(argv[i],"-f") == 0){
			strcat(inst->fileIn, argv[++i]);
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
	if(fileIn == NULL) return myError("Couldn't open the input file!", FILE_OPEN_ERR);
	
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
			assert(inst->pts != NULL);
			inst->costs = malloc(inst->nnodes * inst->nnodes * sizeof(double));
			assert(inst->costs != NULL);
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
	
	fclose(fileIn);
	return 0;
} /*read_fileIn*/

void compute_costs(instance* inst){
	int i, j;
	for(i = 0; i < inst->nnodes; i++)
		for(j = 0; j < inst->nnodes; j++)
			inst->costs[i * inst->nnodes + j] = (i == j) * INFINITE + dist(i, j, inst->pts);
} /*compute_costs*/


/*output elaboration*/

int write_plotting_script(const instance* inst){
	FILE *script = fopen("gnuplot_out.p", "w");
	if (script == NULL) return myError("Couldn't open the script!", FILE_OPEN_ERR);
	
	fprintf(script, "set terminal qt persist size 500,500\n");
	fprintf(script, "plot \"../data/%s\" using 2:3 skip 6 with points\n", inst->fileIn);
	fprintf(script, "set title \"data example\"\n");
	
	fclose(script);
	return 0;
} /*write_plotting_script*/


/*freeing the memory*/

void freeInst(instance* inst) {
	free(inst->pts);
	free(inst->costs);
	free(inst);
} /*freeInst*/