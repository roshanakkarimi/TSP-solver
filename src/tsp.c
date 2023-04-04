#include "tsp.h"
#include <stdbool.h>

/*greedy heur. modes*/
const char heur_modes[2][50] =
{
	"greedy",
    "GRASP"
};

/*ref. heur. and metaheur.*/
const char ref_modes[2][50] =
{
	"two opt.",
    "two opt. plus tabu"
};




/*utilities*/

double dist(int i, int j, const point* pts) {
	return sqrt(sq_dist(i, j, pts));
} /*dist*/

double sq_dist(int i, int j, const point* pts) {
	return (pts[i].x - pts[j].x) * (pts[i].x - pts[j].x) + (pts[i].y - pts[j].y) * (pts[i].y - pts[j].y);
} /*dist*/

void swapInt(int* i1, int* i2) {
	int temp;
	temp = *i2;
	*i2 = *i1;
	*i1 = temp;
} /*swapInt*/

void swapDouble(double* d1, double* d2) {
	double temp;
	temp = *d2;
	*d2 = *d1;
	*d1 = temp;
} /*swapInt*/

/*returns minCost of still possible edges involving node sol[i].
Puts the other node of selected edge in sol[i+1] by a swapping operation.*/
double minCost(int i, int* sol, const instance* inst) {
	int j, indMin;
	double dj;
	int last = sol[i];
	double minCost = inst->costs[last * inst->nnodes + sol[i + 1]];
	indMin = i + 1;
	for(j = i + 2; j < inst->nnodes; j++)
		if((dj = inst->costs[last * inst->nnodes + sol[j]]) < minCost){
			indMin = j;
			minCost = dj;
		} /*if*/
	swapInt(sol + i + 1, sol + indMin);
	return minCost;
} /*minCost*/

double secMinCost(int i, int* sol, const instance* inst){
	int j, iMin, iSecMin;
	double dj;
	double* costs_from_i = inst->costs + inst->nnodes * sol[i];
	double minCost = costs_from_i[sol[i + 1]];
	double secMinCost = costs_from_i[sol[i + 2]];
	iMin = i + 1;
	iSecMin = i + 2;
	if(secMinCost < minCost){
		swapInt(&iMin, &iSecMin);
		swapDouble(&minCost, &secMinCost);
	} /*if*/
	for(j = i + 3; j < inst->nnodes; j++)
		if((dj = costs_from_i[sol[j]]) < minCost){
			iSecMin = iMin;
			secMinCost = minCost;
			iMin = j;
			minCost = dj;
		} /*if*/
		else if(dj < secMinCost){
			iSecMin = j;
			secMinCost = dj;
		} /*if*/
	swapInt(sol + i + 1, sol + iSecMin);
	return secMinCost;
} /*secMinCost*/


/*managing errors and debug*/

int myError(const char* err, int errType){
	if(errType <= 10)
		printf("\nFATAL ERROR:\n%s", err);
	else
		printf("\nWARNING:\n%s", err);
	return errType;
} /*myError*/

double getSolCost(const int* sol, const instance* inst){
	int i;
	double z = 0;
	for(i = 0; i < inst->nnodes - 1; i++)
		z += inst->costs[sol[i] * inst->nnodes + sol[i + 1]];
	z += inst->costs[sol[i] * inst->nnodes + sol[0]];
	return z;
} /*getSolCost*/

bool checkSol(double z, const int* sol, const instance* inst){
	int i, n = inst->nnodes;
	double z_check;
	int* count = calloc(n, sizeof(int));
	assert(count != NULL);
	for(i = 0; i < n; i++)
		count[sol[i]]++;
	for(i = 0; i < n; i++)
		if(count[i] != 1){
			printf("Wrong: node %d counted %d times\n", i , count[i]);
			free(count);
			return false;
		} /*if*/
	free(count);
	z_check = getSolCost(sol, inst);
	if(fabs(z_check - z) > EPSILON){
		printf("Wrong: z_check = %.2f while z = %.2f so diff = %.2f\n", z_check, z, z_check - z);
		return false;
	} /*if*/
	return true;
} /*checkSol*/


/*solving*/

void updateBest(double z, const int* sol, instance* inst){
	int i;
	if(checkSol(z, sol, inst)){ 
		if(z < inst->zbest) {
		inst->zbest = z;
		for(i = 0; i < inst->nnodes; i++)
			inst->best_sol[i] = sol[i];
		} /*if*/
	} /*if*/
} /*updateBest*/

/* input elaboration*/

void cmdHelp() {
	int i;
	printf("\n--Available options:--\n -f to set input file (.tsp, not to be written); if not specified, randomly generated\n");
	printf(" -n to set number of nodes to be randomly generated\n");
	printf(" -rs to set random seed\n -p to set probability\n -tl to set overall time limit\n");
	printf(" -m to set solving algorithm (lowercase)\n");
	printf("    options: ");
	for(i = 0; i < N_HEUR; i++) printf("%s * ", heur_modes[i]);
	printf("\n -r to set refinement and metaheur. algorithm (lowercase)\n");
	printf("    options: ");
	for(i = 0; i < N_REF; i++) printf("%s * ", ref_modes[i]);
	/*printf("\n -two to apply two opt. alg. to refine the solution\n");
	printf(" -tabu to use tabu alg. to escape local mins. (with -two option)\n -tt to set tabu tenure(with -tabu option)\n");*/
    printf("\n -test to run in test mode\n -ns to set number of runs to perform (with -test opt.)\n");
	printf(" -help to see options\n");
} /*help*/

void dispPars(const instance* inst) {
	printf("\n--Current parameters:--\n");
	printf(" time limit: %f\n", inst->timelimit);
	printf(" randseed: %d\n", inst->randseed);
	printf(" verbosity level: %d\n", inst->verbosity);
	printf(" input file: %s.tsp\n", inst->fileIn);
	printf(" probability: %d\n", inst->prob);
	printf(" tabu tenure: %d\n", inst->tabu_tenure);
} /*dispPars*/

void initInst(instance *inst) {
	strcpy(inst->fileIn, "\0");
	inst->verbosity = 1;
	inst->nnodes = N_DEF_NODES;
	inst->randseed = DEFAULT_RAND;
	inst->timelimit = DEFAULT_TL;
	inst->zbest = INFINITE_DBL;
	inst->n_sim = 1;
	inst->prob = 100;
	inst->heur_mode = GRASP;
	inst->ref_mode = NOTHING;
	inst->tabu_tenure = DEFAULT_TT;
} /*initInst*/

bool parse_cmd(int argc, char** argv, instance *inst) {
	int i = 0;
	bool invalid_opt;
	bool test = false;
	while(i < argc - 1) {
		i++;
		invalid_opt = true;
		if(strcmp(argv[i],"-f") == 0){
			strcpy(inst->fileIn, argv[++i]);
			invalid_opt = false;
			continue;
		} /*set input file*/
		if(strcmp(argv[i],"-m") == 0){
			if(strcmp(argv[++i], "greedy") == 0)
				inst->heur_mode = GREEDY;
			else if(strcmp(argv[i], "grasp") == 0)
				inst->heur_mode = GRASP;
			invalid_opt = false;
			continue;
		} /*set algorithm to solve inst.*/
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
		if(strcmp(argv[i],"-n") == 0){
			inst->nnodes = abs(atoi(argv[++i])); 
			invalid_opt = false;
			continue;
		} /*set number of random nodes*/
		if(strcmp(argv[i],"-tt") == 0){
			inst->tabu_tenure = abs(atoi(argv[++i])); 
			invalid_opt = false;
			continue;
		} /*set tabu tenure*/
		if(strcmp(argv[i],"-ns") == 0){
			inst->n_sim = abs(atoi(argv[++i])); 
			invalid_opt = false;
			continue;
		} /*set number of simulations*/
		if(strcmp(argv[i],"-p") == 0){
			inst->prob = abs(atoi(argv[++i])); 
			invalid_opt = false;
			continue;
		} /*set probability*/
		if(strcmp(argv[i],"-v") == 0){
			inst->verbosity = atoi(argv[++i]);
			invalid_opt = false;
			continue;
		} /*set verbosity param.*/
		if(strcmp(argv[i], "-test") == 0){
			test = true;
			invalid_opt = false;
			continue;
		} /*test mod.*/
		if(strcmp(argv[i],"-r") == 0){
			if(strcmp(argv[++i], "two") == 0)
				inst->ref_mode = TWO;
			else if(strcmp(argv[i], "tabu") == 0)
				inst->ref_mode = TWO_TABU;
			invalid_opt = false;
			continue;
		} /*set algorithm to refine sol.*/
		if(strcmp(argv[i],"-help") == 0) invalid_opt = false;
		else if(invalid_opt) printf("\nINVALID OPTION:");
		cmdHelp();
		exit(0);
	} /*while*/
	return test;
} /*parse_cmd*/

int read_fileIn(instance* inst) {
	
	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	bool node_sect = false;
	
	FILE* fileIn;
    sprintf(line, "../data/%s.tsp", inst->fileIn);
	fileIn = fopen(line, "r");
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
			alloc_inst(inst);
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

void alloc_inst(instance* inst){
	inst->pts = malloc(inst->nnodes * sizeof(point));
	assert(inst->pts != NULL);
	inst->best_sol = malloc(inst->nnodes * sizeof(int));
	assert(inst->best_sol != NULL);
	inst->costs = malloc(inst->nnodes * inst->nnodes * sizeof(double));
	assert(inst->costs != NULL);
} /*alloc_inst*/

void rand_points(instance* inst) {
	int i;
	for(i = 0; i < inst->nnodes; i++){
		inst->pts[i].x = rand() % MAX_COORD;
		inst->pts[i].y = rand() % MAX_COORD;
	} /*for*/
} /*rand_points*/

void compute_costs(instance* inst, cost fc){
	int i, j;
	for(i = 0; i < inst->nnodes; i++){
		for(j = 0; j < i; j++)
			inst->costs[i * inst->nnodes + j] = fc(i, j, inst->pts);
		inst->costs[i * inst->nnodes + i] = INFINITE_DBL;
	} /*for*/
	for(i = 0; i < inst->nnodes; i++){
		for(j = i + 1; j < inst->nnodes; j++)
			inst->costs[i * inst->nnodes + j] = inst->costs[j * inst->nnodes + i];
	} /*for*/
} /*compute_costs*/


/*output elaboration*/

int update_plotting_script(FILE* script, const char* fileIn, int ord){
	char fileOut[50];
	sprintf(fileOut, "../out/%s.dat", fileIn);
	if(ord){
		fprintf(script, "replot \"%s\" index %d with lines lc rgbcolor 16777215 lw 3\n", fileOut, ord);
		fprintf(script, "replot \"%s\" index 0 using 2:3:1 with labels\n", fileOut);
	} /*if*/
	fprintf(script, "replot \"%s\" index %d with lines lc rgb 0\n", fileOut, ord + 1);
	fprintf(script, "pause 0.5\n");
	return 0;
} /*update_plotting_script*/

int update_out_file(const instance* inst, const int* sol, FILE *out){
	int i;
	fprintf(out, "\n\n\n#EDGES\n");
	for(i = 0; i < inst->nnodes - 1; i++){
		fprintf(out, "%f %f\n", inst->pts[sol[i]].x, inst->pts[sol[i]].y);
		fprintf(out, "%f %f\n\n", inst->pts[sol[i+1]].x, inst->pts[sol[i+1]].y);
	} /*for*/
	fprintf(out, "%f %f\n", inst->pts[sol[i]].x, inst->pts[sol[i]].y); /*edge closing cycle*/
	fprintf(out, "%f %f", inst->pts[sol[0]].x, inst->pts[sol[0]].y);

	return 0;
} /*update_out_file*/


/*freeing the memory*/

void freeInst(instance* inst) {
	free(inst->pts);
	free(inst->best_sol);
	free(inst->costs);
	free(inst);
} /*freeInst*/