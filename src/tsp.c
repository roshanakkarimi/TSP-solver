#include "tsp.h"
#include <stdbool.h>

/*global*/

double (*pickers[2])(int i, int nearest_prob, int* sol, const instance* inst) = 
{
	greedy_picker,
	grasp_picker
};
	
const char mods[2][50] =
{
	"greedy",
    "GRASP"
};

/*utilities*/

double dist(int i, int j, const point* pts) {
	return sqrt(sq_dist(i, j, pts));
} /*dist*/

double sq_dist(int i, int j, const point* pts) {
	return pow(pts[i].x - pts[j].x, 2) + pow(pts[i].y - pts[j].y, 2);
} /*dist*/

void swapInt(int* i1, int* i2) {
	int temp;
	if(i1 == i2 || (*i2) == (*i1)) return; /*optimization*/
	temp = *i2;
	*i2 = *i1;
	*i1 = temp;
} /*swapInt*/

void swapDouble(double* d1, double* d2) {
	double temp;
	if(d1 == d2 || fabs((*d2) - (*d1)) < EPSILON) return; /*optimization*/
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
	printf("\nFATAL ERROR:\n%s", err);
	return errType;
} /*myError*/

bool checkSol(double z, const int* sol, const instance* inst){
	int i, n = inst->nnodes;
	double z_check = 0;
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
	for(i = 0; i < n - 1; i++)
		z_check += inst->costs[sol[i] * n + sol[i + 1]];
	z_check += inst->costs[sol[i] * n + sol[0]];
	if(fabs(z_check - z) > EPSILON){
		printf("Wrong: z_check = %f while z = %f\n", z_check, z);
		return false;
	} /*if*/
	return true;
} /*checkSol*/


/*solving*/

void updateBest(double z, const int* sol, instance* inst){
	int i;
	if(checkSol(z, sol, inst)){ 
		if(z < inst->zbest || inst->zbest == -1) {
		inst->zbest = z;
		for(i = 0; i < inst->nnodes; i++)
			inst->best_sol[i] = sol[i];
		} /*if*/
	} /*if*/
	/*else printf("Solution is not acceptable!\n");*/
} /*updateBest*/

double greedy_picker(int i, int nearest_prob, int* sol, const instance* inst){ /*just a wrapper to use node_picker standard*/
	return minCost(i, sol, inst);
} /*greedy_picker*/

double grasp_picker(int i, int nearest_prob, int* sol, const instance* inst){
	if(rand() % 100 < nearest_prob)
		return minCost(i, sol, inst);
	else 
		return secMinCost(i, sol, inst);
} /*grasp_picker*/


/* input elaboration*/

void cmdHelp() {
	printf("\n--Available options:--\n -f to set input file (.tsp, not to be written); if notspecified, randomly generated\n");
	printf(" -rs to set random seed\n -test to run in test mode\n -m to set solving algorithm (lowercase)\n");
	printf(" -tl to set overall time limit\n -two to apply two opt. alg. to refine the solution\n");
    printf(" -ns to set number of runs to perform (only with -test opt.)\n -p to set probability\n -help to see options\n");
} /*help*/

void dispPars(const instance* inst) {
	printf("\n--Current parameters:--\n");
	printf(" time limit: %f\n", inst->timelimit);
	printf(" randseed: %d\n", inst->randseed);
	printf(" verbosity level: %d\n", inst->verbosity);
	printf(" input file: %s.tsp\n", inst->fileIn);
	printf(" probability: %d\n", inst->prob);
} /*dispPars*/

void initInst(instance *inst) {
	strcpy(inst->fileIn, "rand");
	inst->verbosity = 1;
	inst->randseed = DEFAULT_RAND;
	inst->timelimit = 300;
	inst->zbest = -1;
	inst->n_sim = 1;
	inst->prob = 100;
	inst->two_opt = false;
	inst->mode = GRASP;
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
				inst->mode = GREEDY;
			else if(strcmp(argv[i], "grasp") == 0)
				inst->mode = GRASP;
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
		if(strcmp(argv[i], "-two") == 0){
			inst->two_opt = true;
			invalid_opt = false;
			continue;
		} /*use two opt. alg.*/
		if(strcmp(argv[i],"-help") == 0) invalid_opt = false;
		else if(invalid_opt) printf("\nINVALID OPTION:");
		cmdHelp();
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
	inst->nnodes = MIN_NODES + rand() % (MAX_NODES - MIN_NODES + 1);
	alloc_inst(inst);
	for(i = 0; i < inst->nnodes; i++){
		inst->pts[i].x = rand() % MAX_COORD;
		inst->pts[i].y = rand() % MAX_COORD;
	} /*for*/
} /*rand_points*/

void compute_costs(instance* inst, cost fc){
	int i, j;
	for(i = 0; i < inst->nnodes; i++)
		for(j = 0; j < inst->nnodes; j++)
			inst->costs[i * inst->nnodes + j] = (i == j) * INFINITE + fc(i, j, inst->pts);
} /*compute_costs*/


/*output elaboration*/

int write_plotting_script(const char* fileOut, int nnodes, int ord){
	FILE *script = fopen("gnuplot_out.p", ((ord == 0) ? "w" : "a"));
	if (script == NULL) return myError("Couldn't create the script!", FILE_OPEN_ERR);
	
	if(ord == 0){
		fprintf(script, "set terminal qt persist size 500,500\n");
		fprintf(script, "set key off\n");
	}
	/*fprintf(script, ", \\\n\"\" skip %d with linespoints lc rgb %d\n", 3 + (3 * ord + 1) * nnodes, 2000 * ord);*/
    /*fprintf(script, "replot \"../out/%s.dat\" skip %d with lines lc rgbcolor %d\n", fileOut, 3 + (3 * ord + 1) * nnodes, 101101 * ord);*/
	if(ord) fprintf(script, "replot \"../out/%s.dat\" index %d with lines lc rgbcolor 16777215 lw 3\n", fileOut, ord);
	fprintf(script, "%splot \"../out/%s.dat\" index 0 using 2:3:1 with labels\n", (ord ? "re" : ""), fileOut);
	fprintf(script, "replot \"../out/%s.dat\" index %d with lines lc rgb 0\n", fileOut, ord + 1);
	fprintf(script, "pause 0.5\n");
	
	fclose(script);
	return 0;
} /*write_plotting_script*/

int write_out_file(const instance* inst, const int* sol, const char* fileOut, const char* mode){
	char fname[100];
	FILE *out;
	int i;
	sprintf(fname, "../out/%s.dat", fileOut);
	out = fopen(fname, mode);
	if(out == NULL) return myError("Couldn't create the output file!", FILE_OPEN_ERR);
	if(strcmp(mode, "w") == 0){
		fprintf(out, "#NODES\n");
		for(i = 0; i < inst->nnodes; i++)
			fprintf(out, "%d %f %f\n", i, inst->pts[i].x, inst->pts[i].y);
	} /*if*/
	fprintf(out, "\n\n\n#EDGES\n");
	for(i = 0; i < inst->nnodes - 1; i++){
		fprintf(out, "%f %f\n", inst->pts[sol[i]].x, inst->pts[sol[i]].y);
		fprintf(out, "%f %f\n\n", inst->pts[sol[i+1]].x, inst->pts[sol[i+1]].y);
	} /*for*/
	fprintf(out, "%f %f\n", inst->pts[sol[i]].x, inst->pts[sol[i]].y); /*edge closing cycle*/
	fprintf(out, "%f %f", inst->pts[sol[0]].x, inst->pts[sol[0]].y);
	
	fclose(out);
	return 0;
} /*write_out_file*/


/*freeing the memory*/

void freeInst(instance* inst) {
	free(inst->pts);
	free(inst->best_sol);
	free(inst->costs);
	free(inst);
} /*freeInst*/