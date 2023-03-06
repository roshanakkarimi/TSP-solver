#include "tsp.h"
#include <stdbool.h>

/*utilities*/

double ptdist(point* pt1, point* pt2) {
	return sqrt(pow(pt1->x - pt2->x, 2) + pow(pt1->y - pt2->y, 2));
} /*ptdist*/

double dist(int i, int j, instance* inst) {
	return ptdist(&inst->pts[i], &inst->pts[j]);
} /*dist*/


/* input elaboration*/

void cmdHelp() {
	printf("\n--Available options:--\n -f to set input file \n -tl to set overall time limit\n -rs to set random seed\n -help to see options\n");
} /*help*/

void dispPars(instance* inst) {
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

int main(int argc, char **argv)
{
	instance inst;
	parse_cmd(argc, argv, &inst);
	return 0;
}