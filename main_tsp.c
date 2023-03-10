#include "tsp.h"

int main(int argc, char **argv)
{
	instance* inst = malloc(sizeof(instance));
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
	greedy_solution(inst);
	free(inst);
	return 0;
}
