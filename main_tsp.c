#include "tsp.h"

int main(int argc, char **argv)
{
	instance inst;
	parse_cmd(argc, argv, &inst);
	read_fileIn(&inst);
	printf("Point 48: %f, %f", inst.pts[47].x, inst.pts[47].y);
	return 0;
}