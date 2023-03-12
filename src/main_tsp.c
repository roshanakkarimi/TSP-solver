#include "tsp.h"

int main(int argc, char **argv)
{
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	initInst(inst);
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
    write_plotting_script(inst);
    system("gnuplot gnuplot_out.p");
    
    /*plot {dat file} with points, {dat file} with linespoints*/

	free(inst);
	return 0;
}