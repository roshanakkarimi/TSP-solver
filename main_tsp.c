#include "tsp.h"

int main(int argc, char **argv)
{
	instance* inst = malloc(sizeof(instance));
	assert(inst != NULL);
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
    write_plotting_script(inst);
    system("gnuplot gnuplot_out.p");
    
    /*set terminal qt size 500,500
    plot {dat file} with points, {dat file} with linespoints
    pause -1*/

	free(inst);
	return 0;
}