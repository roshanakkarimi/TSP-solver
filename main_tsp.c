#include "tsp.h"

int main(int argc, char **argv)
{
	instance* inst = malloc(sizeof(instance));
	parse_cmd(argc, argv, inst);
	read_fileIn(inst);
  
    system("gnuplot gnuplot_out.p");
    
    /*set terminal qt size 500,500
    plot {dat file} with points, {dat file} with linespoints
    pause -1*/

	free(inst);
	return 0;
}