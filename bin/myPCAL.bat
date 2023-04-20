REM Create executable by compiling and linking one single file and tsp.o.
gcc -ansi -pedantic -Wall -Werror -c ../src/tsp.c
gcc -ansi -pedantic -Wall -Werror -c ../src/%1.c
gcc -O3 -o %1.exe %1.o tsp.o