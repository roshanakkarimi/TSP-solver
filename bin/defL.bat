REM Create main_tsp executable by linking every listed file.
REM Specify every .o extension.
gcc -O3 -o main_tsp.exe main_tsp.o tsp.o %1 %2 %3
