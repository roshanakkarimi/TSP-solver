REM DA USARE PER OTTENERE FINO A 9 FILE OGGETTO PARTENDO DA ALTRETTANTI FILE SORGENTI .c.
REM Precompila fino a 9 sorgenti, Compila ogni risultato mediante il C Ansi, Assembla ogni risultato.
REM Compila mediante il C Ansi da uno a 9 file e li assemblea.
gcc -ansi -pedantic -Wall -Werror -c %1.c
gcc -ansi -pedantic -Wall -Werror -c %2.c
gcc -ansi -pedantic -Wall -Werror -c %3.c
gcc -ansi -pedantic -Wall -Werror -c %4.c
gcc -ansi -pedantic -Wall -Werror -c %5.c
gcc -ansi -pedantic -Wall -Werror -c %6.c
gcc -ansi -pedantic -Wall -Werror -c %7.c
gcc -ansi -pedantic -Wall -Werror -c %8.c
gcc -ansi -pedantic -Wall -Werror -c %9.c
