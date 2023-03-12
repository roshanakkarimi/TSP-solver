REM DA USARE PER OTTENERE UN FILE OGGETTO PARTENDO DA UN UNICO SORGENTE .c.
REM Precompila un sorgente, Compila il risultato mediante il C Ansi, Assembla il risultato.
gcc -ansi -pedantic -Wall -Werror -c %1.c
