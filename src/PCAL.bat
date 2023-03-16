REM DA USARE PER OTTENERE UN FIKE ESEGUIBILE PARTENDO DA UN UNICO SORGENTE .c.
REM Precompila un sorgente, Compila il risultato mediante il C Ansi, Assembla il risultato, Linka l'oggetto risultato.
gcc -ansi -pedantic -Wall -Werror -o %1.exe %1.c
