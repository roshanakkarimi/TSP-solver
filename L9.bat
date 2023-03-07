REM DA USARE PER OTTENERE UN FILE EXEGUIBILE PARTENDO DA UNO A 9 FILE OGGETTO .o.
REM Linka fino a 9 file oggetto in un file eseguibile, dandogli il nome del primo file oggetto.
REM Digitate il primo file oggetto senza l'estensione (.o), e ogni altro file oggetto con l'estensione (.o).
REM Es: L List_Tester List.o
gcc -o %1.exe %1.o %2 %3 %4 %5 %6 %7 %8 %9
