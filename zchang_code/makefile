SOUR=analysis_counts.c
OBJ=analysis_counts.o
PROG=analysis_counts.exe

$(PROG) : $(OBJ)
	gcc -o $(PROG) -Wg -lm $(OBJ) ./myscaler_funs.o
$(OBJ): $(SOUR)
	gcc -c $(SOUR) myscaler_funs.c
