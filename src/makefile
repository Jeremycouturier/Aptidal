FF = gcc
OPTF = -O2 

all: aptidal


structure.o: structure.c structure.h parameters.h transformation.h coefficients.h
	$(FF) $(OPTF) -c structure.c
	
coefficients.o: coefficients.c coefficients.h structure.o parameters.h
	$(FF) $(OPTF) -c coefficients.c
	
transformation.o: transformation.c transformation.h structure.o coefficients.o parameters.h
	$(FF) $(OPTF) -c transformation.c
	
calculus.o: transformation.c transformation.h calculus.h structure.o coefficients.o transformation.o parameters.h
	$(FF) $(OPTF) -c calculus.c

aptidal: main.c structure.o coefficients.o transformation.o calculus.o parameters.h
	$(FF) main.c $(OPTF) structure.o coefficients.o transformation.o calculus.o -lgsl -lgslcblas -lm -o aptidal

clean: 
	\rm -f *.o *.log
