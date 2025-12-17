FF = gcc
OPTF = -O2

all: aptidal


structure.o: structure.c structure.h parameters.h transformation.h coefficients.h calculus.h intpla.h
	$(FF) $(OPTF) -c structure.c
	
coefficients.o: coefficients.c coefficients.h structure.o parameters.h
	$(FF) $(OPTF) -c coefficients.c
	
transformation.o: transformation.c transformation.h structure.o coefficients.o parameters.h
	$(FF) $(OPTF) -c transformation.c
	
calculus.o: calculus.c calculus.h structure.o coefficients.o transformation.o parameters.h intpla.h
	$(FF) $(OPTF) -c calculus.c
	
intpla.o: intpla.c intpla.h structure.o coefficients.o transformation.o calculus.o parameters.h
	$(FF) $(OPTF) -c intpla.c

aptidal: main.c structure.o coefficients.o transformation.o calculus.o intpla.o parameters.h
	#$(FF) main.c $(OPTF) structure.o coefficients.o transformation.o calculus.o intpla.o -lgsl -lgslcblas -lm -o aptidal
	$(FF) main.c $(OPTF) structure.o coefficients.o transformation.o calculus.o intpla.o -lm -o aptidal

clean: 
	\rm -f *.o *.log
