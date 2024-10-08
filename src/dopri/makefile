FF = gcc
OPTF = -O0

all: test


dopri.o: dopri.c dopri.h
	$(FF) $(OPTF) -c dopri.c

test: main.c dopri.o
	$(FF) main.c $(OPTF) dopri.o -lm -o test

clean: 
	\rm -f *.o *.log
