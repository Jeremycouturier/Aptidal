#
# Makefile pour intel
# 

############## variables modifiables ################
include libfintsymp/Makefile

#compilateur fortran
#FOR    = ifort 
FOR    = gfortran
#FOR = mpiifort -DUSE_MPI=1

#format fixe du fortran
FIX =    #-fixed

#options d'optimisation
OPT = 


#options de debuggage
DEBUG= 

#options fournis au compilateur
# flags pour debuggage
#FFLAGS= -fp-model precise  -g -O0
# flags pour optimisation
FINCLUDE= -I$(LIBFINTSYMP_DIR)/src -module obj/
FOPT =  -O3 -g -traceback  -xhost 
# FFLAGS=  $(FINCLUDE) $(FOPT) -no-save-temps -nogen-interfaces -warn all -warn notruncated_source -fpp -DTREAL=8 -fp-model precise 

# pour gfortran
FFLAGS= -std=f2008 -std=legacy -ffree-line-length-none -DTREAL=8  -cpp -O3 -I$(LIBFINTSYMP_DIR)/src
# flags pour home 
#FFLAGS= -O3

#options fournis au compilateur pour profile
PROFLAGS= 

#librairies systemes
#LIBS = -limf
LIBS= -lm 

#options fournis au linker
LDFLAGS =   

VPATH=.:./src/:$(LIBFINTSYMP_VPATH)
#
################# Ne rien modifier apres cette ligne ################
FC=$(FOR)

.SUFFIXES : .f .o

OBJS_FULL12 = $(LIBFINTSYMP_OBJS) \
	filepar.o intgsymp2014.o intplabasic.o 


MODUL = *.mod 
EXE = intplabasic.x

$(EXE): $(OBJS_FULL12)
	$(FOR) $(OBJS_FULL12) $(FFLAGS) $(LIBS) -o $(EXE)


clean:
	rm -f libfintsymp/obj/*.o $(OBJS_FULL12) *.i $(MODUL) $(EXE)

