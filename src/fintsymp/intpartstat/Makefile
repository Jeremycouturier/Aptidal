#
# Makefile pour intel
# 

############## variables modifiables ################
include libfintsymp/Makefile

#compilateur fortran
FOR    = ifort 
#FOR    = gfortran
ifeq ($(USE_MPI), 1) 
FOR = mpiifort -DUSE_MPI=1
#FOR = mpif90  -DUSE_MPI=1
else
endif

#format fixe du fortran
FIX =    #-fixed

#options d'optimisation
OPT = 


#options de debuggage
DEBUG= 

#options fournis au compilateur
FDEP= # -gen-dep=dep.make  
FINCLUDE= -I$(LIBFINTSYMP_DIR)/src -module obj/
FOPT =  -O3 -g -traceback  -xhost 
FFLAGS=  $(FINCLUDE) $(FOPT) -no-save-temps -nogen-interfaces -warn all -warn notruncated_source -fpp -DTREAL=8 -fp-model precise 

# pour gfortran
#FFLAGS= -std=f2008 -std=legacy -ffree-line-length-none -DTREAL=8  -Wall  -cpp -O3 -I$(LIBFINTSYMP_DIR)/src -Jobj -Iobj -Wno-unused-dummy-argument

#librairies systemes
#LIBS = -limf
LIBS= -lm 

#options fournis au linker
LDFLAGS =   

VPATH=.:./src/:$(LIBFINTSYMP_VPATH)
#
################# Ne rien modifier apres cette ligne ################
FC=$(FOR)

.SUFFIXES : .f .F90 .o

SRC_FULL12 = ctrl_distpartpla.f \
	filepar.f minmax_mda.f naf_diffang.f ctrl_distpartstar.f intgsympstat.f mainpartstat.f 

OBJS_FULL12 = $(LIBFINTSYMP_OBJS) $(addprefix obj/,$(SRC_FULL12:.f=.o)) 


MODUL = obj/*.mod *.mod
EXE = intpartstat.x

$(EXE): $(OBJS_FULL12)
	$(FOR) $(OBJS_FULL12) $(FFLAGS) $(LIBS) -o $(EXE)


# pour mpi
EXE_MPI = intpartstat_mpi.x

mpi: 
	rm -f libfintsymp/obj/*.o $(OBJS_FULL12) *.i $(MODUL)
	make $(EXE_MPI) USE_MPI=1

$(EXE_MPI) : $(OBJS_FULL12)
	$(FOR) $(OBJS_FULL12) $(FFLAGS) $(LIBS) -o $@


clean:
	rm -f libfintsymp/obj/*.o $(OBJS_FULL12) *.i $(MODUL) $(EXE) $(EXE_MPI)

#regles par defaut pour les .f pour les fichiers sources
obj/%.o : src/%.f 
	$(FOR) $(FFLAGS) -c $< -o $@

### verification de base
check:
	cd tests && ./verif_naf1h.sh
	cd tests && ./verif_minmaxaei1h.sh
	cd tests && ./verif_minmaxalp.sh
	cd tests && ./verif_nafalp.sh
	cd tests && ./verif_nafpart.sh
	cd tests && ./verif_diststarpart.sh
	cd tests && ./verif_distplapart.sh
	cd tests && ./verif_partelln.sh
	cd tests && ./verif_orbpla5.sh

#regles par defaut pour les .f pour la libraririe libfintsymp/src
%.o : %.F90
	$(FOR) $(FFLAGS) -c $< -o $@

 
#############gestion des dependences
ifeq ($(wildcard dep.make),) 
else
 # include dep.make
endif
