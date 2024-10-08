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
FOPT =  -O3 -g -traceback  -axCORE-AVX512
FFLAGS=  $(FINCLUDE) $(FOPT) -no-save-temps -nogen-interfaces -warn all -warn notruncated_source -fpp -DTREAL=8 -fp-model precise 

# pour gfortran
#FFLAGS= -std=f2008 -std=legacy -ffree-line-length-none -DTREAL=8 -Wall  -cpp -O3 -I$(LIBFINTSYMP_DIR)/src -Jobj -Iobj -Wno-unused-dummy-argument

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

OBJS_FULL12 = $(LIBFINTSYMP_OBJS) \
	filepar.o minmax_mda.o minmax_mdc.o naf_diffang.o ctrl_ener.o ctrl_diststar.o naf_difflpm.o \
	minmax_dae2.o ctrl_distpla.o intgsympstat.o intgrestartstat.o mainplastat.o  


MODUL = obj/*.mod *.mod 
EXE = intplastat.x

$(EXE): $(OBJS_FULL12)
	$(FOR) $(OBJS_FULL12) $(FFLAGS) $(LIBS) -o $(EXE)


# pour mpi
EXE_MPI = intplastat_mpi.x

mpi: 
	rm -f libfintsymp/obj/*.o $(OBJS_FULL12) *.i $(MODUL)
	make $(EXE_MPI) USE_MPI=1

$(EXE_MPI) : $(OBJS_FULL12)
	$(FOR) $(OBJS_FULL12) $(FFLAGS) $(LIBS) -o $@


clean:
	rm -f libfintsymp/obj/*.o $(OBJS_FULL12) *.i $(MODUL) $(EXE) $(EXE_MPI)

### verification de base
check:
	$(FOR) tests/verifincl0.f $(LIBFINTSYMP_OBJS) $(FFLAGS) $(LIBS) -o tests/verifincl0.x
	#cd tests && ./verifincl0.x
	cd tests && ./verif_minmaxaei.sh
	cd tests && ./verif_minmaxalp.sh
	cd tests && ./verif_minmaxalpj.sh
	cd tests && ./verif_minmaxalch.sh
	cd tests && ./verif_minmaxalcj.sh
	cd tests && ./verif_minmaxae2h.sh
	cd tests && ./verif_minmaxae2j.sh
	cd tests && ./verif_naf1h.sh
	cd tests && ./verif_naf2h.sh
	cd tests && ./verif_naf3h.sh
	cd tests && ./verif_naf4h.sh
	cd tests && ./verif_naf1j.sh
	cd tests && ./verif_naf2j.sh
	cd tests && ./verif_naf3j.sh
	cd tests && ./verif_naf4j.sh
	cd tests && ./verif_ctrl.sh

#regles par defaut pour les .f pour la libraririe libfintsymp/src
%.o : %.F90
	$(FOR) $(FFLAGS) -c $< -o $@


#############gestion des dependences
ifeq ($(wildcard dep.make),) 
else
 # include dep.make
endif
