
########################################################################
#pour utiliser ce Makefile :
# faire dans le Makefile principal:
# 1. definir TREAL=8 pour double
#            TREAL=10 pour etendu (uniquement gcc)
#            TREAL=16 pour quadruple
#
#
# 3. ajouter : 
#     "include libfintsymp/Makefile" apres ces define
#     "$(LIBFINTSYMP_OBJS)" dans la liste des objets
#     "$(LIBFINTSYMP_VPATH)" dans la liste de VPATH
########################################################################


########################################################################
#exemple :
#
# include libfintsymp/Makefile
# FFLAGS= -I$(LIBFINTSYMP_VPATH) 
# 
# OBJS=main.o   $(LIBFINTSYMP_OBJS)
#
# VPATH=$(LIBFINTSYMP_VPATH):src 
#
#
########################################################################


#nom du dossier de libasdctools
LIBFINTSYMP_DIR=libfintsymp

LIBASDFTOOLS_SRCSF = restart.f kahansum.f arret.f kepsaut.f  kepsaut_hyper.F90 elliptiD3.f coord1.1.f \
  modnaffN.f ode.f dopri8t.f odex1Dt.f odex2Dt.f

# fichiers pour les particules
LIBFINTSYMP_PART_SRCSF= coordpart.f elliptiDpart.f ci_partnewt.f ciread_partnewt.f converterpart.f \
  syspart.f sysparth.f io_txt_carpart.f  io_txt_ellpart.f ciwrite_partnewt.f

LIBFINTSYMP_PART_OBJS=  $(addprefix $(LIBFINTSYMP_DIR)/obj/,$(LIBFINTSYMP_PART_SRCSF:.f=.o)) 

# fichiers pour les circumbinaires
LIBFINTSYMP_CIRCUMBIN_SRCSF= coordcircumbin.f syscircumnewt.f convertercircumbin.f 

LIBFINTSYMP_CIRCUMBIN_OBJS=  $(addprefix $(LIBFINTSYMP_DIR)/obj/,$(LIBFINTSYMP_CIRCUMBIN_SRCSF:.f=.o)) 

#nom des fichiers de src
LIBFINTSYMP_CMN_SRCSF=$(LIBASDFTOOLS_SRCSF) gmsun.f buffer.f io_dump.F90 syssympbase.F90 syspla.f sysplanewth.f \
 schema.f intsympbase.f io_int.f converter.f io_txt_ncomp.f sysplanewtj.f sysplanewtj_adapt.f sysplanewth_adapt.f ci_base.f \
 ci_planewt.f ciread_base.f ciread_planewt.f ciwrite_base.f ciwrite_planewt.f  cidist_seq.f cidist_mpi.f \
 moncin2.f minmax.f naf_base.f io_txt_car.f finalconsumer.f sysextension.f  sysplaexth.f \
 io_txt_mligncomp.f sysplaode.f intodebase.f  forplaexth.f readfile.f interherm.f interdriver.f forplatab.f \
 ctrl_base.f ci_pladistmin.f ciread_pladistmin.f ciwrite_pladistmin.f ci_multibase.f ciread_multibase.f ciwrite_multibase.f

LIBFINTSYMP_CMN_OBJS= $(addprefix $(LIBFINTSYMP_DIR)/obj/,$(patsubst %.F90,%.o, $(LIBFINTSYMP_CMN_SRCSF:.f=.o))) 


#ensemble de tous les objets
LIBFINTSYMP_OBJS= $(LIBFINTSYMP_CMN_OBJS) $(LIBFINTSYMP_PART_OBJS) $(LIBFINTSYMP_CIRCUMBIN_OBJS)



#chemin de recherche pour VPATH
LIBFINTSYMP_VPATH=$(LIBFINTSYMP_DIR)/src:$(LIBFINTSYMP_DIR)/src_part:$(LIBFINTSYMP_DIR)/src_circumbin

############################
#regles par defaut pour les .f pour la libraririe libfintsymp/src
$(LIBFINTSYMP_DIR)/obj/%.o : $(LIBFINTSYMP_DIR)/src/%.f 
	$(FOR) $(FFLAGS) -c $< -o $@

#regles par defaut pour les .f pour la libraririe libfintsymp/src
$(LIBFINTSYMP_DIR)/obj/%.o : $(LIBFINTSYMP_DIR)/src/%.F90
	$(FOR) $(FFLAGS) -c $< -o $@

#regles par defaut pour les .f pour la libraririe libfintsymp/src_part
$(LIBFINTSYMP_DIR)/obj/%.o : $(LIBFINTSYMP_DIR)/src_part/%.f 
	$(FOR) $(FFLAGS) -c $< -o $@

#regles par defaut pour les .f pour la libraririe libfintsymp/src_circum
$(LIBFINTSYMP_DIR)/obj/%.o : $(LIBFINTSYMP_DIR)/src_circumbin/%.f 
	$(FOR) $(FFLAGS) -c $< -o $@
