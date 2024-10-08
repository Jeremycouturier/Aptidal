#include "realc.f"
!********************************************************
!
!               Integration  generale d un systeme avec 2 etoiles
!               et n-1 planetes  avec differents integrateurs
! version sequentielle ou mpi
! 
! la version mpi est activee si le prefined USE_MPI=1
!
! Derniere revision: 19/11/2015
! M. GASTINEAU
!********************************************************

            
!*******************************************************
! Programme principal  
!*******************************************************
      program maincircumbin2015
#if USE_MPI
      use mpi
#endif
      use mod_cidist_mpi
      use mod_ciread_planewt
      use mod_intg_symp
      use mod_filepar
      implicit none

      type(t_cidist_mpi) :: distributeur ! distributeur de conditions initiales
      type(t_extfilepar) :: filepar ! fichier de parametres lus + donnees derivees
      type(t_ciread_planewt) :: reader ! lecteur des conditions initiales
      procedure(t_cidist_userfunc), pointer :: pfcninteg
      
      integer cmd_count
      character(len=512) :: fileparname, arg
      
#if USE_MPI
      integer nbci
      integer code
      call MPI_INIT(code)
#endif
      ! recupere le nombre d'arguments
      cmd_count = COMMAND_ARGUMENT_COUNT()
      if(cmd_count.ne.1) then 
         write(*,*) 
         write(*,*) "Usage : <nom du programme> <fich.par>" 
         write(*,*) "e.g. : intcircumbinbasic.x essai.par" 
        stop
      endif  
      call GET_COMMAND_ARGUMENT(1, arg)
      fileparname = trim(arg)
      write (*,*) "utilisation de  ", trim(fileparname)


! ---  Lecture et ecriture des parametres      
      call filepar%fread(fileparname)
      
! --- Ecriture des parametres
      if (distributeur%m_ismaster)  then
       call filepar%fwrite(6)
       call filepar%fwritepar()
      endif
      
! --- ouvre le fichier de conditions initiales
        call reader%set_name(filepar%nf_initext, 8)
#if USE_MPI
       if (distributeur%m_ismaster)  then
         nbci = reader%get_countline()
         call distributeur%set_nb_ci(nbci)
         call distributeur%set_sendblocksize(1)
        endif
#endif
        call reader%fopen()
        pfcninteg => cidist_integfunc
        call distributeur%run(reader, filepar, pfcninteg)
        call reader%fclose()
        
#if USE_MPI
      call MPI_FINALIZE(code)
#endif
      stop     
      end 

      