#include "realc.f"
!********************************************************
!
!               Integration  generale d un systeme
!               planetaire  avec differents integrateurs
! version sequentielle ou mpi
! 
! la version mpi est activee si le prefined USE_MPI=1
!
! derive de integ2010.f by J. Laskar and P. Robutel
!
! Derniere revision: 17/03/2014
! M. GASTINEAU
!  19/10/2012 : Reecriture complete du code pour une meilleure genericite 
!********************************************************

            
!*******************************************************
! Programme principal  
!*******************************************************
      program maininteg2014
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
         write(*,*) "e.g. :intplabasic essai.par" 
        stop
      endif  
      call GET_COMMAND_ARGUMENT(1, arg)
      fileparname = trim(arg)
!      call getparfilename('integ2010.par', filepar)
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

      
