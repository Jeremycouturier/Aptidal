!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file mainplastat.f 
!!  \brief Integration  generale d un systeme planetaire  avec differents integrateurs
!!
!!     Realise des statistiques sur les simulations effectuees
!!
!!
!! version sequentielle ou mpi
!! 
!! la version mpi est activee si le prefined USE_MPI=1
!!
! Derniere revision: 24/03/2014
!***********************************************************************
#include "realc.f"
            
!*******************************************************
!> Programme principal  
!*******************************************************
      program mainplastat
#if USE_MPI
      use mpi
#endif
      use mod_cidist_mpi
      use mod_ciread_planewt
      use mod_ciwrite_planewt
      use mod_ciread_pladistmin
      use mod_ciwrite_pladistmin
      use mod_ciread_multibase
      use mod_ciwrite_multibase
      use mod_intg_symp
      use mod_filepar
      use mod_intg_symp_start
      use mod_intg_symp_restart
      implicit none

      type(t_cidist_mpi) :: distributeur ! distributeur de conditions initiales
      type(t_extfilepar) :: filepar ! fichier de parametres lus + donnees derivees
      type(t_ciread_multibase) :: reader ! lecteur des conditions initiales
      type(t_ciwrite_multibase) :: writer ! ecrivain des conditions initiales
      type(t_ciread_planewt), target :: readernewt ! lecteur des conditions initiales
      type(t_ciwrite_planewt), target :: writernewt ! ecrivain des conditions initiales
      type(t_ciread_pladistmin), target :: readerdistmin ! lecteur des distances minimales entre planetes
      type(t_ciwrite_pladistmin), target :: writerdistmin ! ecrivain des distances minimales entre planetes
      procedure(t_cidist_userfunc), pointer :: pfcninteg
      type(t_ciread_base_ptr), dimension(:), allocatable :: arreader
      type(t_ciwrite_base_ptr), dimension(:), allocatable :: arwriter
      integer cmd_count
      character(len=512) :: fileparname, arg
      integer rank
      
#if USE_MPI
      integer nbci
      integer code
      call MPI_INIT(code)

      call distributeur%init()
      rank = distributeur%get_rank()
#else
      rank = 0      
#endif
      ! recupere le nombre d'arguments
      cmd_count = COMMAND_ARGUMENT_COUNT()
      if(cmd_count.ne.1) then 
         write(*,*) 
         write(*,*) "Usage : <nom du programme> <fich.par>" 
         write(*,*) "e.g. : intplastat essai.par" 
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
       call filepar%fwritepartrip();
      endif
      
! --- create the filename of output files  
      call filepar%create_filename(rank)
      call filepar%open_file()

!---- cas du redemarrage
      if (filepar%if_restart.gt.0) then
  
      
! --- ouvre le(s) fichier(s) de conditions initiales
        if (filepar%ctrl_distpla%m_compute.ne.0) then
          call readerdistmin%set_name(filepar%ctrl_distpla%m_nfdistmin,   &
     &     9)
        endif
        call cidist_integfunc_restart(filepar, readerdistmin)

      else
!---- cas du demarrage standard      

! --- ouvre le(s) fichier(s) de conditions initiales
        if (filepar%ctrl_distpla%m_compute.ne.0) then
          allocate(arreader(1:2))
          call readerdistmin%set_name(filepar%ctrl_distpla%m_nfdistmin,   &
     &     9)
          arreader(2)%m_ptr => readerdistmin
        else
          allocate(arreader(1:1))
        endif
        call readernewt%set_name(filepar%nf_initext, 8)
        arreader(1)%m_ptr => readernewt
        call reader%set_readers(arreader)
        
        

      if (distributeur%m_ismaster)  then
! --  ecrit les conditions lues
        if (filepar%ctrl_distpla%m_compute.ne.0) then
          allocate(arwriter(1:2))
          call writerdistmin%set_name(filepar%nf_cictrlpladistmin, 67)
          arwriter(2)%m_ptr => writerdistmin
        else
          allocate(arwriter(1:1))
        endif
        call writernewt%set_name(filepar%nf_ci, 66)
        arwriter(1)%m_ptr => writernewt
        call writer%set_writers(arwriter)
        call writer%fwrite(reader)
        
#if USE_MPI
         nbci = reader%get_countline()
         call distributeur%set_nb_ci(nbci)
         call distributeur%set_sendblocksize(1)
#endif
       endif
       
! distribue toutes les conditions initiales  
        call reader%fopen()
        pfcninteg => cidist_integfunc
        call distributeur%run(reader, filepar, pfcninteg)
        call reader%fclose()
        
      endif  
        
#if USE_MPI
      call MPI_FINALIZE(code)
#endif
      stop     
      end 

      
