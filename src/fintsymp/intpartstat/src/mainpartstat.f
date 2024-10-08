!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file mainpartstat.f 
!!  \brief Integration  generale d un systeme planetaire  avec differents integrateurs
!! ayant un ensemble de particule-tests.
!!
!!     Realise des statistiques sur les simulations effectuees des particules
!!
!!
!! version sequentielle ou mpi
!! 
!! la version mpi est activee si le prefined USE_MPI=1
!!
! Derniere revision: 15/09/2014
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
      use mod_ciread_partnewt
      use mod_ciwrite_partnewt
      use mod_intg_symp
      use mod_filepar
      implicit none

      type(t_cidist_mpi) :: distributeur ! distributeur de conditions initiales
      type(t_extfilepar) :: filepar ! fichier de parametres lus + donnees derivees
      type(t_ciread_planewt) :: readerpla ! lecteur des conditions initiales des planetes
      type(t_ciread_partnewt) :: readerpart ! lecteur des conditions initiales des particules
      type(t_ciwrite_planewt) :: writerpla ! ecrivain des conditions initiales des planetes
      type(t_ciwrite_partnewt) :: writerpart ! ecrivain des conditions initiales des particules
      type(t_cidist_integfunc_arg) :: multipledata ! donnee utilisateur : fichier de parametres + autres
      procedure(t_cidist_manyuserfunc), pointer :: pfcninteg
      integer cmd_count
      character(len=512) :: fileparname, arg
      integer rank, nbcipla
      class(t_ci_base), pointer :: cipla !condition initiale du systeme planetaire
      
#if USE_MPI
      integer nbcipart
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
         write(*,*) "e.g. : intpartstat essai.par" 
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

! --- ouvre le fichier de conditions initiales du systeme planetaire
        call readerpla%set_name(filepar%nf_initext, 8)
        nbcipla = readerpla%get_countline()
        if (nbcipla.ne.1) then 
         write(*,*) "Erreur : Un seul systeme planetaire supporte !!!"
         write(*,*) "Nombre de systemes lus : ", nbcipla
         write(*,*) "Arret du programme !!!"
         stop
        endif
        call readerpla%fopen()
        nbcipla = readerpla%freadline(cipla)
        call readerpla%fclose()

! --- ouvre le fichier de conditions initiales des particules
        call readerpart%set_name(filepar%nf_initpart, 9)
        call distributeur%set_manyblocksize(filepar%part_blocksize)

      if (distributeur%m_ismaster)  then
! --  ecrit les conditions lues
        call writerpla%set_name(filepar%nf_ci_pla, 66)
        call writerpla%fwrite(readerpla)
        call writerpart%set_name(filepar%nf_ci_part, 66)
        call writerpart%fwrite(readerpart)
        
#if USE_MPI
         nbcipart = readerpart%get_countline()
         call distributeur%set_nb_ci(nbcipart)
         call distributeur%set_sendblocksize(filepar%part_blocksize)
#endif
       else
          
       endif
       
!desactive pour tous les esclaves sauf  le noeud 1 la sortie pour les planetes 
#if USE_MPI
        multipledata%m_boutpla = .false.
       if (rank.eq.1) then
        multipledata%m_boutpla = .true.       
       endif
#else
        multipledata%m_boutpla = .true.       
#endif

! distribue toutes les conditions initiales  
        call readerpart%fopen()
        multipledata%m_filepar = filepar
        multipledata%m_cipla => cipla
        pfcninteg => cidist_integfunc
        call distributeur%runmany(readerpart, multipledata, pfcninteg)
        call readerpart%fclose()
        
        deallocate(cipla)
        
#if USE_MPI
      call MPI_FINALIZE(code)
#endif
      stop     
      end 

      
