#include "realc.f"
!********************************************************
!
!               Integration  generale d un seul systeme planetaire
!                avec des particules avec differents integrateurs
! version sequentielle ou mpi
! 
! la version mpi est activee si le prefined USE_MPI=1
!
! derive de intplabasic
!
! Derniere revision: 16/06/2014
!********************************************************

            
!*******************************************************
! Programme principal  
!*******************************************************
      program mainbasicpart2014
#if USE_MPI
      use mpi
#endif
      use mod_cidist_mpi
      use mod_ciread_planewt
      use mod_ciread_partnewt
      use mod_intg_symp
      use mod_filepar
      implicit none

      type(t_cidist_mpi) :: distributeur ! distributeur de conditions initiales
      type(t_extfilepar) :: filepar ! fichier de parametres lus + donnees derivees
      type(t_ciread_planewt) :: readerpla ! lecteur des conditions initiales des planetes
      type(t_ciread_partnewt) :: readerpart ! lecteur des conditions initiales des particules
      class(t_ci_base), pointer :: cipla !condition initiale du systeme planetaire
      type(t_cidist_integfunc_arg) :: multipledata ! donnee utilisateur : fichier de parametres + autres
      procedure(t_cidist_manyuserfunc), pointer :: pfcninteg
      
      integer nbcipla

      integer cmd_count
      character(len=512) :: fileparname, arg
      
#if USE_MPI
      integer nbcipart
      integer code
      call MPI_INIT(code)
#endif
      ! recupere le nombre d'arguments
      cmd_count = COMMAND_ARGUMENT_COUNT()
      if(cmd_count.ne.1) then 
         write(*,*) 
         write(*,*) "Usage : <nom du programme> <fich.par>" 
         write(*,*) "e.g. :intpartbasic essai.par" 
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
      
! --- lit le systeme planetaire       
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

#if USE_MPI
       if (distributeur%m_ismaster)  then
         nbcipart = readerpart%get_countline()
         call distributeur%set_nb_ci(nbcipart)
         call distributeur%set_sendblocksize(filepar%part_blocksize)
        endif
#endif
        call readerpart%fopen()
        multipledata%m_filepar = filepar
        multipledata%m_cipla => cipla
        pfcninteg => cidist_integfunc
        call distributeur%runmany(readerpart, multipledata, pfcninteg)
        call readerpart%fclose()
        
#if USE_MPI
      call MPI_FINALIZE(code)
#endif
      stop     
      end 

      
