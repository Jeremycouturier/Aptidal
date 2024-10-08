!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file filepar.f 
!!  \brief gestion du fichier de parametres 
!! les donnees dans le fichier de conditions initiales nf_initext doivent etre en UA et an
!!
!!
! history : creation 17/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour gerer le fichier de parametres
!***********************************************************************
      module mod_filepar

!***********************************************************************
!> @class t_filepar
!! classe de gestion du fichier de parametres : 
!! contient uniquement les donnees du fichier de parametres
!!  
!***********************************************************************
      type  :: t_filepar

        character(100)  :: chemin,nf_rad,int_type
        character(100) :: nf_initext
        integer :: if_int,if_ell,if_car
        integer(8) :: n_iter,n_out
        integer ::    out_ell,if_invar
        real(TREAL) :: tinit,dt,dmax
        integer :: ref_gmsun
        integer :: type_pas
        

       contains
       
       procedure ::fread => filepar_fread ! lecture du fichier de parametres
       procedure ::fwrite => filepar_fwrite ! ecriture du fichier de parametres avec son numero 
       procedure ::fwritepar => filepar_fwritepar ! ecriture du fichier de parametres dont le nom est genere
                              
      end type t_filepar  

!***********************************************************************
!> @class t_extfilepar
!! classe de gestion du fichier parametres : 
!! contient les donnees du fichier de parametres et les variables derivees
!!  
!***********************************************************************
      type, extends(t_filepar)  :: t_extfilepar

        character(512) :: nf_int !< nom du fichier de sortie des integrales premieres
        character(512) :: nf_car !< nom du fichier de sortie des elements cartesiens
        character(512) :: nf_ell !< nom du fichier de sortie des elements elliptiques
        character(512) :: nf_sch !< nom du fichier de sortie des donnees du schema adaptatif
        
        integer :: int_var   !< variables de l'integrateur (: 0 : Jacobi  1 :helio canonique)

       contains
                      
        procedure ::create_filename => extfilepar_create_filename  ! creation des noms de fichiers de sortie                      
      end type t_extfilepar  

      contains
      
!***********************************************************************
!> @brief lit le fichier de parametres name et remplit this
!***********************************************************************
      subroutine filepar_fread(this, fname)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres
       character(len=*), intent(in) :: fname     !< nom du fichier a lire
       
        character(100)  :: chemin,nf_rad,int_type
        character(100) :: nf_initext
        integer :: if_int,if_ell,if_car
        integer(8) :: n_iter,n_out
        integer :: out_ell,if_invar, ref_gmsun
        real(TREAL) :: tinit,dt
        integer :: type_pas

      namelist/lect/chemin,nf_rad,nf_initext,
     &   int_type,tinit,dt,n_iter,n_out,
     &   out_ell,int_type, ref_gmsun,
     &   if_invar,if_int,if_ell,if_car, type_pas
    
      ref_gmsun = 1
      type_pas = 0
      open(10,file=fname,status='old')
      read(10,lect)
      close(10)
      
      this%chemin = chemin
      this%nf_rad = nf_rad
      this%nf_initext = nf_initext
      this%int_type = int_type
      this%tinit = tinit
      this%dt = dt
      this%n_iter = n_iter
      this%n_out = n_out
      this%out_ell = out_ell
      this%int_type = int_type
      this%if_invar = if_invar
      this%if_int = if_int
      this%if_ell = if_ell
      this%if_car = if_car
      this%ref_gmsun = ref_gmsun
      this%type_pas = type_pas
       
      end  subroutine filepar_fread  
          
!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier  nupar
!***********************************************************************
      subroutine filepar_fwrite(this, nupar)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres
       integer, intent(in) :: nupar     !< numero de fichier

      write(nupar,*)' chemin     = ',                                    &
     & this%chemin(1:len_trim(this%chemin))
      write(nupar,*)' nf_rad     = ',                                    &
     & this%nf_rad(1:len_trim(this%nf_rad))
      write(nupar,*)' nf_initext = ',                                    &
     & this%nf_initext(1:len_trim(this%nf_initext))
      write(nupar,*)' int_type   = ',trim(this%int_type)
      write(nupar,*)' ref_gmsun  = ',this%ref_gmsun
      write(nupar,*)' tinit      = ',this%tinit
      write(nupar,*)' dt         = ',this%dt
      write(nupar,*)' n_iter     = ',this%n_iter
      write(nupar,*)' n_out      = ',this%n_out
      write(nupar,*)' out_ell    = ',this%out_ell
      write(nupar,*)' if_invar   = ',this%if_invar
      write(nupar,*)' if_int     = ',this%if_int
      write(nupar,*)' if_ell     = ',this%if_ell
      write(nupar,*)' if_car     = ',this%if_car
      write(nupar,*)' type_pas     = ',this%type_pas

      end  subroutine filepar_fwrite  

!***********************************************************************
!> @brief ecrit  des parametres de this dans le fichier 
!! dont le nom est genere a partir de chemin  et nf_rad
!***********************************************************************
      subroutine filepar_fwritepar(this)
       implicit none
       class(t_filepar), intent(inout):: this    !< fichier de parametres

      character(100) nf_par
      integer nu_par

      nf_par = trim(this%chemin)//trim(this%nf_rad)//'.par'
      nu_par = 11
      open(nu_par,file=nf_par,status='unknown')
      call filepar_fwrite(this, nu_par)
      close(nu_par)
      
      end  subroutine filepar_fwritepar        

!***********************************************************************
!> creation des noms de fichiers de sortie
!***********************************************************************
      subroutine extfilepar_create_filename(this, supprad)
      implicit none
      class(t_extfilepar), intent(inout):: this    !< fichier de parametres
      character(len=*), intent(in) :: supprad !< extension a ajouter au radical
      integer :: l_ch,l_rad,l_nf,ll
      character(512) nf
    
      l_ch = len_trim(this%chemin)
      l_rad = len_trim(this%nf_rad)
      ll = len_trim(this%int_type)
      l_nf = l_ch+l_rad
      nf = trim(this%chemin)//trim(this%nf_rad)
      nf = trim(nf)//'_'//trim(supprad)//'_'//trim(this%int_type)
      this%nf_int = trim(nf)//'.int'
      this%nf_ell = trim(nf)//'.ell'
      this%nf_car = trim(nf)//'.car'
      this%nf_sch = trim(nf)//'.sch'      
      write(*,*) trim(this%nf_int)
      write(*,*) trim(this%nf_ell)
      write(*,*) trim(this%nf_car)
      write(*,*) trim(this%nf_sch)
     
      end  subroutine extfilepar_create_filename    

      end module mod_filepar
