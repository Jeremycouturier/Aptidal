!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file syspart.f 
!!  \brief systeme avec M particules 
!!         pouvant etre integre par un schema symplectique
!!
!!     cette classe t_syspart fournit les routines elementaires communes a l'ensemble 
!!      des systemes possibles
!!
!!
!
! history : creation 26/06/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec M particules 
!***********************************************************************
      module mod_syspart
      use mod_buffer
      use mod_arret
      use mod_sysextension
      
!***********************************************************************
!> @class t_syspart
!! classe de base decrivant  M particules
!!
!!  
!***********************************************************************
      type, abstract, extends(t_sysextension) ::  t_syspart
          ! ces champs ne doivent etre acceder en lecture que par les fonctions exterieures au module
                   
          integer :: m_part_nb  !< nombre de particules 
          integer :: m_plan_nb !< nombre  de planetes
          real(TREAL), dimension(:,:), allocatable :: m_part_coordX !< coordonnees (position) des particules: (3, m_part_nb)
          real(TREAL), dimension(:,:), allocatable :: m_part_coordXP !< coordonnees (vitesse) des particules: (3, m_part_nb)
          real(TREAL), dimension(:,:), allocatable :: m_part_errX !< erreur sur les coordonnees (position) des particules : (3, m_part_nb)
          real(TREAL), dimension(:,:), allocatable :: m_part_errXP !< erreur sur les  coordonnees (vitesse) des particules : (3, m_part_nb)
          real(TREAL), dimension(:), allocatable :: m_part_dataout !< tableau temporaire pour copier les particules+le temps dans le buffer
         
          real(TREAL) :: m_star_mass !< masse de l'etoile
          real(TREAL) :: m_cG !< constante de Gauss AU**3/an**2

          
          type(t_arret), dimension(:), allocatable :: m_part_mce  !< tableau de controles d'erreurs de chaque particule : (1:m_part_nb)

       contains
          
          
          procedure:: set_nb => syspart_set_nb ! fixe le nombre de particules 
          procedure:: get_nb => syspart_get_nb ! retourne le nombre de particules 
          procedure:: set_coord => syspart_set_coord ! fixe les coordonnees des particules (meme que celle integree)
          
          procedure:: set_star_mass => syspart_set_star_mass ! fixe la masse de l'etoile
          
          procedure:: create_output    => syspart_create_output ! creation de this dans la copie de sortie
          procedure:: copy_output    => syspart_copy_output ! copie de this dans la copie de sortie
          procedure :: write_output    => syspart_write_output ! ecriture de la sortie (appellee par l'integrateur)
          procedure:: copy_output_feedback=>syspart_copy_output_feedback ! copie de la copie de sortie vers this des informations necessaires 
         
          procedure :: add_part_out_step => syspart_add_part_out_step ! fixe le pas de sortie, la taille du buffer des planetes

          ! controle des erreurs
          procedure :: set_error_time => syspart_set_error_time ! fixe le temps pour le controle des erreurs
!          procedure :: check_error => syspart_check_error ! retourne true si une erreur s'est produite
!          procedure :: get_error  => syspart_get_error ! recupere l'erreur
!          procedure :: set_error  => syspart_set_error ! fixe l'erreur
         
      end type t_syspart  
     
      contains
         

!***********************************************************************
!> @brief fixe le nombre de particules et le nombre de planetes
!***********************************************************************
      subroutine syspart_set_nb(this, nplan, npart)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       integer, intent(in) :: nplan                !< nombre de planetes
       integer, intent(in) :: npart                !< nombre de particules
       
       call this%set_graphdotname("PARTICULES")
       this%m_part_nb = npart
       this%m_plan_nb = nplan
       
      end  subroutine syspart_set_nb  
      
!***********************************************************************
!> @brief retourne le nombre de particules 
!***********************************************************************
      function syspart_get_nb(this)
       implicit none
       class(t_syspart), intent(in):: this  !< dummy argument
       integer :: syspart_get_nb       
       
       syspart_get_nb = this%m_part_nb
      end  function syspart_get_nb  

!***********************************************************************
!> @brief fixe les coordonnees des particules. 
!! Il n'y a aucune conversion des coordonnees.
!! Ces coordoonnees doivent etre compatibles avec l'integrateur
!***********************************************************************
      subroutine syspart_set_coord(this, X, XP)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:,:), intent(in) :: X  !< position des particules
       real(TREAL), dimension(:,:), intent(in) :: XP  !< moment des particules
       integer npart
       
       npart = this%m_part_nb
       allocate(this%m_part_coordX(3, npart))
       if( any( shape(this%m_part_coordX) /= shape(X) ) ) then
          stop "Erreur : set_part_coord : X : dimensions differentes"
       endif
       this%m_part_coordX = X
       allocate(this%m_part_coordXP(3, npart))
       if( any( shape(this%m_part_coordXP) /= shape(XP) ) ) then
          stop "Erreur : set_part_coord : XP : dimensions differentes"
       endif
       this%m_part_coordXP = XP
       
       ! tableau pour la sommation compensee
       allocate(this%m_part_errX(3, npart))
       this%m_part_errX = REALC(0.E0)
       allocate(this%m_part_errXP(3, npart))
       this%m_part_errXP = REALC(0.E0)
       
       ! tableau pour les arrets
       allocate(this%m_part_mce(1:npart))

      end  subroutine syspart_set_coord  


!***********************************************************************
!> @brief fixe la masse de l'etoile 
!***********************************************************************
      subroutine syspart_set_star_mass(this, cG, m0)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: m0          !< masse de l'etoile
       real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
       
       this%m_star_mass = m0
       this%m_cG = cG
      end  subroutine syspart_set_star_mass  
 
!***********************************************************************
!> @brief met a jour le temps d'une erreur.
!! indique le comportement du programme qui l'appelle par la suite
!***********************************************************************
      subroutine syspart_set_error_time(this, t)
       implicit none
       class(t_syspart), intent(inout):: this  !< systeme de particules
       real(TREAL), intent(in):: t    !< temps a laquelle s'est produite l'erreur
       
      ! met a jour si les temps d'arret ne sont pas initialises
      call array_arret_update_time(this%m_part_mce, t)
                    
      end subroutine syspart_set_error_time  

#if 0     

!***********************************************************************
!> @brief controle des erreurs du systeme planetaire
!! @return retourne false  si aucune erreur \n
!! retourne true si une erreur s'est produite depuis le dernier appel a cette fonction
!***********************************************************************
      function syspla_check_error(this) result(ret)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       
       logical ret ! code de retour
       ret = .false.
       
         if (this%m_plan_arret%m_stop.ne.0) then
                ret = .true.
         endif
         
      end function  syspla_check_error  

!***********************************************************************
!> @brief recupere l'erreur du systeme
!***********************************************************************
      function syspart_get_error(this) result(arret)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       
       type(t_arret) :: arret ! code de retour
       arret = this%m_plan_arret
         
      end function  syspart_get_error  

!***********************************************************************
!> @brief fixe l'erreur du systeme
!***********************************************************************
      subroutine syspart_set_error(this, arret)
       implicit none
       class(t_syspla), intent(inout):: this  !< dummy argument
       type(t_arret), intent(in) :: arret !< erreur
       
       this%m_plan_arret = arret
       
      end subroutine  syspart_set_error  
#endif

!***********************************************************************
!> @brief creation d'une copie pour effectuer les pas de sortie
!! en general, duplique simplement this
!***********************************************************************
      subroutine syspart_create_output(this, sysout)
       implicit none
       class(t_syspart), intent(in):: this  !< dummy argument
       class(t_sysextension), intent(inout) :: sysout  !< systeme cree pour la sortie
       integer npart
       
       npart = this%m_part_nb
       
       select type(sysout)
        class is(t_syspart)
         if (allocated(sysout%m_part_coordX)) then
           deallocate(sysout%m_part_coordX);
           deallocate(sysout%m_part_coordXP);
           deallocate(sysout%m_part_errX);
           deallocate(sysout%m_part_errXP); 
         endif
         allocate(sysout%m_part_coordX(3, npart),                        &
     &    SOURCE=this%m_part_coordX)
         allocate(sysout%m_part_coordXP(3, npart),                       &
     &    SOURCE=this%m_part_coordXP)
         allocate(sysout%m_part_errX(3, npart),                          &
     &    SOURCE=this%m_part_errX)
         allocate(sysout%m_part_errXP(3, npart),                         &
     &    SOURCE=this%m_part_errXP)
        class default
         stop 'Erreur: this ne derive pas de syspart_create_output' 
      end select
      end  subroutine syspart_create_output
        
!***********************************************************************
!> @brief copie de this dans la copie de sortie (sysout)
!! en general, sysout = this 
!***********************************************************************
      subroutine syspart_copy_output(this, sysout)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       class(t_sysextension), intent(inout) :: sysout  !< systeme cree pour la sortie

       select type(sysout)
        class is(t_syspart)
         sysout%m_part_coordX(:,:) = this%m_part_coordX(:,:)
         sysout%m_part_coordXP(:,:) = this%m_part_coordXP(:,:)
         sysout%m_part_errX(:,:) = this%m_part_errX(:,:)
         sysout%m_part_errXP(:,:) = this%m_part_errXP(:,:)
        class default
         stop 'Erreur: this ne derive pas de syspart_copy_output' 
      end select
      
      end  subroutine syspart_copy_output
 

!***********************************************************************
!> @brief ecriture de la sortie (appellee par l'integrateur)
!!   a l'iteration donne et au temps donne
!***********************************************************************
      subroutine syspart_write_output(this, syspla, it, t)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       class(t_syspla), intent(inout):: syspla  !< systeme planetaire associee a l'extension
       integer(8), intent(in) :: it !< iteration actuelle
       real(TREAL), intent(in) :: t  !< temps correspondant a it
       integer :: j
       integer npart, nplan
       integer ideb, ifin

       npart = this%get_nb()
       nplan = syspla%get_plan_nb()

       do j=1, size(this%m_stepout)

        if ((mod(it,this%m_stepout(j)).eq.0).and.                            &
     &     (syspla%m_localstepkind.eqv..true.)) then
!        if (syspla%isvalidoutput(it,j)) then
        
         ! cas d'un buffer de particules
         if (this%m_kind(j).eq.1) then
          ! write(*,*) " sortie part ",it,t
          this%m_part_dataout(1) = t
          ! pour les planetes
          ideb = 2
          ifin = 1+3*nplan
          this%m_part_dataout(ideb:ifin)=reshape(                           &
     &                 syspla%m_plan_coordX,                                &
     &                 shape(this%m_part_dataout(1:3*nplan)))
          ideb = 2+3*nplan
          ifin = 1+6*nplan
          this%m_part_dataout(ideb:ifin)=reshape(                           &
     &                 syspla%m_plan_coordXP,                               &
     &                 shape(this%m_part_dataout(1:3*nplan)))
          ! pour les particules
          ideb = 2+6*nplan
          ifin = 1+6*nplan+3*npart
           this%m_part_dataout(ideb:ifin)=reshape(this%m_part_coordX,      &
     &                 shape(this%m_part_dataout(1:3*npart)))
          ideb = 2+6*nplan+3*npart
          ifin = 1+6*nplan+6*npart
          this%m_part_dataout(ideb:ifin)=reshape(this%m_part_coordXP,      &
     &                 shape(this%m_part_dataout(1:3*npart)))

          call this%set_error_time(t)
          call this%m_buffer(j)%set_multictrlerr(this%m_part_mce)
          
          call this%m_buffer(j)%writedata(this%m_part_dataout)
          
          call this%m_buffer(j)%get_multictrlerr(this%m_part_mce)
          if (this%m_buffer(j)%update_count_multictrlerr()) then
            !call syspla%m_plan_arret%set_error(-8)
            call syspla%m_plan_arret%set_time(t)
            call syspla%m_plan_arret%set_value(real(0.,kind=TREAL))
            call syspla%m_plan_arret%set_body(-1)
          endif

         ! verifie si une erreur s'est produite dans le successeur
         if (this%m_buffer(j)%check_successor_error()) then
         call this%m_buffer(j)%get_successor_error(syspla%m_plan_arret)
         endif

         endif
        endif
       enddo
       
      end  subroutine syspart_write_output

!***********************************************************************
!> @brief copie de this dans dst des informations transmises 
!! par les buffers de sortie
!! en general, les erreurs
!***********************************************************************
      subroutine syspart_copy_output_feedback(this, dst)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       class(t_sysextension), intent(inout) :: dst  !< systeme distant

       select type(dst)
        class is(t_syspart)
          dst%m_part_mce(:) = this%m_part_mce(:)
        class default
         stop 'Erreur: this ne derive pas de syspart_copy_output_fe.' 
      end select

      end  subroutine syspart_copy_output_feedback


!***********************************************************************
!> @brief fixe le pas de sortie, la taille du buffer des particules
!! cela cree un nouveau buffer
!! lorsque le buffer est plein, la fonction usercallback est appellee avec userdata
!***********************************************************************
      subroutine syspart_add_part_out_step(this, stepout, buffersize,           &
     &                             buffercs)
       implicit none
       class(t_syspart), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les particules
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 

       integer ntot, ind
       
       ntot = this%get_nb()+this%m_plan_nb
       call this%add_out_step(stepout, buffersize,                              &
     &                               buffercs, 1+ntot*6, 1)
       if (.not.allocated(this%m_part_dataout)) then
           allocate(this%m_part_dataout(1:1+ntot*6))
       endif    
       
       ind = size(this%m_stepout)
       call this%m_buffer(ind)%init_multictrlerr(this%get_nb())
       
      end  subroutine syspart_add_part_out_step  

      end module mod_syspart
