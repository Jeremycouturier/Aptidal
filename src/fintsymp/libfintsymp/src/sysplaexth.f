!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysplaexth.f 
!!  \brief systeme avec N planetes
!!         pouvant etre integre par un schema symplectique.
!!         les coordonnees des corps sont exprimees en heliocentriques (position/vitesse).
!!         les interactions newtonniennes sont prises en compte ainsi que des effets complementaires.
!!
!!         ces effets ou objets complementaires doivent derivees de t_sysextension.
!!         ils sont ajoutes par set_extension_nb et add_extension(j,...)
!!
!
! history : creation 20/08/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec N planetes 
! exprime en variables heliocentriques
! les interactions newtonniennes sont prises en compte ainsi que des effets complementaires
!***********************************************************************
      module mod_sysplaextH
       use mod_sysplanewtH
       use mod_sysextension

!***********************************************************************
!> @class t_sysplaextH
!! classe decrivant 1 systeme planetaire avec N planetes
!! exprimee en variables heliocentriques.
!! les interactions newtonniennes sont prises en compte ainsi que des effets complementaires
!!  
!!  
!***********************************************************************
      type, extends(t_sysplanewtH) :: t_sysplaextH 

          integer ::m_nbextension = 0 !< nombre d'extensions 
          type(t_sysextension_ptr), dimension(:), allocatable :: m_ext !< tableau des extensions
          
       contains
          
          procedure, public:: set_extension_nb =>                        &
     &       sysplaextH_set_extension_nb ! fixe le nombre d'extensions
          
          procedure, public:: add_extension => sysplaextH_add_extension ! ajoute l'extension
          procedure, public:: get_extension => sysplaextH_get_extension ! retourne l'extension d'indice fixe

          procedure :: pasA =>sysplaextH_pasA
          procedure :: pasB =>sysplaextH_pasB
          procedure :: pasC =>sysplaextH_pasC
          
          procedure:: create_output  => sysplaextH_create_output! creation d'une copie pour effectuer les pas de sortie
          procedure:: copy_output    => sysplaextH_copy_output ! copie de this dans la copie de sortie
          procedure, NOPASS:: delete_output =>sysplaextH_delete_output ! destruction de l'objet cree par create_output
          procedure :: write_output    => sysplaextH_write_output ! ecriture de la sortie (appellee par l'integrateur)
          procedure:: copy_output_feedback=>                             &
     &        sysplaextH_copy_output_feedback ! copie de la copie de sortie vers this des informations necessaires 
          procedure :: flush_output    => sysplaextH_flush_output ! forcer a transmettre les buffers (appellee par l'integrateur)
          
          procedure :: get_outputsteps => sysplaextH_get_outputsteps ! fourni les pas de sorties  (appellee par l'integrateur)

          procedure :: set_error_time => sysplaextH_set_error_time ! fixe le temps pour le controle des erreurs

          procedure, public :: graphdot => sysplaextH_graphdot ! affiche le graphique dot

          final :: sysplaextH_destructor ! destructeur
          
      end type t_sysplaextH  
      
      contains
         
!***********************************************************************
!> @brief destructeur
!***********************************************************************
      subroutine sysplaextH_destructor(this)
       implicit none
       type(t_sysplaextH), intent(inout):: this  !< dummy argument
       integer j
       
       do j=1, this%m_nbextension
          if (associated(this%m_ext(j)%m_ptr)) then
            deallocate (this%m_ext(j)%m_ptr)
          endif  
       enddo
      end  subroutine sysplaextH_destructor  



!***********************************************************************
!> @brief fixe le nombre d'extensions (effets ou objets)
!***********************************************************************
      subroutine sysplaextH_set_extension_nb(this, nb)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       integer, intent(in) :: nb  !< nombre d'extensions a prevoir
       
       this%m_nbextension = nb
       allocate(this%m_ext(nb))
       call this%set_graphdotname("SYSEXT Helio")
      end  subroutine sysplaextH_set_extension_nb  

!***********************************************************************
!> @brief ajoute l'extension a l'indice j (j commence a 1)
!! l'extension est duplique
!***********************************************************************
      subroutine sysplaextH_add_extension(this, j, ext)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       integer, intent(in) :: j  !< indice de l'extension
       class(t_sysextension), intent(in):: ext  !< extension a inserer
       
       allocate(this%m_ext(j)%m_ptr, source=ext)
      end  subroutine sysplaextH_add_extension  

!***********************************************************************
!> @brief retourne un pointeur vers l'extension a l'indice j (j commence a 1)
!***********************************************************************
      subroutine sysplaextH_get_extension(this, j, ext)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       integer, intent(in) :: j  !< indice de l'extension
       class(t_sysextension), pointer, intent(inout):: ext  !< extension a inserer
       
       ext => this%m_ext(j)%m_ptr
      end  subroutine sysplaextH_get_extension  


!***********************************************************************
!> @brief execute le pas A de l'integrateur
!***********************************************************************
      subroutine sysplaextH_pasA(this, mdt)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= tau
       
       integer j
       
        call sysplanewtH_pasA(this,mdt)
        do j=1, this%m_nbextension
          call this%m_ext(j)%m_ptr%pasA(this, mdt)
        enddo
     
      end  subroutine sysplaextH_pasA  


!***********************************************************************
!> @brief execute le pas B de l'integrateur = B1(tau)B2(tau)B1(tau) 
!***********************************************************************
      subroutine sysplaextH_pasB(this, mdt)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
        integer j

       call sysplanewtH_pasB1(this,mdt/REALC(2.E0))
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasB1(this, mdt/REALC(2.E0))
       enddo
       call sysplanewtH_pasB2(this,mdt)
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasB2(this, mdt)
       enddo
       call sysplanewtH_pasB1(this,mdt/REALC(2.E0))
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasB1(this, mdt/REALC(2.E0))
       enddo

      end  subroutine sysplaextH_pasB  

!***********************************************************************
!> @brief execute le pas C de l'integrateur
!***********************************************************************
      subroutine sysplaextH_pasC(this, mdt)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt

       integer j
       
       call sysplanewtH_pasC(this,mdt)
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasC(this, mdt)
       enddo
       
      end  subroutine sysplaextH_pasC  

!***********************************************************************
!> @brief creation d'une copie pour effectuer les pas de sortie
!! en general, duplique simplement this
!***********************************************************************
      subroutine sysplaextH_create_output(this, sysout)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       class(t_syssymp), allocatable, intent(out) :: sysout  !< systeme cree pour la sortie
       
       integer j
       call syspla_create_output(this, sysout)
       select type(sysout)
        class is(t_sysplaextH)
         if (allocated(sysout%m_ext)) then
          deallocate(sysout%m_ext)
         endif
         call sysout%set_extension_nb(this%m_nbextension)
         do j=1, this%m_nbextension
          allocate(sysout%m_ext(j)%m_ptr, source=this%m_ext(j)%m_ptr)
          call this%m_ext(j)%m_ptr%create_output(sysout%m_ext(j)%m_ptr)
         enddo
        class default
         stop 'Erreur: this ne derive pas de sysplaextH_create_output' 
      end select

      end  subroutine sysplaextH_create_output
        
!***********************************************************************
!> @brief destruction de l'objet cree par create_output
!! en general, detruit  this
!***********************************************************************
      subroutine sysplaextH_delete_output(this)
       implicit none
       class(t_syssymp), allocatable, intent(inout):: this  !< dummy argument
       
       deallocate(this)

      end  subroutine sysplaextH_delete_output
      
!***********************************************************************
!> @brief copie de this dans la copie de sortie (sysout)
!! en general, sysout = this 
!***********************************************************************
      subroutine sysplaextH_copy_output(this, sysout)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       class(t_syssymp), intent(inout) :: sysout  !< systeme cree pour la sortie

       integer j
       call syspla_copy_output(this, sysout)
       select type(sysout)
        class is(t_sysplaextH)
         do j=1, this%m_nbextension
          call this%m_ext(j)%m_ptr%copy_output(sysout%m_ext(j)%m_ptr)
         enddo
              
        class default
         stop 'Erreur: this ne derive pas de sysplaextH_copy_output' 
      end select
      end  subroutine sysplaextH_copy_output
 
!***********************************************************************
!> @brief ecriture de la sortie (appellee par l'integrateur)
!!   a l'iteration donne et au temps donne
!***********************************************************************
      subroutine sysplaextH_write_output(this, it, t)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       integer(8), intent(in) :: it !< iteration actuelle
       real(TREAL), intent(in) :: t  !< temps correspondant a it
       integer :: j

       call syspla_write_output(this,it,t)
       do j=1, this%m_nbextension
        call this%m_ext(j)%m_ptr%write_output(this, it,t)
       enddo
       
      end  subroutine sysplaextH_write_output

!***********************************************************************
!> @brief force a transmettre le contenu des buffers (appellee par l'integrateur)
!***********************************************************************
      subroutine sysplaextH_flush_output(this)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       integer :: j

       call syspla_flush_output(this)
       do j=1, this%m_nbextension
            call this%m_ext(j)%m_ptr%flush_output()
       enddo
       
      end  subroutine sysplaextH_flush_output

!***********************************************************************
!> @brief retourne les pas de sortie ou l'integrateur doit appeller la fonction 
!!        write_output
!***********************************************************************
      subroutine sysplaextH_get_outputsteps(this,nstepscall,arstepscall)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       integer, intent(out) :: nstepscall  !< taille de arstepscall
       integer(8), dimension(:), allocatable, intent(out) :: arstepscall !< tableau des pas de sortie
       integer j, index

       !comptage
       nstepscall = 0
       if (allocated(this%m_stepout)) then
            nstepscall = size(this%m_stepout)
       endif
       
       do j=1, this%m_nbextension
          nstepscall=nstepscall+this%m_ext(j)%m_ptr%get_noutputsteps()
       enddo
       
       ! allocation et remplissage
       allocate(arstepscall(1:nstepscall))
       index = 0
       if (allocated(this%m_stepout)) then
            index = size(this%m_stepout)
            arstepscall(1:index) = this%m_stepout
       endif
       index=index+1
       do j=1, this%m_nbextension
          call this%m_ext(j)%m_ptr%get_outputsteps(index, arstepscall)
       enddo
       
      end  subroutine sysplaextH_get_outputsteps

!***********************************************************************
!> @brief copie de this dans dst des informations transmises 
!! par les buffers de sortie
!! en general, les erreurs
!***********************************************************************
      subroutine sysplaextH_copy_output_feedback(this, dst)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< dummy argument
       class(t_syssymp), intent(inout) :: dst  !< systeme distant
       integer j

       select type(dst)
        class is(t_sysplaextH)
        do j=1, this%m_nbextension
          call this%m_ext(j)%m_ptr%copy_output_feedback(                   &
     &              dst%m_ext(j)%m_ptr)
        enddo
              
        class default
         stop 'Erreur: this ne derive pas de sysplaextH_copy_outp....' 
      end select
      end  subroutine sysplaextH_copy_output_feedback

!***********************************************************************
!> @brief met a jour le temps d'une erreur.
!! indique le comportement du programme qui l'appelle par la suite
!! @return retourne 0 si le programme doit se poursuivre \n
!! retourne 1 si on doit s'arreter d'integrer (mais poursuite de l'execution)
!***********************************************************************
      function sysplaextH_set_error_time(this, t)  result(ret)
       implicit none
       class(t_sysplaextH), intent(inout):: this  !< systeme planetaire
       real(TREAL), intent(in):: t    !< temps a laquelle s'est produite l'erreur
       integer ret ! code de retour
       integer j
       
       ret = syspla_set_error_time(this, t)
        do j=1, this%m_nbextension
          call this%m_ext(j)%m_ptr%set_error_time(t)
        enddo
                    
      end function sysplaextH_set_error_time  

!***********************************************************************
!> @brief  fonction pour generer le graph au format dot
!***********************************************************************
       subroutine sysplaextH_graphdot(this,nodesrc,nodecounter,num_file)
        implicit none
        class(t_sysplaextH), intent(inout) :: this  !< dummy argument
        integer, intent(in) :: nodesrc !< numero du noeud  source 
        integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
        integer, intent(in) :: num_file !<numero du fichier de sortie
        
        integer curnode, j

        curnode = nodecounter+1
        call buffer_dot_writelabel(num_file, curnode,                         &
     &         trim(this%m_dotname))
        call buffer_dot_writenode(num_file, nodesrc,curnode)
        nodecounter = curnode

        do j=1, size(this%m_buffer)
           call this%m_buffer(j)%graphdot(curnode,nodecounter,num_file)
        enddo
        do j=1, this%m_nbextension
        call this%m_ext(j)%m_ptr%graphdot(curnode,nodecounter,num_file)
        enddo
        
       end subroutine sysplaextH_graphdot

      end module mod_sysplaextH
      