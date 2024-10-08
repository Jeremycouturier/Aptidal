!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysextension.f 
!!  \brief classe abstraite pour l'ajout d'extensions aux systemes planetaires de base
!!
!
! history : creation 20/08/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ajout d'extensions aux systemes planetaires de base
!***********************************************************************
      module mod_sysextension
         use mod_syspla
         
!***********************************************************************
!> @class t_sysextension
!! classe abstraite dont doivent deriver les extensions a ajouter 
!!  
!***********************************************************************
      type, abstract :: t_sysextension 
          
          character(len=30) :: m_dotname = "EXTENSION" !< nom de l'extension pour le graphique dot

          ! buffers de sortie
          integer(8), dimension(:), allocatable     :: m_stepout !< pas de sortie (particule)
          type(t_buffer), dimension(:), allocatable :: m_buffer !< buffer de sortie (particule)
          !> buffer de sortie
          !! = 1 : buffer pour les particules
          integer(8), dimension(:), allocatable     :: m_kind 

       contains
          
          procedure(t_extpasA), deferred :: pasA
          procedure(t_extpasA), deferred :: pasB1
          procedure(t_extpasA), deferred :: pasB2
          procedure(t_extpasA), deferred :: pasC
          
          procedure :: flush_output    => sysext_flush_output ! forcer a transmettre les buffers (appellee par l'integrateur)
          procedure :: add_out_step => sysext_add_out_step ! fixe le pas de sortie, la taille du buffer des particules
          procedure :: clear_out_step => sysext_clear_out_step ! supprime tous les buffers lies en sortie au systeme
          procedure :: get_noutputsteps => sysext_get_noutputsteps ! fourni le nombre de sortie
          procedure :: get_outputsteps => sysext_get_outputsteps ! fourni les parties de sortie

          procedure:: create_output    => sysext_create_output ! creation de this dans la copie de sortie
          procedure:: copy_output    => sysext_copy_output ! copie de this dans la copie de sortie
          procedure :: write_output    => sysext_write_output ! ecriture de la sortie (appellee par l'integrateur)
          procedure:: copy_output_feedback=>sysext_copy_output_feedback ! copie de la copie de sortie vers this des informations necessaires 

          procedure :: set_error_time => sysext_set_error_time ! fixe le temps pour le controle des erreurs

          procedure :: graphdot => sysext_graphdot ! affiche le graphique dot
          procedure :: set_graphdotname => sysext_set_graphdotname ! fixe le nom du graphique dot

      end type t_sysextension  

      type :: t_sysextension_ptr 

          class(t_sysextension),pointer :: m_ptr => NULL() !< pointeur vers une seule extension
          
      end type t_sysextension_ptr 


           abstract interface
!***********************************************************************
! type de fonction pour realiser un pas d'integration A, B, ou C
!***********************************************************************
            subroutine t_extpasA(this, syspla, tau)
             import t_sysextension
             import t_syspla
             class(t_sysextension), intent(inout) :: this
             class(t_syspla), intent(inout) :: syspla
             real(TREAL), intent(in) :: tau 
            end subroutine t_extpasA
      
          end interface

      contains
         
!***********************************************************************
!> @brief force a transmettre le contenu des buffers (appellee par l'integrateur)
!***********************************************************************
      subroutine sysext_flush_output(this)
       implicit none
       class(t_sysextension), intent(inout):: this  !< dummy argument
       integer :: j

       do j=1, size(this%m_stepout)
          call this%m_buffer(j)%flushdata()
       enddo
       
      end  subroutine sysext_flush_output

!***********************************************************************
!> @brief ajoute un nouveau buffer au tableau de ceux existants
!***********************************************************************
      subroutine sysext_add_out_step(this, stepout, buffersize,           &
     &                             buffercs, ncomp, pkind)
       implicit none
       class(t_sysextension), intent(inout):: this       !< dummy argument
       integer(8), intent(in) :: stepout           !< pas de sortie pour les planetes
       integer(8), intent(in) :: buffersize        !< nombre de pas de sortie avant d'appeller la fonction usercallback
       class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
       integer, intent (in) :: pkind !< =0 => buffer de planete, = 1 =buffer de particule         
       integer, intent (in) ::  ncomp !< nombre de composantes par pas de sorties
       integer newind
       integer(8), dimension(:), allocatable       :: newstepout 
       type(t_buffer), dimension(:), allocatable :: newbuffer
       integer(8), dimension(:), allocatable       :: newkind
       
       newind = 1
       
       if (allocated(this%m_stepout)) then
           newind = size(this%m_stepout)+1
           allocate(newstepout(1:newind))
           allocate(newbuffer(1:newind))
           allocate(newkind(1:newind))
           newstepout(1:newind-1) = this%m_stepout
           newbuffer(1:newind-1) = this%m_buffer
           newkind(1:newind-1) = this%m_kind
           call move_alloc(TO=this%m_stepout, FROM=newstepout)
           call move_alloc(TO=this%m_buffer, FROM=newbuffer)
           call move_alloc(TO=this%m_kind, FROM=newkind)
       else
           allocate(this%m_stepout(1:1))
           allocate(this%m_buffer(1:1))
           allocate(this%m_kind(1:1))
       endif
      
       this%m_stepout(newind) = stepout
       call this%m_buffer(newind)%init(buffersize, ncomp, buffercs) 
       this%m_kind(newind) = pkind
       
      end  subroutine sysext_add_out_step  


!***********************************************************************
!> @brief routine pour supprimer tous les buffers lies au systeme en sortie
!***********************************************************************
      subroutine sysext_clear_out_step(this)
       implicit none
       class(t_sysextension), intent(inout):: this       !< dummy argument
       
       if (allocated(this%m_stepout)) then
        deallocate(this%m_stepout)
        deallocate(this%m_buffer)
        deallocate(this%m_kind)
       endif
       
      end  subroutine sysext_clear_out_step  

!***********************************************************************
!> @brief creation d'une copie pour effectuer les pas de sortie
!! en general, duplique simplement this
!***********************************************************************
      subroutine sysext_create_output(this, sysout)
       implicit none
       class(t_sysextension), intent(in):: this  !< dummy argument
       class(t_sysextension),  intent(inout) :: sysout  !< systeme cree pour la sortie
       
      end  subroutine sysext_create_output
        
!***********************************************************************
!> @brief copie de this dans la copie de sortie (sysout)
!! en general, sysout = this 
!***********************************************************************
      subroutine sysext_copy_output(this, sysout)
       implicit none
       class(t_sysextension), intent(inout):: this  !< dummy argument
       class(t_sysextension), intent(inout) :: sysout  !< systeme cree pour la sortie

      end  subroutine sysext_copy_output


!***********************************************************************
!> @brief ecriture de la sortie (appellee par l'integrateur)
!!   a l'iteration donne et au temps donne
!***********************************************************************
      subroutine sysext_write_output(this, syspla, it, t)
       implicit none
       class(t_sysextension), intent(inout):: this  !< dummy argument
       class(t_syspla), intent(inout):: syspla  !< systeme planetaire associee a l'extension
       integer(8), intent(in) :: it !< iteration actuelle
       real(TREAL), intent(in) :: t  !< temps correspondant a it
       
      end  subroutine sysext_write_output

!***********************************************************************
!> @brief retourne les pas de sortie ou l'integrateur doit appeller la fonction 
!!        write_output
!***********************************************************************
      subroutine sysext_get_outputsteps(this, index, arstepscall)
       implicit none
       class(t_sysextension), intent(inout):: this  !< dummy argument
       integer, intent(inout) :: index  !< indice ou ecrire les pas de sortie (en sortie doit etre incremntee)
       integer(8), dimension(:), intent(inout) :: arstepscall !< tableau des pas de sortie

       integer nstepscall
       nstepscall = size(this%m_stepout)
       if (nstepscall.gt.0) then
        arstepscall(index:index+nstepscall-1) = this%m_stepout
        index=index+nstepscall
       endif
       
      end  subroutine sysext_get_outputsteps

!***********************************************************************
!> @brief copie de this dans dst des informations transmises 
!! par les buffers de sortie
!! en general, les erreurs
!***********************************************************************
      subroutine sysext_copy_output_feedback(this, dst)
       implicit none
       class(t_sysextension), intent(inout):: this  !< dummy argument
       class(t_sysextension), intent(inout) :: dst  !< systeme distant

      end  subroutine sysext_copy_output_feedback

!***********************************************************************
!> @brief retourne le nombre de pas de sortie ou l'integrateur doit appeller la fonction 
!!        write_output
!***********************************************************************
      function sysext_get_noutputsteps(this)
       implicit none
       class(t_sysextension), intent(inout):: this  !< dummy argument
       integer :: sysext_get_noutputsteps  !< taille de arstepscall

       sysext_get_noutputsteps = size(this%m_stepout)

      end  function sysext_get_noutputsteps

!***********************************************************************
!> @brief met a jour le temps d'une erreur.
!! indique le comportement du programme qui l'appelle par la suite
!***********************************************************************
      subroutine sysext_set_error_time(this, t)
       implicit none
       class(t_sysextension), intent(inout):: this  !< systeme de particules
       real(TREAL), intent(in):: t    !< temps a laquelle s'est produite l'erreur
       
                    
      end subroutine sysext_set_error_time  

!***********************************************************************
!> @brief  fonction pour generer le graph au format dot
!***********************************************************************
       subroutine sysext_graphdot(this, nodesrc,nodecounter,num_file)
        implicit none
        class(t_sysextension), intent(inout) :: this  !< dummy argument
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
        
       end subroutine sysext_graphdot

!***********************************************************************
!> @brief specifie le nom pour le graphique DOT utilise par this
!***********************************************************************
       subroutine sysext_set_graphdotname (this, dotname)
        implicit none
        class(t_sysextension), intent(inout) :: this  !< dummy argument
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        this%m_dotname = dotname
       end subroutine sysext_set_graphdotname

      end module mod_sysextension
      