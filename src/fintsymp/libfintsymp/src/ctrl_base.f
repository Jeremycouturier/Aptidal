!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ctrl_base.f 
!!  \brief Controle d'une quantite 
!!
!!    --------------   ctrl_base        ---------------------------\n
!!    | buffer src |  ----------------> | buffer dst (optionnel)  |\n
!!    --------------                    ---------------------------\n
!!
!!   on envoie les donnees uniquement en cas d'erreur au buffer de sortie.
!!   le buffer de sortie est optionnel. 
!!
! history : creation 25/01/2018
!***********************************************************************
#include "realc.f"


!***********************************************************************
! module pour le controle de base
!***********************************************************************
      module mod_ctrl_base
       use mod_finalconsumer

!***********************************************************************
!> @class t_ctrlbase
!! Classe de base pour les controles provoquant un arret.
!! Elle peut avoir un buffer de sortie optionnel de 1 seul element
!! On derive les classes utiles de celle-ci. \n
!!  
!***********************************************************************
      type, extends(t_finalconsumer), abstract :: t_ctrl_base
        type(t_buffer), public, allocatable :: m_output_buffer !< consommateur du buffer en sortie
        integer, public :: m_ndata !< nombre de donnees a transmettre 
          
      contains
         
          
        procedure,public :: init => ctrl_base_init  ! initialise le controleur
        procedure,public :: set_output => ctrl_base_set_output ! fixe le buffer de sortie et la fonction a appeller lorsque le buffer est plein

        ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
        procedure(i_processinputbuffer),public, DEFERRED  ::                &
     &           processinputbuffer    ! fonction appelee pour le traitement du buffer


      end type t_ctrl_base  

           abstract interface

!***********************************************************************
! interface pour le traitement du buffer :
! retourne l'indice de la derniere erreur dans le buffer (0 si pas d'erreur)
!***********************************************************************
            function i_processinputbuffer(this, buffer, poswrite)
             import t_ctrl_base
             import t_buffer
             class(t_ctrl_base), intent(inout) :: this !< destructeur
             class(t_buffer), intent(inout) :: buffer !< buffer de sortie
             integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
             integer(8) i_processinputbuffer
            
            end function i_processinputbuffer

          end interface

      contains

!***********************************************************************
!> @brief intialise le controleur avec le nom du graphe et le buffer d'input
!***********************************************************************
         subroutine ctrl_base_init(this, ndata, dotgraphname)
          implicit none
          class(t_ctrl_base), intent(inout):: this  !< dummy argument
          integer, intent(in) :: ndata !< nombre de donnees a transmettre
          character(len=*), intent(in) :: dotgraphname !< nom du controleur 
          procedure(t_procbufferfull),pointer ::pfcnfull 
          procedure(t_procbufferfull),pointer ::pfcnflush 
          procedure(t_procgraphdot), pointer :: pfcndot
          
          this%m_ndata = ndata
          call this%set_graphdotname(dotgraphname)
          pfcnfull=>ctrl_base_onbufferfull_cb
          pfcnflush=>ctrl_base_onbufferflush_cb
          pfcndot=>ctrl_base_ongraphdot_cb
          call this%setconsumer(pfcnfull,pfcnflush, pfcndot)

         end  subroutine ctrl_base_init

!***********************************************************************
!> @brief fixe la taille du buffer de sortie et son consommateur (optionnel)
!***********************************************************************
         subroutine ctrl_base_set_output(this, buffercs)
          implicit none
          class(t_ctrl_base), intent(inout):: this  !< dummy argument
          class(t_buffer_consumer), pointer, intent(in) :: buffercs !< consommateur du buffer 
          
          allocate(this%m_output_buffer)
          call buffer_init(this%m_output_buffer, 1_8, this%m_ndata,     &
     &             buffercs) 
         end  subroutine ctrl_base_set_output

!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!***********************************************************************
          subroutine ctrl_base_onbufferfull(this,buffer,poswrite)
           implicit none
           class(t_ctrl_base), intent(inout) :: this    !< de type t_ctrl_base
           class(t_buffer), intent(inout) :: buffer !< buffer d'entree
           integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)
           integer(8) lasterror
           real(TREAL), dimension(1:this%m_ndata) :: R

           lasterror = this%processinputbuffer(buffer, poswrite)
           if ((lasterror.ne.0).and.                                     &
     &           allocated(this%m_output_buffer)) then
                 call buffer%readdata(lasterror, R)
                 call this%m_output_buffer%writedata(R)
           endif

          end subroutine ctrl_base_onbufferfull

!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree est plein
!***********************************************************************
          subroutine ctrl_base_onbufferfull_cb(this,userdata,poswrite)
           implicit none
           class(t_buffer), intent(inout) :: this !< buffer de sortie
           class(*), intent(inout) :: userdata    !< de type t_ctrl_base
           integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)

           select type(userdata)
            class is(t_ctrl_base)
             call ctrl_base_onbufferfull(userdata, this, poswrite)

            class default
             stop 'ctrl_base_onbufferfull_cb : bad class'
           end select 
           
           call this%empty()
          
          end subroutine ctrl_base_onbufferfull_cb

!***********************************************************************
!> @brief fonction appellee lorsque le buffer d'entree doit etre flushe
!***********************************************************************
          subroutine ctrl_base_onbufferflush_cb(this,userdata,poswrite)
           implicit none
           class(t_buffer), intent(inout) :: this !< buffer de sortie
           class(*), intent(inout) :: userdata    !< de type t_ctrl_base
           integer(8), intent(in) :: poswrite !< derniere position ecrite dans buffer d'ecriture (les indices a lire vont de 1 a poswrite inclus)

           call ctrl_base_onbufferfull_cb(this, userdata, poswrite)
           select type(userdata)
            class is(t_ctrl_base)
             if (allocated(userdata%m_output_buffer)) then
                 call userdata%m_output_buffer%flushdata()
             endif

            class default
             stop 'ctrl_base_onbufferflush_cb : bad class'
           end select 
           
           call this%empty()
          
          end subroutine ctrl_base_onbufferflush_cb

!***********************************************************************
!> @brief  fonction appellee pour generer le graph au format dot
!***********************************************************************
       subroutine ctrl_base_ongraphdot_cb(this,userdata, nodesrc,          &
     &                   nodecounter, num_file)
        implicit none
        class(t_buffer), intent(inout) :: this !< buffer ayant genere l'appel
        class(*), intent(inout) :: userdata !< donnee utilisateur de type t_plan_converter
        integer, intent(in) :: nodesrc !< numero du noeud  source 
        integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
        integer, intent(in) :: num_file !<numero du fichier de sortie
             
        integer curnode

        curnode = nodecounter+1

          select type(userdata)
           class is(t_ctrl_base)
             call buffer_dot_writelabel(num_file,curnode,                 &
     &          trim(userdata%m_dotname))
             call buffer_dot_writenode(num_file, nodesrc,curnode)
             nodecounter = curnode
             if (allocated(userdata%m_output_buffer)) then
                call userdata%m_output_buffer%graphdot(curnode,         &
     &          nodecounter, num_file)
             endif
           class default
            stop 'ctrl_base_ongraphdot_cb : bad class'
          end select 

       end subroutine ctrl_base_ongraphdot_cb

      end module mod_ctrl_base

