!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file finalconsumer.f 
!!  \brief Classe de base pour les consommateurs finaux de buffer.
!! Ces classes n'envoient pas leu donnee vers d'autres buffers.
!! On derive les classes utiles de celle-ci. 
!!
!!
! history : creation 09/04/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour les consommateurs finaux de buffer
!***********************************************************************
      module mod_finalconsumer
       use mod_buffer
       
!***********************************************************************
!> @class t_finalconsumer
!! Classe de base pour les consommateurs finaux de buffer.
!! Ces classes n'envoient pas leurs donnees vers d'autres buffers.
!! On derive les classes utiles de celle-ci. \n
!!  
!***********************************************************************
      type :: t_finalconsumer
          type(t_buffer_consumer), public, pointer :: m_input_buffer     &
     &          =>NULL()!< consommateur du buffer en entree
          
          character(len=40) :: m_dotname !< nom du convertisseur pour le graphique dot

      contains
          
          generic   :: setconsumer  => finalconsumer_setconsumer,        &
     &          finalconsumer_setconsumerdot   ! fixe les fonctions a appeller pour le buffer entrant
          procedure :: finalconsumer_setconsumer   
          procedure :: finalconsumer_setconsumerdot  
          procedure :: getconsumer  => finalconsumer_getconsumer   ! retourne le consommateur du buffer

          procedure, NON_OVERRIDABLE :: set_graphdotname =>              &
     &              finalconsumer_set_graphdotname ! affiche le graphique dot

          final ::  finalconsumer_destructor ! destructor       

      end type t_finalconsumer  

     
      contains
         
!***********************************************************************
!> @brief destructor
!***********************************************************************
       subroutine finalconsumer_destructor (this)
        implicit none
        type ( t_finalconsumer ) :: this  !< dummy argument
             if (associated(this%m_input_buffer)) then
              deallocate(this%m_input_buffer)
             endif
       end subroutine finalconsumer_destructor

!***********************************************************************
!> @brief specifie le nom pour le graphique DOT utilise par this
!***********************************************************************
       subroutine finalconsumer_set_graphdotname (this, dotname)
        implicit none
        class(t_finalconsumer), intent(inout) :: this  !< dummy argument
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        this%m_dotname = dotname
       end subroutine finalconsumer_set_graphdotname

!***********************************************************************
!> @brief fixe les fonctions callback pour ce consommateur
!***********************************************************************
      subroutine finalconsumer_setconsumerdot(this, callbackonfull,          &
     &     callbackonflush, pfcndot)
       implicit none
       class(t_finalconsumer), intent(inout):: this  !< consommateur du buffer : dummy argument
       procedure(t_procbufferfull),pointer,intent(in)::callbackonfull !< fonction a appeler lorsque le buffer est plein
       procedure(t_procbufferfull),pointer,intent(in)::callbackonflush !< fonction a appeler lorsque le contenu du buffer doit etre obilgatoirement propage
       procedure(t_procgraphdot), pointer,intent(in) :: pfcndot !< fonction a appeler pour generer le graphique
       
        allocate(this%m_input_buffer)
        call this%m_input_buffer%set_owner(this, callbackonfull,             &
     &         callbackonflush,pfcndot)
       
      end subroutine  finalconsumer_setconsumerdot

!***********************************************************************
!> @brief fixe les fonctions callback pour ce consommateur
!***********************************************************************
      subroutine finalconsumer_setconsumer(this, callbackonfull,           &
     &     callbackonflush)
       implicit none
       class(t_finalconsumer), intent(inout):: this  !< consommateur du buffer : dummy argument
       procedure(t_procbufferfull),pointer,intent(in)::callbackonfull !< fonction a appeler lorsque le buffer est plein
       procedure(t_procbufferfull),pointer,intent(in)::callbackonflush !< fonction a appeler lorsque le contenu du buffer doit etre obilgatoirement propage
       procedure(t_procgraphdot), pointer :: pfcndot !< fonction a appeler pour generer le graphique
       
       pfcndot => finalconsumer_ongraphdot
       allocate(this%m_input_buffer)
       call this%m_input_buffer%set_owner(this, callbackonfull,             &
     &         callbackonflush,pfcndot)
       
      end subroutine  finalconsumer_setconsumer

!***********************************************************************
!> @brief retourne le consommateur du buffer d'entree
!***********************************************************************
            function finalconsumer_getconsumer(this) result(cs)
             implicit none
             class(t_finalconsumer), intent(inout):: this  !< dummy argument
#if __INTEL_COMPILER <= 1600
             type(t_buffer_consumer), pointer :: cs 
#else
             class(t_buffer_consumer), pointer :: cs 
#endif
             
              cs => this%m_input_buffer
            end  function finalconsumer_getconsumer

!***********************************************************************
!> @brief  fonction appellee pour generer le graph au format dot
!***********************************************************************
       subroutine finalconsumer_ongraphdot(this,userdata, nodesrc,          &
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
           class is(t_finalconsumer)
             call buffer_dot_writelabel(num_file,curnode,                 &
     &          trim(userdata%m_dotname))
             call buffer_dot_writenode(num_file, nodesrc,curnode)
             nodecounter = curnode
           class default
            stop 'finalconsumer_ongraphdot : bad class'
          end select 

       end subroutine finalconsumer_ongraphdot

      end module mod_finalconsumer

