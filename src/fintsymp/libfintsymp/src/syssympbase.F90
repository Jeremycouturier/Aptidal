!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file syssympbase.f 
!!  \brief systeme pouvant etre integre par un schema symplectique
!!
!!     cette classe t_syssymp est une classe abstraite dont doit deriver 
!!     les systemes integrees avec intsympbase
!!
!!
! history : creation 20/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme integre avec un integrateur symplectique
!***********************************************************************
      module mod_syssympbase
             use mod_arret
      
!***********************************************************************
!> @class t_syssymp
!! classe de base decrivant 1 systeme integre avec un integrateur symplectique
!!  
!***********************************************************************
      type, abstract :: t_syssymp 
      
          logical :: sys_pas_fixe = .true. ! systeme a pas fixe
          
       contains         
          
          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          
          procedure(create_output), deferred :: create_output ! creation d'une copie pour effectuer les pas de sortie
          procedure(copy_output), deferred :: copy_output   ! copie de this dans la copie de sortie
          procedure(idelete_output), deferred , NOPASS:: delete_output ! destruction de l'objet cree par create_output
          
          
          procedure(get_outputsteps), deferred  :: get_outputsteps ! fourni les pas de sorties  (appellee par l'integrateur)
          procedure(iwriteoutput), deferred  :: write_output    ! ecriture de la sortie (appellee par l'integrateur)
          procedure(pflushoutput), deferred  :: flush_output    ! forcer a transmettre les buffers a la fin de l'integration (appellee par l'integrateur)
          procedure(copy_output_feedback), deferred ::                  &
     &            copy_output_feedback ! copie de la copie de sortie vers this des informations necessaires 


          procedure(t_pasA), deferred :: pasA
          procedure(t_pasA), deferred :: pasB
          procedure(t_pasA), deferred :: pasC
          
          procedure(t_set_error_time), deferred :: set_error_time  ! fixe le temps pour le controle des erreurs
          procedure(t_check_error), deferred :: check_error  ! controle si une erreur s'est produite
          procedure(t_get_error_time), deferred :: get_error_time  ! recupere le temps de l'erreur
          procedure(t_get_error_code), deferred :: get_error_code  ! recupere le code d'erreur
          procedure(t_get_error_body), deferred :: get_error_body  ! recupere le corps generant l'erreur
          procedure(t_get_error), deferred ::  get_error ! recupere l'erreur
          procedure(t_set_error), deferred ::  set_error ! recupere l'erreur

          procedure(t_graphdot), deferred :: graphdot  ! procedure appellee pour generer le graphique dot
         
          procedure :: get_adaptativetime =>  syssympbase_get_adaptativetime ! retourne le temps dans le cas d'un schemea adaptatif

          procedure :: dumpbin =>  syssympbase_dumpbin ! dump binaire vers un fichier
          procedure :: restorebin => syssympbase_restorebin ! restauration binaire depuis un fichier

      end type t_syssymp  

           abstract interface
!***********************************************************************
! type de fonction pour realiser un pas d'integration A, B, ou C
!***********************************************************************
            subroutine t_pasA(this,mdt)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this
             real(TREAL), intent(in) :: mdt 
            end subroutine t_pasA

!***********************************************************************
! ecriture dans le buffer de sortie
!***********************************************************************
            subroutine iwriteoutput(this,it,t)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this !< dummy argument
             integer(8), intent(in) :: it !< iteration actuelle
             real(TREAL), intent(in) :: t  !< temps correspondant a it
            end subroutine iwriteoutput

!***********************************************************************
! force le flush du buffer de sortie
!***********************************************************************
            subroutine pflushoutput(this)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this !< dummy argument
             end subroutine pflushoutput

!***********************************************************************
!> @brief creation d'une copie pour effectuer les pas de sortie
!! en general, duplique simplement this
!***********************************************************************
            subroutine create_output(this, sysout)
             import t_syssymp
             implicit none
             class(t_syssymp), intent(inout):: this  !< dummy argument
             class(t_syssymp), allocatable, intent(out) :: sysout  !< systeme cree pour la sortie
             end  subroutine create_output
        
!***********************************************************************
!> @brief destruction de l'objet cree par create_output
!! en general, detruit  this
!***********************************************************************
            subroutine idelete_output(this)
             import t_syssymp
             implicit none
             class(t_syssymp), allocatable, intent(inout):: this  !< dummy argument
            end  subroutine idelete_output
      
!***********************************************************************
!> @brief copie de this dans la copie de sortie (sysout)
!! en general, sysout = this 
!***********************************************************************
            subroutine copy_output(this, sysout)
             import t_syssymp
             implicit none
             class(t_syssymp), intent(inout):: this  !< dummy argument
             class(t_syssymp), intent(inout) :: sysout  !< systeme cree pour la sortie
            end  subroutine copy_output

!***********************************************************************
!> @brief copie des informations fournis en feedback des buffers de  this dans dst
!***********************************************************************
            subroutine copy_output_feedback(this, dst)
             import t_syssymp
             implicit none
             class(t_syssymp), intent(inout):: this  !< dummy argument
             class(t_syssymp), intent(inout) :: dst  !< destination (meme type que this)
            end  subroutine copy_output_feedback

!***********************************************************************
!> @brief retourne les pas de sortie ou l'integrateur doit appeller la fonction 
!!        write_output
!***********************************************************************
            subroutine get_outputsteps(this, nstepscall, arstepscall)
             import t_syssymp
             implicit none
             class(t_syssymp), intent(inout):: this  !< dummy argument
             integer, intent(out) :: nstepscall  !< taille de arstepscall
             integer(8), dimension(:), allocatable, intent(out) ::       &
     &          arstepscall !< tableau des pas de sortie
            end  subroutine get_outputsteps

!***********************************************************************
! type de fonction pour controler si une erreur s'est produite
!***********************************************************************
           function t_check_error(this) result(ret)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             logical ret ! code de retour
            end function t_check_error

!***********************************************************************
! type de fonction pour recuperer une erreur
!***********************************************************************
           function t_get_error(this) result(arret)
             import t_syssymp
             import t_arret
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             type(t_arret) :: arret ! code de retour
            end function t_get_error

!***********************************************************************
! type de fonction pour fixer une erreur
!***********************************************************************
           subroutine t_set_error(this,arret)
             import t_syssymp
             import t_arret
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             type(t_arret), intent(in) :: arret ! code de retour
            end subroutine t_set_error

!***********************************************************************
!> type de fonction pour fixer le temps de l'erreur dans les donnees
!! indique le comportement du programme qui l'appelle par la suite
!! @return retourne 0 si le programme doit se poursuivre \n
!! retourne 1 si on doit s'arreter d'integrer (mais poursuite de l'execution)
!***********************************************************************
            function t_set_error_time(this,t) result(ret)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             real(TREAL), intent(in) :: t !< temps de l'erreur 
             integer ret ! code de retour
            end function t_set_error_time

!***********************************************************************
!> type de fonction pour recuperer le temps de l'erreur dans les donnees
!***********************************************************************
            function t_get_error_time(this) result(t)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             real(TREAL) :: t ! temps de l'erreur 
            end function t_get_error_time

!***********************************************************************
!> type de fonction pour recuperer le corps ayant genere l'erreur dans les donnees
!***********************************************************************
            function t_get_error_body(this) result(body)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             integer :: body ! corps
            end function t_get_error_body

!***********************************************************************
!> type de fonction pour recuperer le code de l'erreur
!***********************************************************************
            function t_get_error_code(this) result(code)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             integer :: code ! code
            end function t_get_error_code

!***********************************************************************
!>  type de fonction pour generer le graph au format dot
!***********************************************************************
       subroutine t_graphdot(this, nodesrc, nodecounter,  num_file)
             import t_syssymp
             class(t_syssymp), intent(inout) :: this  !< dummy argument
             integer, intent(in) :: nodesrc !< numero du noeud  source 
             integer, intent(inout) :: nodecounter !< numero du noeud pour les creations (en sortie, contient le dernier numero utilise)
             integer, intent(in) :: num_file !<numero du fichier de sortie
        
       end subroutine t_graphdot
 

          end interface

     
      contains
         

!***********************************************************************
!> @brief retourne le temps du systeme dans le cas d'un schemea adaptatif
!***********************************************************************
      function syssympbase_get_adaptativetime(this) result(t)
       implicit none
       class(t_syssymp), intent(in):: this  !< dummy argument
       
       real(TREAL) ::t ! temps du systeme
       t  = 0.
       write(*,*) 'systeme nutilisant pas un temps adaptatif'
       stop 1
         
      end function  syssympbase_get_adaptativetime  

!***********************************************************************
!> @brief dump binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (in) objet a sauvegarder
!***********************************************************************
      subroutine syssympbase_dumpbin(this, file) 
       use mod_io_dump
       implicit none
       class(t_syssymp), intent(in):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       write(file%m_nf) this%sys_pas_fixe

      end subroutine  syssympbase_dumpbin  

!***********************************************************************
!> @brief restauration binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (inout) objet a sauvegarder
!***********************************************************************
      subroutine syssympbase_restorebin(this, file) 
       use mod_io_dump
       implicit none
       class(t_syssymp), intent(inout):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       read(file%m_nf) this%sys_pas_fixe

      end subroutine  syssympbase_restorebin  

      end module mod_syssympbase
