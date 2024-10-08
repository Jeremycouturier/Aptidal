!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciwrite_base.f 
!!  \brief ecriture d'un fichier texte de conditions initiales quelconques basees sur un stockage par ligne
!!
!!
! history : creation 25/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture d'un fichier texte de conditions initiales quelconques
!***********************************************************************
      module mod_ciwrite_base
       use mod_ci_base
!***********************************************************************
!> @class t_ciwrite_base
!! classe d'ecriture d'un fichier texte de conditions initiales quelconques
!!  
!***********************************************************************
      type, abstract :: t_ciwrite_base
          
          integer, private :: m_nf !< numero fortran du fichier 
          character(len=512), private :: m_name !< nom du fichier

       contains
          
          
          procedure :: set_name => ciwrite_base_set_name ! fixe le nom du fichier et le numero de fichier
          procedure :: get_nf => ciwrite_base_get_nf ! recupere le numero de fichier fortran 
          procedure :: fopen => ciwrite_base_fopen ! ouverture d'un fichier de conditions initiales
          procedure :: fclose => ciwrite_base_fclose ! fermeture d'un fichier de conditions initiales
          procedure :: fwrite => ciwrite_base_fwrite ! ecriture de toutes les conditions initiales issues d'un reader


          ! les procedures suivantes doivent etre surchargees par les types derivant de cette classe
          procedure(fwriteline), deferred :: fwriteline  ! ecriture d'une seule ligne de conditions initiales
          
      end type t_ciwrite_base  

           abstract interface

            ! interface pour l'ecriture d'une ligne
            subroutine fwriteline(this, ci)
             import t_ciwrite_base
             import t_ci_base
             class(t_ciwrite_base), intent(inout) :: this !< lecteur du fichier de ci
             class(t_ci_base),  intent(in) :: ci !<  condition initiale a ecrire
            end subroutine fwriteline


          end interface
    
!***********************************************************************
!> @class t_ciwrite_base_ptr
!! type pour les tableaux de pointeurs d'ecrivains
!***********************************************************************
      type :: t_ciwrite_base_ptr 
          class(t_ciwrite_base),pointer :: m_ptr => NULL() !< pointeur vers un seul ecrivain
      end type t_ciwrite_base_ptr 

      contains
      
!***********************************************************************
!> @brief fixe le nom du fichier et le numero fortran
!***********************************************************************
            subroutine ciwrite_base_set_name(this, name, nf)
             implicit none
             class(t_ciwrite_base), intent(inout):: this  !< ecrivain de ci
             character(len=*), intent(in) :: name  !< nom du fichier
             integer, intent(in) :: nf  !< numero fortran du fichier
             this%m_name = name
             this%m_nf = nf
            end  subroutine ciwrite_base_set_name

!***********************************************************************
!> @brief ouvre le fichier en ecriture
!***********************************************************************
            subroutine ciwrite_base_fopen(this)
             implicit none
             class(t_ciwrite_base), intent(inout):: this  !< ecrivain de ci
             open(unit=this%m_nf,form="formatted",file=this%m_name,      &
     &        status='unknown')
            end  subroutine ciwrite_base_fopen

!***********************************************************************
!> @brief ferme le fichier
!***********************************************************************
            subroutine ciwrite_base_fclose(this)
             implicit none
             class(t_ciwrite_base), intent(inout):: this  !< ecrivain de ci
             close(this%m_nf)
            end  subroutine ciwrite_base_fclose
    
!***********************************************************************
!> @brief recupere le numero de fichier fortran.
!! Cette fonction ne doit etre utilisee que par les types derivees de ce type this.
!***********************************************************************
            function ciwrite_base_get_nf(this) result(nf)
             implicit none
             class(t_ciwrite_base), intent(inout):: this  !< ecrivain de ci
             integer :: nf
             nf = this%m_nf
            end  function ciwrite_base_get_nf

!***********************************************************************
!> @brief ouvre le fichier du lecteur, lit toutes les conditions intiiales via le lecteur puis ferme le fichier du lecteur.
!! Chaque condition initiale est ensuite ecrite par l'ecrivain en ouvrant le fichier puis le fermant a la fin.
!***********************************************************************
            subroutine ciwrite_base_fwrite(this, reader) 
             use mod_ciread_base
             implicit none
             class(t_ciwrite_base), intent(inout):: this  !< ecrivain de ci
             class(t_ciread_base), intent(inout):: reader  !< lecteur de ci
             
             class(t_ci_base), pointer :: ligneci ! condition initiale lue

             call reader%fopen()
             call this%fopen()
             
             do while(reader%freadline(ligneci))
               call this%fwriteline(ligneci)
               deallocate(ligneci)
             enddo
             
             call reader%fclose()
             call this%fclose()
             
            end  subroutine ciwrite_base_fwrite
      
      end module mod_ciwrite_base
