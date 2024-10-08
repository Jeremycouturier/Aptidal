!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ci_base.f 
!!  \brief representation abstraite des conditions initiales
!!
!!     Cette classe t_ci_base est utilisee comme classe de base 
!!     pour toutes les autres classes.
!!
! history : creation 14/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour stocker les conditions initiales des planetes
!***********************************************************************
      module mod_ci_base
       
!***********************************************************************
!> @class t_ci_base
!! classe de base pour le stockage des conditions initiales des planetes
!!      contient seulement son identifiant
!!  
!***********************************************************************
      type :: t_ci_base  !, abstract

          character(len=64) :: m_id !< identifiant 
          

       contains
                    
          procedure :: set_id => ci_base_set_id ! fixe l'identifiant

          procedure :: debug => ci_base_debug ! affiche les informations de debug

      end type t_ci_base  

!***********************************************************************
!> @class t_ci_base_ptr
!! type pour les tableaux de pointeurs de conditions initiales
!***********************************************************************
      type :: t_ci_base_ptr 
          class(t_ci_base),pointer :: m_ptr => NULL() !< pointeur vers une seule condition initiale
      end type t_ci_base_ptr 

           abstract interface

!***********************************************************************
! interface pour la lecture d'une ligne
!***********************************************************************
            subroutine destructor(this)
             import t_ci_base
             class(t_ci_base), intent(inout) :: this !< destructeur
            end subroutine destructor

          end interface


       contains

!***********************************************************************
!> @brief fixe l'identifiant de la condition initiale
!***********************************************************************
            subroutine ci_base_set_id(this, id)
             implicit none
             class(t_ci_base), intent(inout):: this  !< conditions initiales
             character(len=*), intent(in) :: id  !< nom de l'identifiant
             this%m_id = id
            end  subroutine ci_base_set_id
     
!***********************************************************************
!> @brief affiche le debug de la structure 
!***********************************************************************
      subroutine ci_base_debug(this)
       implicit none
       class(t_ci_base), intent(in):: this    !< conteneur de conditions initiales
       
         write(*,*) "ci_base : id", this%m_id
     
      end  subroutine ci_base_debug  

      end module mod_ci_base
