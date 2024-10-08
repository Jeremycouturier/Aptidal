!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciwrite_pladistmin.f 
!!  \brief ecriture d'un fichier texte de distances minimales des planetes
!!
!!     Cette classe t_ciwrite_pladistmin gere l'ecriture du fichier des distances minimales
!!     pour des planetes  (1 composante par corps) .
!!     
!!     Le format du fichier est decrit dans ciread_pladistmin.f
!!
! history : creation 25/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture d'un fichier  de distances minimales des planetes
!***********************************************************************
      module mod_ciwrite_pladistmin
       use mod_ciwrite_base
       use mod_ci_pladistmin
       
!***********************************************************************
!> @class t_ciwrite_pladistmin
!! classe d'ecriture d'un fichier texte de distances minimales des planetes
!!  
!***********************************************************************
      type, extends(t_ciwrite_base)  :: t_ciwrite_pladistmin
          
       contains
                    
          procedure :: fwriteline => ciwrite_pladistmin_fwriteline ! ecriture d'une seule ligne de conditions initiales         
          
      end type t_ciwrite_pladistmin  

     
      contains

!***********************************************************************
!> @brief ecriture d'une ligne du fichier de distances minimales (1 composantes par planetes)
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a ecrire.
!***********************************************************************
       subroutine ciwrite_pladistmin_fwriteline(this, ci) 
        implicit none
        class(t_ciwrite_pladistmin), intent(inout) :: this !< lecteur du fichier de ci
        class(t_ci_base), intent(in) :: ci !< condition initiale a ecrire
        
        integer :: nbplan
        integer k
        character(len=1002) fmt1002
             
        select type(ci)
         class is(t_ci_pladistmin)
          nbplan = ci%m_plan_nb
          ! ecriture des conditions initiales
#if TREAL == 8 
          write (fmt1002,'(A,I4,A,I4,A)') "(A,X,I3,",nbplan,             &
     &    "(1X,D23.16))"
#else
          write (fmt1002,'(A,I4,A,I4,A)') "(A,X,I3,",nbplan,             &
     &    "(1X,D33.26))"
#endif
          write(this%get_nf(),fmt1002) trim(ci%m_id), nbplan,            &   
     &        (ci%m_plan_dmin(k),k=1,nbplan)     
         
         class default
          stop 'ligneci in ciwrite_pladistmin_fwriteline : bad class'
        end select 

       end subroutine ciwrite_pladistmin_fwriteline

      end module mod_ciwrite_pladistmin
