!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ci_partnewt.f 
!!  \brief representation des conditions initiales d'une particule
!!      contenant seulement ses coordonnees initiales (6 composantes)
!!
!!     Cette classe t_ci_partnewt est utilisee pour stockee les informations lues par t_ciread_partnewt.
!!
! history : creation 17/06/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour stocker les conditions initiales des particules
!***********************************************************************
      module mod_ci_partnewt
       use mod_ci_base
!***********************************************************************
!> @class t_ci_partnewt
!! classe de stockage des conditions initiales de la particule
!!      contenant seulement ses coordonnees initiales (6 composantes)
!!  
!***********************************************************************
      type, extends(t_ci_base) :: t_ci_partnewt

          integer :: m_ci_type  !< type de conditions initiales
          
          real(TREAL), dimension(1:6) :: m_plan_coord_ci !< coordonnees initiales de la particule

       contains
                    
          procedure:: set_plan_ci_type => ci_partnewt_set_plan_ci_type ! fixe le type de conditions initiales

          procedure:: debug => ci_partnewt_debug ! debug de la structure
          
      end type t_ci_partnewt  

      contains
          
!***********************************************************************
!> @brief fixe le type de conditions initiales
!***********************************************************************
      subroutine ci_partnewt_set_plan_ci_type(this, ci_type)
       implicit none
       class(t_ci_partnewt), intent(inout):: this    !< conteneur de conditions initiales
       integer, intent(in) :: ci_type               !< type de conditions initiales
       
       this%m_ci_type = ci_type
       
      end  subroutine ci_partnewt_set_plan_ci_type  

!***********************************************************************
!> @brief affiche le debug de la structure 
!***********************************************************************
      subroutine ci_partnewt_debug(this)
       implicit none
       class(t_ci_partnewt), intent(in):: this    !< conteneur de conditions initiales
       
         integer j
         write(*,*) "ci_partnewt_debug - start"
         call ci_base_debug(this)
         write(*,*) "m_ci_type", this%m_ci_type
         write(*,*) "ci :" 
         write(*,*) (this%m_plan_coord_ci(j),j=1,6)
         write(*,*) "-----------------------"
           
      end  subroutine ci_partnewt_debug  


      end module mod_ci_partnewt
