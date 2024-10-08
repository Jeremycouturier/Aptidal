!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ci_planewt.f 
!!  \brief representation des conditions initiales d'un systeme planetaire
!!      contenant seulement ses masses et ses coordonnees initiales (6 composantes par corps)
!!      attention:  le corps central est en position 0 pour les masses
!!
!!     Cette classe t_ci_planewt est utilisee pour stockee les informations lues par ti_ciread_planewt.
!!
! history : creation 14/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour stocker les conditions initiales des planetes
!***********************************************************************
      module mod_ci_planewt
       use mod_ci_base
!***********************************************************************
!> @class t_ci_planewt
!! classe de stockage des conditions initiales des planetes
!!      contenant seulement ses masses et ses coordonnees initiales (6 composantes par corps)
!!  
!***********************************************************************
      type, extends(t_ci_base) :: t_ci_planewt

          integer :: m_plan_nb  !< nombre de planetes (sans tenir compte du corps central)
          integer :: m_ci_type  !< type de conditions initiales
          
          real(TREAL), dimension(:), allocatable :: m_plan_mass !< masse des planetes(0:m_plan_nb)  (commence a 0: corps central)
          real(TREAL), dimension(:,:), allocatable :: m_plan_coord_ci !< coordonnees initiales des planetes(1:6,1:m_plan_nb)

       contains
                    
          procedure:: set_plan_nb => ci_planewt_set_plan_nb ! fixe le nombre de planetes 
          procedure:: set_plan_ci_type => ci_planewt_set_plan_ci_type ! fixe le type de conditions initiales

          procedure:: debug => ci_planewt_debug ! debug de la structure
          
      end type t_ci_planewt  

      contains
      
!***********************************************************************
!> @brief fixe le nombre de planetes et alloue la memoire
!***********************************************************************
      subroutine ci_planewt_set_plan_nb(this, nplan)
       implicit none
       class(t_ci_planewt), intent(inout):: this    !< conteneur de conditions initiales
       integer, intent(in) :: nplan                 !< nombre de planetes (sans le corps central)
       
       this%m_plan_nb = nplan
       
       allocate(this%m_plan_mass(0:nplan))
       allocate(this%m_plan_coord_ci(1:6,1:nplan))
       
      end  subroutine ci_planewt_set_plan_nb  
          
!***********************************************************************
!> @brief fixe le type de conditions initiales
!***********************************************************************
      subroutine ci_planewt_set_plan_ci_type(this, ci_type)
       implicit none
       class(t_ci_planewt), intent(inout):: this    !< conteneur de conditions initiales
       integer, intent(in) :: ci_type               !< type de conditions initiales
       
       this%m_ci_type = ci_type
       
      end  subroutine ci_planewt_set_plan_ci_type  

!***********************************************************************
!> @brief affiche le debug de la structure 
!***********************************************************************
      subroutine ci_planewt_debug(this)
       implicit none
       class(t_ci_planewt), intent(in):: this    !< conteneur de conditions initiales
       
         integer k,j
         write(*,*) "ci_planewt_debug - start"
         call ci_base_debug(this)
         write(*,*) "m_plan_nb", this%m_plan_nb
         write(*,*) "masses :"
         do k=0, this%m_plan_nb
          write(*,*) k," : ",  this%m_plan_mass(k)
         enddo
         write(*,*) "m_ci_type", this%m_ci_type
         write(*,*) "ci :" 
         do k=1, this%m_plan_nb
          write(*,*) k," : ",  (this%m_plan_coord_ci(j,k),j=1,6)
         enddo
         write(*,*) "-----------------------"
           
      end  subroutine ci_planewt_debug  


      end module mod_ci_planewt
