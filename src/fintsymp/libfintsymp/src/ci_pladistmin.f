!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ci_pladistmin.f 
!!  \brief representation des distances minimales d'un systeme planetaire
!!      contenant seulement une distance minimale (1 composante par corps)
!!
!!     Cette classe t_ci_pladistmin est utilisee pour stockee les informations lues par t_ciread_pladistmin.
!!
! history : creation 29/01/2018
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour stocker les distances minimales d'un systeme planetaire
!***********************************************************************
      module mod_ci_pladistmin
       use mod_ci_base
!***********************************************************************
!> @class t_ci_pladistmin
!! classe de stockage des distances minimales d'un systeme planetaire
!!      contenant seulement les distances minimales (1 composante par corps)
!!  
!***********************************************************************
      type, extends(t_ci_base) :: t_ci_pladistmin

          integer :: m_plan_nb  !< nombre de planetes (sans tenir compte du corps central)
          
          real(TREAL), dimension(:), allocatable :: m_plan_dmin !<  distance miniale (1:m_plan_nb)  

       contains
                    
          procedure:: set_plan_nb => ci_pladistmin_set_plan_nb ! fixe le nombre de planetes 

          procedure:: debug => ci_pladistmin_debug ! debug de la structure
          
      end type t_ci_pladistmin  

      contains
      
!***********************************************************************
!> @brief fixe le nombre de planetes et alloue la memoire
!***********************************************************************
      subroutine ci_pladistmin_set_plan_nb(this, nplan)
       implicit none
       class(t_ci_pladistmin), intent(inout):: this    !< conteneur de conditions initiales
       integer, intent(in) :: nplan                 !< nombre de planetes (sans le corps central)
       
       this%m_plan_nb = nplan
       
       allocate(this%m_plan_dmin(1:nplan))
       
      end  subroutine ci_pladistmin_set_plan_nb  
          
!***********************************************************************
!> @brief affiche le debug de la structure 
!***********************************************************************
      subroutine ci_pladistmin_debug(this)
       implicit none
       class(t_ci_pladistmin), intent(in):: this    !< conteneur de conditions initiales
       
         integer k
         write(*,*) "ci_pladistmin_debug - start"
         call ci_base_debug(this)
         write(*,*) "m_plan_nb", this%m_plan_nb
         write(*,*) "distance minimale :"
         do k=1, this%m_plan_nb
          write(*,*) k," : ",  this%m_plan_dmin(k)
         enddo
         write(*,*) "-----------------------"
           
      end  subroutine ci_pladistmin_debug  


      end module mod_ci_pladistmin
