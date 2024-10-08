!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Boue
! 
!>  \file ci_circummaree.f 
!!  \brief representation des conditions initiales d'un systemes planetaires
!!      contenant seulement ses masses et ses coordonnees initiales (6 composantes par corps)
!!      plus les parametres physiques de la planete interne : rayon equatorial, moment d'inertie,
!!      second nombre de Love, temps de relaxation (5 parametres) et le vecteur instantane de 
!!      rotation (3 composantes)
!!      attention:  le corps central est en position 0 pour les masses
!!
!!     Cette classe t_ci_circummaree est utilisee pour stockee les informations lues par t_ciread_circummaree.
!!     Elle etend la classe t_ci_planewt
!!
! history : creation 14/12/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour stocker les conditions initiales des planetes 
! ainsi que le moment cinetique de la premiere planete
!***********************************************************************
      module mod_ci_circummaree
       use mod_ci_planewt
!***********************************************************************
!> @class t_ci_circummaree
!! classe de stockage des conditions initiales des planetes
!!      contenant seulement ses masses et ses coordonnees initiales (6 composantes par corps)
!!      plus le moment cinetique de rotation initial de la planete
!!  
!***********************************************************************
      type, extends(t_ci_planewt) :: t_ci_circummaree

          real(TREAL) :: m_Rpla !< rayon equatorial de la planete interne
          real(TREAL) :: m_xi   !< moment d'inertie de la planete interne
          real(TREAL) :: m_k20  !< second nombre de Love fluide de la planete interne
          real(TREAL) :: m_taue !< temps de relaxation elastique de Maxwell de la planete interne
          real(TREAL) :: m_tau2 !< temps de relaxation global de la planete interne
          real(TREAL), dimension(3) :: m_vrot_ci !< vecteur instantane de rotation de la planete interne

       contains

          procedure:: debug => ci_circummaree_debug ! debug de la structure
          
      end type t_ci_circummaree

      contains

!***********************************************************************
!> @brief affiche le debug de la structure 
!***********************************************************************
      subroutine ci_circummaree_debug(this)
       implicit none
       class(t_ci_circummaree), intent(in):: this    !< conteneur de conditions initiales
       
         integer j
         write(*,*) "ci_circummaree_debug - start"
         call ci_planewt_debug(this)
         write(*,*) "Rpla    :", this%m_Rpla
         write(*,*) "xi      :", this%m_xi
         write(*,*) "k20     :", this%m_k20
         write(*,*) "taue    :", this%m_taue
         write(*,*) "tau2    :", this%m_tau2
         write(*,*) "vrot    :", (this%m_vrot_ci(j),j=1,3)
         write(*,*) "-----------------------"
           
      end  subroutine ci_circummaree_debug  


      end module mod_ci_circummaree
