!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciwrite_planewt.f 
!!  \brief ecriture d'un fichier texte de conditions initiales des planetes
!!
!!     Cette classe t_ciwrite_planewt gere l'ecriture du fichier de conditions initiales
!!     pour des planetes (masse + coordonnees initiales seulement (6 composantes par corps) ).
!!     
!!     Le format du fichier est decrit dans ciread_planewt.f
!!
! history : creation 25/03/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture d'un fichier de conditions initiales ne contenant que des masses et coordonnees
!***********************************************************************
      module mod_ciwrite_planewt
       use mod_ciwrite_base
       use mod_ci_planewt
       
!***********************************************************************
!> @class t_ciwrite_planewt
!! classe d'ecriture d'un fichier texte de conditions initiales des planetes
!! comprenant masses + coordonnees initiales
!!  
!***********************************************************************
      type, extends(t_ciwrite_base)  :: t_ciwrite_planewt
          
       contains
                    
          procedure :: fwriteline => ciwrite_planewt_fwriteline ! ecriture d'une seule ligne de conditions initiales         
          
      end type t_ciwrite_planewt  

     
      contains

!***********************************************************************
!> @brief ecriture d'une ligne du fichier de conditions initiales contenant uniquement masses + ci (6 composantes par corps)
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a ecrire.
!***********************************************************************
       subroutine ciwrite_planewt_fwriteline(this, ci) 
        implicit none
        class(t_ciwrite_planewt), intent(inout) :: this !< lecteur du fichier de ci
        class(t_ci_base), intent(in) :: ci !< condition initiale a ecrire
        
        integer :: nbplan
        integer j,k
        character(len=1002) fmt1002
             
        select type(ci)
         class is(t_ci_planewt)
          nbplan = ci%m_plan_nb
          ! ecriture des conditions initiales
#if TREAL == 8 
          write (fmt1002,'(A,I4,A,I4,A)') "(A,X,I3,",nbplan+1,           &
     &    "(1X,D23.16),X,I2,",6*nbplan,"(1X,D23.16))"
#else
          write (fmt1002,'(A,I4,A,I4,A)') "(A,X,I3,",nbplan+1,           &
     &    "(1X,D33.26),X,I2,",6*nbplan,"(1X,D33.26))"
#endif
          write(this%get_nf(),fmt1002) trim(ci%m_id), nbplan,            &   
     &        (ci%m_plan_mass(k),k=0,nbplan), ci%m_ci_type,              &   
     &        ((ci%m_plan_coord_ci(j,k),j=1,6),k=1,nbplan)     
         
!#if TREAL == 8 
!1002  format(A,X,I3,<nbplan+1>(1X,D23.16),X,I2,<6*nbplan>(1X,D23.16))       !(double)
!#else
!1002  format(A,X,I3,<nbplan+1>(1X,D33.26),X,I2,<6*nbplan>(1X,D33.26))       !(etendu)
!#endif
         class default
          stop 'ligneci in ciwrite_planewt_fwriteline : bad class'
        end select 

       end subroutine ciwrite_planewt_fwriteline

      end module mod_ciwrite_planewt
