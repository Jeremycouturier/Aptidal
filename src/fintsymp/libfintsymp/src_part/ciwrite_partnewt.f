!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciwrite_partnewt.f 
!!  \brief ecriture d'un fichier texte de conditions initiales des planetes
!!
!!     Cette classe t_ciwrite_partnewt gere l'ecriture du fichier de conditions initiales
!!     pour des particules (coordonnees initiales seulement (6 composantes par corps) ).
!!     
!!     Le format du fichier est decrit dans ciread_partnewt.f
!!
! history : creation 16/09/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture d'un fichier de conditions initiales ne contenant que des masses et coordonnees
!***********************************************************************
      module mod_ciwrite_partnewt
       use mod_ciwrite_base
       use mod_ci_partnewt
       
!***********************************************************************
!> @class t_ciwrite_partnewt
!! classe d'ecriture d'un fichier texte de conditions initiales des particules
!! comprenant les coordonnees initiales
!!  
!***********************************************************************
      type, extends(t_ciwrite_base)  :: t_ciwrite_partnewt
          
       contains
                    
          procedure :: fwriteline => ciwrite_partnewt_fwriteline ! ecriture d'une seule ligne de conditions initiales         
          
      end type t_ciwrite_partnewt  

     
      contains

!***********************************************************************
!> @brief ecriture d'une ligne du fichier de conditions initiales contenant uniquement les ci (6 composantes par corps)
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a ecrire.
!***********************************************************************
       subroutine ciwrite_partnewt_fwriteline(this, ci) 
        implicit none
        class(t_ciwrite_partnewt), intent(inout) :: this !< lecteur du fichier de ci
        class(t_ci_base), intent(in) :: ci !< condition initiale a ecrire
        
        integer j
             
        select type(ci)
         class is(t_ci_partnewt)
          ! ecriture des conditions initiales
          write(this%get_nf(),1002) trim(ci%m_id), ci%m_ci_type,          &   
     &        (ci%m_plan_coord_ci(j),j=1,6)     
         
#if TREAL == 8 
1002  format(A,X,I3,6(1X,D23.16))       !(double)
#else
1002  format(A,X,I3,6(1X,D33.26))       !(etendu)
#endif
         class default
          stop 'ligneci in ciwrite_partnewt_fwriteline : bad class'
        end select 

       end subroutine ciwrite_partnewt_fwriteline

      end module mod_ciwrite_partnewt
