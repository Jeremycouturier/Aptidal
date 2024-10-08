!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciread_partnewt.f 
!!  \brief lecture d'un fichier texte de conditions initiales des particules
!!
!!     Cette classe t_ciread_partnewt gere la lecture du fichier de conditions initiales
!!     pour des particules (coordonnees initiales seulement (6 composantes par corps) ).
!!     Ce fichier de conditions initiales contient :
!!     Chaque ligne du fichier est composee de :
!!             ch  ci_type  ci_part_1[1:6]
!!
!!       avec :
!!             ch : chaine de caracteres (sans espaces) qui peut servir d'identifiant pour la condition initiale
!!                  ne doit contenir au plus que 63 caracteres
!!                  e.g., N002 ou 002 ou 2 ou CERES0001
!!
!!             ci_type :  type des coordonnees des particules
!!             ci_part_1(1:6) :   coordonnees (6 composantes) pour chaque particule
!!
! history : creation 17/06/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la lecture d'un fichier de conditions initiales ne contenant que des  coordonnees
!***********************************************************************
      module mod_ciread_partnewt
       use mod_ciread_base
       use mod_ci_partnewt
       
!***********************************************************************
!> @class t_ciread_partnewt
!! classe de lecture d'un fichier texte de conditions initiales des particules
!! comprenant coordonnees initiales
!!  
!***********************************************************************
      type, extends(t_ciread_base)  :: t_ciread_partnewt
          
       contains
                    
          procedure :: freadline => t_ciread_partnewt_freadline ! lecture d'une seule ligne de conditions initiales         
          
      end type t_ciread_partnewt  

     
      contains

!***********************************************************************
!> @brief lecture d'une ligne du fichier de conditions initiales contenant uniquement ci (6 composantes par corps)
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a lire.
!! @return retourne .true. en cas de succes, retourne .false. en cas d'eof
!***********************************************************************
            function t_ciread_partnewt_freadline(this, ci) result(ret)
             implicit none
             class(t_ciread_partnewt), intent(inout) :: this !< lecteur du fichier de ci
             class(t_ci_base), pointer, intent(out) :: ci !< condition initiale lue en sortie
             logical :: ret !< =true => condition lue, =false => aucune condition lue
             
             class(t_ci_partnewt), pointer :: cipart
             character(len=64) :: id
             character(len=4096) :: ligne
             integer :: ci_type
             integer j
             
             ! allocation de la condition initiale
             allocate(cipart)
             
             ! lecture de la ligne sous forme de chaine
             ret = this%fread_strline(ligne)
             if (ret) then
             !write(*,*) 'ligne=',ligne
             
              ! lecture du nombre de corps
              read(ligne,*) id 
              write(*,*) 'id=',id
             
              call cipart%set_id(id)
                          
             ! lecture du type de conditions initiales
            read(ligne,*)id,ci_type,(cipart%m_plan_coord_ci(j),j=1,6)     

              call cipart%set_plan_ci_type(ci_type)
             
              ci => cipart
             else
              deallocate(cipart)
              ci =>NULL()
             endif
             
            end function t_ciread_partnewt_freadline

      end module mod_ciread_partnewt
