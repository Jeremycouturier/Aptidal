!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file ciread_pladistmin.f 
!!  \brief lecture d'un fichier texte de distance minimale des planetes
!!
!!     Cette classe t_ciread_pladistmin gere la lecture du fichier de distance minimale 
!!     pour des planetes (1 composante par planete) .
!!     Ce fichier de conditions initiales contient :
!!     Chaque ligne du fichier est composee de :
!!             ch  nplan distmin_1 ....distmin_nplan
!!
!!       avec :
!!             ch : chaine de caracteres (sans espaces) qui peut servir d'identifiant pour la condition initiale
!!                  ne doit contenir au plus que 63 caracteres
!!                  e.g., N002 ou 002 ou 2 ou CERES0001
!!
!!             nplan : nombre de corps =0:nplan  (donc nplan+1 en realite)
!!             ecrit sur 2 caracteres
!!             distmin_1 ,..., distmin_nplan : distance minimale des planetes (nplan masses)
!!
! history : creation 29/01/2018
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la lecture d'un fichier de conditions initiales ne contenant que des masses et coordonnees
!***********************************************************************
      module mod_ciread_pladistmin
       use mod_ciread_base
       use mod_ci_pladistmin
       
!***********************************************************************
!> @class t_ciread_pladistmin
!! classe de lecture d'un fichier texte de conditions initiales des planetes
!! comprenant masses + coordonnees initiales
!!  
!***********************************************************************
      type, extends(t_ciread_base)  :: t_ciread_pladistmin
          
       contains
                    
          procedure :: freadline => ciread_pladistmin_freadline ! lecture d'une seule ligne de conditions initiales         
          
      end type t_ciread_pladistmin  

     
      contains

!***********************************************************************
!> @brief lecture d'une ligne du fichier de distances minimales (1 composantes par planetes)
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a lire.
!! @return retourne .true. en cas de succes, retourne .false. en cas d'eof
!***********************************************************************
            function ciread_pladistmin_freadline(this, ci) result(ret)
             implicit none
             class(t_ciread_pladistmin), intent(inout) :: this !< lecteur du fichier de ci
             class(t_ci_base), pointer, intent(out) :: ci !< condition initiale lue en sortie
             logical :: ret !< =true => condition lue, =false => aucune condition lue
             
             class(t_ci_pladistmin), pointer :: cipla
             character(len=64) :: id
             character(len=4096) :: ligne
             integer :: nbplan
             integer k
             
             ! allocation de la condition initiale
             allocate(cipla)
             
             ! lecture de la ligne sous forme de chaine
             ret = this%fread_strline(ligne)
             if (ret) then
              !write(*,*) 'ligne=',trim(ligne)
             
              ! lecture du nombre de corps
              read(ligne,*) id, nbplan 
              write(*,*) 'id pladistmin=',id
              write(*,*) 'nbplan pladistmin=',nbplan
             
              call cipla%set_id(id)
              call cipla%set_plan_nb(nbplan)
                          
              ! lecture des distances minimales des planetes
              read(ligne,*) id,nbplan,(cipla%m_plan_dmin(k),k=1,nbplan)  

             
              ci => cipla
             else
               deallocate(cipla)
               ci =>NULL()
             endif
             
            end function ciread_pladistmin_freadline

      end module mod_ciread_pladistmin
