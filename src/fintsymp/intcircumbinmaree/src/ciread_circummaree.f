!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Boue
! 
!>  \file ciread_circummaree.f 
!!  \brief lecture d'un fichier texte de conditions initiales des planetes
!!
!!     Cette classe t_ciread_circummaree gere la lecture du fichier de conditions initiales
!!     pour des planetes (masse + coordonnees initiales seulement (6 composantes par corps) ),
!!     plus le moment cinetique de la premiere planete (3 composantes a la fin de la ligne
!!     Ce fichier de conditions initiales contient :
!!     Chaque ligne du fichier est composee de :
!!             ch  nplan m_0 m_1 ....m_nplan ci_type  ci_pla_1[1:6] .... ci_pla_nplan[1:6] ci_lcin[1:3]
!!
!!       avec :
!!             ch : chaine de caracteres (sans espaces) qui peut servir d'identifiant pour la condition initiale
!!                  ne doit contenir au plus que 63 caracteres
!!                  e.g., N002 ou 002 ou 2 ou CERES0001
!!
!!             nplan : nombre de corps =0:nplan  (donc nplan+1 en realite)
!!             ecrit sur 2 caracteres
!!             m_0,..., m_nplan : masse des corps (nplan+1 masses)
!!             ci_type :  type des coordonnees des planetes
!!             ci_pla_i(1:6) :   coordonnees (6 composantes) pour chaque planete i
!!             ci_lcin(1:3)  :   moment cinetique de la premiere planete (3 composantes) 
!!
! history : creation 14/12/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour la lecture d'un fichier de conditions initiales ne contenant que des masses et coordonnees
! et le moment cinetique de rotation de la premiere planete
!***********************************************************************
      module mod_ciread_circummaree
       use mod_ciread_base
       use mod_ci_circummaree
       
!***********************************************************************
!> @class t_ciread_circummaree
!! classe de lecture d'un fichier texte de conditions initiales des planetes
!! comprenant masses + coordonnees initiales
!!  
!***********************************************************************
      type, extends(t_ciread_base)  :: t_ciread_circummaree
          
       contains
                    
          procedure :: freadline => t_ciread_circummaree_freadline ! lecture d'une seule ligne de conditions initiales         
          
      end type t_ciread_circummaree

     
      contains

!***********************************************************************
!> @brief lecture d'une ligne du fichier de conditions initiales contenant uniquement masses + ci (6 composantes par corps)
!! + le moment cinetique de rotation de la premiere planete (3 composantes)
!! Elle suppose que le fichier ouvert est bien au debut de la ligne a lire.
!! @return retourne .true. en cas de succes, retourne .false. en cas d'eof
!***********************************************************************
            function t_ciread_circummaree_freadline(this, ci)            &
     &               result(ret)
             implicit none
             class(t_ciread_circummaree), intent(inout) :: this !< lecteur du fichier de ci
             class(t_ci_base), pointer, intent(out) :: ci !< condition initiale lue en sortie
             logical :: ret !< =true => condition lue, =false => aucune condition lue
             
             class(t_ci_circummaree), pointer :: cipla
             character(len=64) :: id
             character(len=4096) :: ligne
             integer :: nbplan, ci_type
             integer j,k
             
             ! allocation de la condition initiale
             allocate(cipla)
             
             ! lecture de la ligne sous forme de chaine
             ret = this%fread_strline(ligne)
             if (ret) then
             !write(*,*) 'ligne=',ligne
             
              ! lecture du nombre de corps
              read(ligne,*) id, nbplan 
              write(*,*) 'id=',id
              write(*,*) 'nbplan=',nbplan
             
              call cipla%set_id(id)
              call cipla%set_plan_nb(nbplan)
                          
              ! lecture du type de conditions initiales
              read(ligne,*) id, nbplan                                   &
     &          ,  (cipla%m_plan_mass(k),k=0,nbplan)                     &
     &          ,   ci_type                                              &
     &          , ((cipla%m_plan_coord_ci(j,k),j=1,6),k=1,nbplan)        &
     &          ,   cipla%m_Rpla                                         &
     &          ,   cipla%m_xi                                           &
     &          ,   cipla%m_k20                                          &
     &          ,   cipla%m_taue                                         &
     &          ,   cipla%m_tau2                                         &
     &          ,  (cipla%m_vrot_ci(j),j=1,3)

              call cipla%set_plan_ci_type(ci_type)
             
              ci => cipla
             else
               deallocate(cipla)
               ci =>NULL()
             endif
             
            end function t_ciread_circummaree_freadline

      end module mod_ciread_circummaree
