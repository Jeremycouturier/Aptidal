!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file io_txt_car.f 
!!  \brief ecriture en ascii sur fichier de coordonnees cartesiennes sur 1 ligne
!!
!!     cette classe t_io_txt_car gere la sortie de n corps ayant 6 composantes sous  la forme de 6*n+1 colonnes.\n
!!  Chaque ligne peut etre prefixee par un radical fourni.
!!     
!!
! history : creation 07/04/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture des integrales premieres dans un fichier
!***********************************************************************
      module mod_io_txt_car
       use mod_io_txt_ncomp
       
!***********************************************************************
!> @class t_io_txt_car
!! classe d'ecriture des integrales premieres dans un fichier
!!  
!***********************************************************************
      type, extends(t_io_txt_ncomp) :: t_io_txt_car

       contains
          
          procedure :: set_plan_nb =>io_txt_car_set_plan_nb ! fixe le nombre de planetes 
          
      end type t_io_txt_car  
      
       interface t_io_txt_car
            module procedure io_txt_car_dotname 
       end interface  t_io_txt_car

      contains
         
!***********************************************************************
!> @brief constructeur qui specifie le nom pour le graphique DOT
!***********************************************************************
       function io_txt_car_dotname (dotname) result (iotxt)
        implicit none
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        type ( t_io_txt_car ) :: iotxt
        iotxt%m_dotname = dotname
       end function io_txt_car_dotname

!***********************************************************************
!> @brief fixe le nombre de planetes
!***********************************************************************
            subroutine io_txt_car_set_plan_nb(this, nbplan)
             implicit none
             class(t_io_txt_car), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nbplan  !< nombre de planetes
             
             integer, dimension(6*nbplan+1) :: newindex
             integer j,k
             
             newindex(1) = 1 ! pour le temps
             do k=0, nbplan-1
              do j=1,3
               newindex(6*k+j+1) = 3*k+j+1
              enddo 
              do j=1,3
               newindex(6*k+j+3+1) = 3*nbplan+3*k+j+1
              enddo 
             enddo
             call  io_txt_ncomp_reorder(this,newindex)
            end  subroutine io_txt_car_set_plan_nb

      end module mod_io_txt_car
