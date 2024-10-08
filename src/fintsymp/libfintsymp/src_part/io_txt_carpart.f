!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file io_txt_carpart.f 
!!  \brief ecriture en ascii sur fichier de coordonnees cartesiennes sur 1 ligne
!!         pour les particules
!!
!!     cette classe t_io_txt_car_part gere la sortie de n particules ayant 6 composantes
!!  sous  la forme de 7 colonnes par ligne avec le temps au debut.\n
!!  Chaque ligne peut etre prefixee par le radical fourni pour chaque particule.
!!     
!!
! history : creation 29/08/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour l'ecriture des coordonnees des particules dans un fichier
!***********************************************************************
      module mod_io_txt_carpart
       use mod_io_txt_mligncomp
       
!***********************************************************************
!> @class t_io_txt_car_part
!! classe d'ecriture des coordonnees des particules dans un fichier
!!  
!***********************************************************************
      type, extends(t_io_txt_mligncomp) :: t_io_txt_car_part

       contains
          
          procedure :: set_part_nb =>io_txt_car_set_part_nb ! fixe le nombre de planetes et de particules
          
      end type t_io_txt_car_part  
      
       interface t_io_txt_car_part
            module procedure io_txt_car_part_dotname 
       end interface  t_io_txt_car_part

      contains
         
!***********************************************************************
!> @brief constructeur qui specifie le nom pour le graphique DOT
!***********************************************************************
       function io_txt_car_part_dotname (dotname) result (iotxt)
        implicit none
        character(len=*), intent (in) :: dotname !< nom pour le graphique dot
        type ( t_io_txt_car_part ) :: iotxt
        iotxt%m_dotname = dotname
       end function io_txt_car_part_dotname

!***********************************************************************
!> @brief fixe le nombre de planetes et de particules
!! le nombre de planetes est utilise pour ignore cette colonne
!***********************************************************************
            subroutine io_txt_car_set_part_nb(this, nbplan, nbpart)
             implicit none
             class(t_io_txt_car_part), intent(inout):: this  !< dummy argument
             integer, intent(in) :: nbplan  !< nombre de planetes
             integer, intent(in) :: nbpart  !< nombre de particules
             
             integer, dimension(7, nbpart) :: newindex
             integer j,k
             
             do k=0, nbpart-1
              newindex(1,k+1) = 1 ! pour le temps
              do j=1,3
               newindex(j+1,k+1) = 3*k+j+1+6*nbplan
              enddo 
              do j=1,3
               newindex(j+3+1,k+1) = 3*nbpart+3*k+j+1+6*nbplan
              enddo 
             enddo
             write(*,*) 'newindex', newindex
             call  io_txt_mligncomp_reorder(this,newindex)
            end  subroutine io_txt_car_set_part_nb

      end module mod_io_txt_carpart
