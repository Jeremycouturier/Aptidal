!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file forplatab.f 
!!  \brief definition du forcage en utilisant une fonction tabulee pour des planetes.
!!
!
! history : creation 15/09/2016
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour le forcage en utilisant une fonction tabulee pour des planetes
!***********************************************************************
      module mod_forplatab
       use mod_forplaextH
       use mod_interdriver

!***********************************************************************
!> @class t_forcingplatab
!! classe de forcage en utilisant une fonction tabulee pour des planetes
!! fournit la methode virtuelle "compute" qui retourne les donnees au temps t
!!  
!***********************************************************************
      type, extends(t_forcing) ::  t_forcingplatab
      
          type (t_interdriver) :: m_driver !< driver d'interpolation
          integer :: m_typecoord !<= type de cooordonnees
          integer :: m_nbpla !<= nombre de planetes
       contains
          procedure :: load => forcingplatab_load !< icharge le fichiers et calcule les coeffiicents d'interpolation
          
          procedure :: compute => forcingplatab_compute !< interpolation des positions/vitesses
          
      end type t_forcingplatab  
      
      
      contains
         
!***********************************************************************
! charge le fichier nomfich pour nbpla qui sont dans un jeu de coordonnees
! le fichier contient donc ncol+1 colonnes dans ce fichier
!***********************************************************************
         subroutine forcingplatab_load(this, typecoord, nomfich, nbpla)
          use modreadfile
          implicit none
          class(t_forcingplatab), intent(inout):: this  !< dummy argument
          integer, intent(in) :: typecoord !< type de coordonnees
          character(len=*) , intent(in) :: nomfich !<= nom du fichier
          integer, intent(in) :: nbpla !<= nombre de planetes
          
          this%m_typecoord = typecoord
          if (this%m_typecoord.ne.5) then
            stop 'typecoord doit etre 5 pour la solution forcee'
          endif

          this%m_nbpla = nbpla
          call this%m_driver%load(nomfich, nbpla*6)
          
         end subroutine  forcingplatab_load

!***********************************************************************
!> @brief interpole la position heliocentrique /vitesse heliocentrique des planetes
! X1  contiendra la position heliocentrique
! X2  conteindra la vitesse heliocentrique
!***********************************************************************
      subroutine forcingplatab_compute(this, t , X1, X2)
       implicit none
       class(t_forcingplatab), intent(inout) :: this
       real(TREAL), intent(in) :: t  !<= temps
       real(TREAL), dimension(:,:), intent(out) :: X1 !<= position
       real(TREAL), dimension(:,:), intent(out) :: X2 !<= vitesse
       real(TREAL), dimension(1:6*this%m_nbpla) :: y
       integer j, offpla

       call this%m_driver%calc(t, y)
       do j=1, this%m_nbpla
           offpla = 6*(j-1)
           X1(:,j) = y(offpla+1:offpla+3)
           X2(:,j) = y(offpla+4:offpla+6)
       enddo

      end  subroutine forcingplatab_compute  

      end module mod_forplatab
      