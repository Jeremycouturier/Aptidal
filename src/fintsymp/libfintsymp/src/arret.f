!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file arret.f 
!!  \brief Implementation of :
!!    gestion des arrets avec un code specifique
!!
!!    ce fichier ne doit rien contenir de specifique a un projet 
!
!
!
! history : creation 2014/03/19
!***********************************************************************
#include "realc.f"

!***********************************************************************
!> module pour la gestion des arrets
!***********************************************************************
      module mod_arret
      
!***********************************************************************
!> @class t_arret
!! type stockant les causes de l'erreur
!! les valeurs negatives de m_cause sont reservees a libasdftools et libfintsymp
!! les valeurs positives sont reservees pour les applications
!!
!! liste des codes d'erreur (champ m_cause) : 
!!  * -3 : probleme de convergence dans kepsaut
!!  * -4 : cas non elliptique
!!  * -5 : trop grande variation de l'energie 
!!  * -6 : distance trop proche de l'etoile
!!  * -7 : distance trop loin de l'etoile
!!  * -8 : plus aucune particule (ou autre) a integrer
!!  * -9 : distance trop proche d'une planete
!***********************************************************************
      type :: t_arret
       integer, public :: m_stop  = 0 !< = 0 => pas d'erreur, = 1 => erreur presente (il faut tester sur le 0)
       integer, public :: m_cause = 0 !< enumeration de la raison de l'erreur
       real(TREAL), public  :: m_time  !< temps ou s'est produite l'erreur
       logical, public  :: m_time_initialized = .false.  !< = true si le champ m_time est valide
       real(TREAL), public  :: m_value = 0 !< valeur associee a l'erreur
       integer, public :: m_body = 0 !< numero du corps generant l'erreur 
       
      contains
       procedure, NON_OVERRIDABLE :: set_error => arret_set_error !< indique qu'une erreur s'est produite
       procedure, NON_OVERRIDABLE :: set_time => arret_set_time !< indique que le temps auquel l'erreur s'est produite
       procedure, NON_OVERRIDABLE :: set_body => arret_set_body !< indique quel corps a genere l'erreur
       procedure, NON_OVERRIDABLE :: set_value => arret_set_value !< indique que le temps auquel l'erreur s'est produite
       
      end type t_arret
      
      contains
      
!***********************************************************************
!> @brief indique qu'une erreur s'est produite en mettant a jour son code
!***********************************************************************
      subroutine arret_set_error(this, code)
       implicit none
       class(t_arret), intent(inout):: this    !<  causes des erreurs
       integer, intent(in):: code    !<  code de l'erreur produite
       
       this%m_stop = 1
       this%m_cause = code
       
      end subroutine arret_set_error  
         
!***********************************************************************
!> @brief indique qu'une erreur s'est produite en mettant a jour son temps
!***********************************************************************
      subroutine arret_set_time(this, temps)
       implicit none
       class(t_arret), intent(inout):: this    !<  causes des erreur
       real(TREAL), intent(in):: temps    !< temps a laquelle s'est produite l'erreur
       
       this%m_stop = 1
       this%m_time = temps
       this%m_time_initialized = .true.
       
      end subroutine arret_set_time  

!***********************************************************************
!> @brief indique qu'une erreur s'est produite en mettant a jour sa valeur
!***********************************************************************
      subroutine arret_set_value(this, value)
       implicit none
       class(t_arret), intent(inout):: this    !<  causes des erreur
       real(TREAL), intent(in):: value    !< valeur associee a l'erreur
       
       this%m_stop = 1
       this%m_value = value
       
      end subroutine arret_set_value  

!***********************************************************************
!> @brief indique quel corps a genere l'erreur
!***********************************************************************
      subroutine arret_set_body(this, corps)
       implicit none
       class(t_arret), intent(inout):: this    !<  causes des erreurs
       integer, intent(in):: corps    !<  numero du corps
       
       this%m_stop = 1
       this%m_body = corps
       
      end subroutine arret_set_body  

!***********************************************************************
!> @brief met a jour le temps des erreurs de mce qui ont un stop!=0 
!! et qui ont leur temps non encore initialise
!***********************************************************************
      subroutine array_arret_update_time(mce, t)
       implicit none
       type(t_arret), dimension(:), intent(inout):: mce  !< tableau d'arrets dont on veut mettre a jour les temps
       real(kind=TREAL), intent(in) :: t !< time 
       
       integer k

          do k=1, size(mce)
            if ((mce(k)%m_stop.ne.0).and.                                 &
     &          (.not.mce(k)%m_time_initialized)) then
                call mce(k)%set_time(t)
            endif
          enddo
           
      end  subroutine array_arret_update_time  


      end module mod_arret
