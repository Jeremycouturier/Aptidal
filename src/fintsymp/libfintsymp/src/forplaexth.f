!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file forplaexth.f 
!!  \brief systeme avec N planetes
!!         donnees par un modele "force"
!!         les coordonnees des corps sont exprimees en heliocentriques (position/vitesse).
!!
!
! history : creation 12/09/2016
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec N planetes 
! exprime en variables heliocentriques
! le modele des planetes est un modele "force"
!***********************************************************************
      module mod_forplaextH
       use mod_sysplaextH

!***********************************************************************
!> @class t_forcing
!! classe decrivant le forcage
!! fournit la methode virtuelle "compute" qui retourne les donnees au temps t
!! les donnees retournees par "compute" doivent etre en Ph/Vh 
!!  
!***********************************************************************
      type, abstract ::  t_forcing
       contains
          
          procedure(forcing_compute), deferred :: compute
          
      end type t_forcing  
      
           abstract interface

            subroutine forcing_compute(this, t , X1, X2)
             import t_forcing
             class(t_forcing), intent(inout) :: this
             real(TREAL), intent(in) :: t 
             real(TREAL), dimension(:,:), intent(out) :: X1
             real(TREAL), dimension(:,:), intent(out) :: X2
            end subroutine forcing_compute

          end interface

!***********************************************************************
!> @class t_forplaextH
!! classe decrivant 1 systeme planetaire avec N planetes
!! exprimee en variables heliocentriques.
!! le modele des planetes est un modele "force"
!!  
!!  
!***********************************************************************
      type, extends(t_sysplaextH) ::  t_forplaextH
      
          class(t_forcing), pointer :: m_forcing => NULL() !< methode qui controle le forcage
          
          real(TREAL) :: time !< temps modifie dnas le pas A
          
       contains
          
          procedure, public:: set_forcing => forplaextH_set_forcing !< fixe le modele du forcage

          procedure :: pasA =>forplaextH_pasA
          procedure :: pasB =>forplaextH_pasB
          procedure :: pasC =>forplaextH_pasC
          
          procedure :: getpla => forplaextH_getpla !< recupere la position/vitesse des planetes

          procedure:: copy_output    => forplaextH_copy_output ! copie de this dans la copie de sortie
          
      end type t_forplaextH  
      
      contains
         
!***********************************************************************
!> @brief fixe la methode de forcage et change le nom du graph dot
!***********************************************************************
      subroutine forplaextH_set_forcing(this, f)
       implicit none
       class(t_forplaextH), intent(inout):: this  !< dummy argument
       class(t_forcing), pointer, intent(in) :: f  !< forcage
       
       this%m_forcing => f
       this%time = 0
       call this%set_graphdotname("FOREXT Helio")
      end  subroutine forplaextH_set_forcing  

!***********************************************************************
!> @brief recupere la position/vitesse des planetes au temps t
! modifie les coordonnees des planetes (coordX/XP)
!***********************************************************************
      subroutine forplaextH_getpla(this, t)
       implicit none
       class(t_forplaextH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: t !<= temps 
      real(TREAL),dimension(3,this%m_plan_nb)  :: Vh
       
       call this%m_forcing%compute(t,this%m_plan_coordX,Vh)
       call coord_vh2vb(this%m_plan_nb, this%m_plan_mass0, Vh,             &
     &                  this%m_plan_coordXP)

     
      end  subroutine forplaextH_getpla  

!***********************************************************************
!> @brief execute le pas A de l'integrateur
!***********************************************************************
      subroutine forplaextH_pasA(this, mdt)
       implicit none
       class(t_forplaextH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= tau
       
       integer j
       
        this%time = this%time+mdt
        call this%getpla(this%time)
        do j=1, this%m_nbextension
          call this%m_ext(j)%m_ptr%pasA(this, mdt)
        enddo
     
      end  subroutine forplaextH_pasA  


!***********************************************************************
!> @brief execute le pas B de l'integrateur = B1(tau)B2(tau)B1(tau) 
!***********************************************************************
      subroutine forplaextH_pasB(this, mdt)
       implicit none
       class(t_forplaextH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
       integer j

       !write(*,*) 'pas B t0=',this%time
       call this%getpla(this%time)
       !write(*,*) 't0= X', this%m_plan_coordX
       !write(*,*) 't0= XP', this%m_plan_coordXP
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasB1(this, mdt/REALC(2.E0))
       enddo
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasB2(this, mdt)
       enddo
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasB1(this, mdt/REALC(2.E0))
       enddo

      end  subroutine forplaextH_pasB  

!***********************************************************************
!> @brief execute le pas C de l'integrateur
!***********************************************************************
      subroutine forplaextH_pasC(this, mdt)
       implicit none
       class(t_forplaextH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt

       integer j
       
       call this%getpla(this%time)
       do j=1, this%m_nbextension
         call this%m_ext(j)%m_ptr%pasC(this, mdt)
       enddo
       
      end  subroutine forplaextH_pasC  

!***********************************************************************
!> @brief copie de this dans la copie de sortie (sysout)
!! en general, sysout = this 
!***********************************************************************
      subroutine forplaextH_copy_output(this, sysout)
       implicit none
       class(t_forplaextH), intent(inout):: this  !< dummy argument
       class(t_syssymp), intent(inout) :: sysout  !< systeme cree pour la sortie

       select type(sysout)
        class is(t_forplaextH)
          call sysplaextH_copy_output(this, sysout)
          sysout%time = this%time
              
        class default
         stop 'Erreur: this ne derive pas de forplaextH_copy_output' 
      end select
      end  subroutine forplaextH_copy_output

      end module mod_forplaextH
      