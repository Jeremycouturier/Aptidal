!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysparth.f 
!!  \brief systeme avec M particules 
!!         pouvant etre integre par un schema symplectique
!!         les coordonnees des corps sont exprimees en heliocentriques (position/vitesse).
!!         seules les interactions newtonniennes sont prises en compte.
!!
!!
!
! history : creation 26/06/2014
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec M particules exprimee en heliocentrique
!***********************************************************************
      module mod_syspartnewtH
      use mod_syspart
      
!***********************************************************************
!> @class t_syspartnewtH
!! enesemble de  M particules dont 
!!  les coordonnees des corps sont exprimees en heliocentriques (position/vitesse).
!!   seules les interactions newtonniennes sont prises en compte.
!!
!!  
!***********************************************************************
      type, extends(t_syspart) ::  t_syspartnewtH

       contains
          
          procedure :: pasA =>syspartnewtH_pasA
          procedure :: pasB1 =>syspartnewtH_pasB1
          procedure :: pasB2 =>syspartnewtH_pasB2
          procedure :: pasC =>syspartnewtH_pasC
         
      end type t_syspartnewtH  
     
      contains
         

!***********************************************************************
!> @brief execute le pas A de l'integrateur pour les particules
!***********************************************************************
      subroutine syspartnewtH_pasA(this, syspla, tau)
       use kahansum
       use mod_kepsaut
       implicit none
       class(t_syspartnewtH), intent(inout):: this  !< dummy argument
       class(t_syspla), intent(inout) :: syspla !< systeme planetaire
       real(TREAL), intent(in) :: tau !<= tau
       
       real(TREAL) :: mu_helio
       
       real(TREAL) ,dimension(3):: X, Xp, err, errp
       real(TREAL) ,dimension(3):: DX, DXp 
       integer i, npart, npartinteg
       npart = this%m_part_nb
       
       mu_helio = this%m_cG*this%m_star_mass
       npartinteg = 0
      do i= 1,npart
         if (this%m_part_mce(i)%m_stop.eq.0) then
             npartinteg = npartinteg+1
             X = this%m_part_coordX(:,i)
             err = this%m_part_errX(:,i)
             Xp = this%m_part_coordXP(:,i)
             errp = this%m_part_errXP(:,i)
             call keps_new_aut(X,Xp,mu_helio,tau,DX,DXp,                    &
     &              this%m_part_mce(i))

             call  kahansum3(X,err, DX)
             call  kahansum3(Xp,errp, DXp)
             this%m_part_coordX(:,i) = X
             this%m_part_coordXP(:,i) = Xp
             this%m_part_errX(:,i) = err
             this%m_part_errXP(:,i) = errp
         endif
      end do
      
      ! verifie si au moins  une particule integree
      ! si acune particule => forcer l'arret de l'integration
      if (npartinteg.eq.0) then
        call syspla%m_plan_arret%set_error(-8)
        call syspla%m_plan_arret%set_time(real(0.,kind=TREAL))
        call syspla%m_plan_arret%set_value(real(0.,kind=TREAL))
        call syspla%m_plan_arret%set_body(7)
      endif

     
      end  subroutine syspartnewtH_pasA  

!***********************************************************************
!> @brief execute le pas B2 de l'integrateur pour les particules
!***********************************************************************
      subroutine syspartnewtH_pasB2(this, syspla,tau)
       use kahansum
       implicit none
       class(t_syspartnewtH), intent(inout):: this  !< dummy argument
       class(t_syspla), intent(inout) :: syspla !< systeme planetaire
       real(TREAL), intent(in) :: tau !< pas de temps
       
       !write(*,*) "syspartnewtH_pasB2", tau
       
      real(TREAL),dimension(3) :: TD
      real(TREAL),dimension(3) :: Acc
      real(TREAL),dimension(3) :: Vect
      real(TREAL), dimension(3):: errp, Vb
      real(TREAL) :: mplj
      integer i, j, nplan, npart
      
      nplan = syspla%m_plan_nb
      npart = this%m_part_nb
      
      do i=1,npart
         if (this%m_part_mce(i)%m_stop.eq.0) then
             Acc = REALC(0.E0)
             do j=1,nplan
                Vect=this%m_part_coordX(:,i)-syspla%m_plan_coordX(:,j)
                TD = Vect/sqrt(dot_product(Vect,Vect))**3
                mplj = syspla%m_plan_mass(j)
                Acc(:) = Acc(:) + mplj*TD(:)
             end do
             Acc = this%m_cG*Acc
         
    ! -- step with the compensated sommation !! 
    !!      XPpart = XPpart - tau * Acc
             Vb = this%m_part_coordXP(:,i)
             errp = this%m_part_errXP(:,i)
             call  kahansum3(Vb,errp, - tau*Acc)
             this%m_part_coordXP(:,i) = Vb
             this%m_part_errXP(:,i) = errp
         endif 

      end do
       
      end  subroutine syspartnewtH_pasB2  

!***********************************************************************
!> @brief execute le pas B1 de l'integrateur pour les particules
!***********************************************************************
      subroutine syspartnewtH_pasB1(this, syspla, tau)
       use kahansum
       use mod_sysplanewtH
       implicit none
       class(t_syspartnewtH), intent(inout):: this  !< dummy argument
       class(t_syspla), intent(inout) :: syspla !< systeme planetaire
       real(TREAL), intent(in) :: tau !< pas de temps
 
       real(TREAL), dimension(3):: aux
       real(TREAL), dimension(3):: delt,err,Xh
       integer iv, ip, nplan, npart
       
        select type(syspla)
         class is(t_sysplanewtH)
         
!        X(:,ip) = X(:,ip) 
!     &             + tau*(aux(:,ip)-mpsm0(ip)*XP(:,ip))

      nplan = syspla%m_plan_nb
      npart = this%m_part_nb
      do iv = 1,3
         aux(iv)=dot_product(syspla%m_plan_mpsm0,                       &
     &                         syspla%m_plan_coordXP(iv,:))
      enddo 

      delt = tau*aux
      do ip=1,npart
         if (this%m_part_mce(ip)%m_stop.eq.0) then
             Xh = this%m_part_coordX(:,ip)
             err = this%m_part_errX(:,ip)
             call  kahansum3(Xh,err, delt)
             this%m_part_coordX(:,ip) = Xh
             this%m_part_errX(:,ip) = err
         endif 
      end do

         class default
          stop 'syspla in syspartnewtH_pasB1 : bad class'
        end select 

      end  subroutine syspartnewtH_pasB1 


!***********************************************************************
!> @brief execute le pas C de l'integrateur pour les particules
!***********************************************************************
      subroutine syspartnewtH_pasC(this, syspla, tau)
       implicit none
       class(t_syspartnewtH), intent(inout):: this  !< dummy argument
       class(t_syspla), intent(inout) :: syspla !< systeme planetaire
       real(TREAL), intent(in) :: tau !<= tau
       
       write(*,*) "syspartnewtH_pasC", tau
       write(*,*) "pas de correcteur en heliocentrique"
       stop
       
      end  subroutine syspartnewtH_pasC  


      end module mod_syspartnewtH
