!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysplanewth.f 
!!  \brief systeme avec N planetes
!!         pouvant etre integre par un schema symplectique.
!!         les coordonnees des corps sont exprimees en heliocentriques (position/vitesse).
!!         seules les interactions newtonniennes sont prises en compte.
!!
!!
!
! history : creation 20/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec N planetes 
! exprime en variables heliocentriques
! seules les interactions newtonniennes sont prises en compte
!***********************************************************************
      module mod_sysplanewtH
       use mod_syspla
       private mon_cinH
       private Keplerian_EnergyH, Perturb_EnergyH
       
!***********************************************************************
!> @class t_sysplanewtH
!! classe decrivant 1 systeme planetaire avec N planetes
!! exprimee en variables heliocentriques.
!! seules les interactions newtonniennes sont prises en compte.
!!  
!!  
!!  
!***********************************************************************
      type, extends(t_syspla) :: t_sysplanewtH 

          real(TREAL), dimension(:), allocatable ::   m_plan_mu_helio !< mu heliocentrique
          real(TREAL), dimension(:), allocatable ::   m_plan_mpsm0    !< masse planete / masse etoile
          
       contains
          private
          
          procedure, public:: set_plan_mass => sysplanewtH_set_plan_mass ! fixe la masse des planetes 
          
          procedure, public :: pasA =>sysplanewtH_pasA
          procedure, public :: pasB =>sysplanewtH_pasB
          procedure, public :: pasC =>sysplanewtH_pasC
          
          procedure, public :: energie   =>sysplanewtH_energie
          procedure, public :: energie_with_err   =>                    &
     &                         sysplanewtH_energie_with_err
          procedure :: sysplanewtH_keplerian_energy
          procedure ::  sysplanewtH_perturbation_energy
          procedure, public :: mom_cin   =>sysplanewtH_mom_cin 
          procedure, public :: mom_cin_with_err   =>                    &
     &                         sysplanewtH_mom_cin_with_err 

          procedure, public :: dumpbin =>  sysplanewtH_dumpbin ! dump binaire vers un fichier
          procedure, public :: restorebin => sysplanewtH_restorebin ! restauration binaire depuis un fichier

      end type t_sysplanewtH  
      
      contains
         
!***********************************************************************
!> @brief fixe la masse des planetes 
!!  la fonction set_star_mass doit etre appellee auparavant
!***********************************************************************
      subroutine sysplanewtH_set_plan_mass(this, m)
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:), intent(in) :: m  !< masse des planetes
       real(TREAL) :: m0
       
       call syspla_set_plan_mass(this, m)
       call this%set_graphdotname("SYSPLA Helio")

       allocate(this%m_plan_mu_helio(this%m_plan_nb))
       if( any( shape(this%m_plan_mu_helio) /= shape(m) ) ) then
          stop "Erreur : set_plan_mass : dimensions differentes"
       endif
       allocate(this%m_plan_mpsm0(this%m_plan_nb))

       m0 = this%m_star_mass
       this%m_plan_mu_helio = this%m_cG*(m+m0)
       this%m_plan_mpsm0 = m/m0

       
      end  subroutine sysplanewtH_set_plan_mass  

!***********************************************************************
!> @brief execute le pas A de l'integrateur
!***********************************************************************
      subroutine sysplanewtH_pasA(this, mdt)
       use kahansum
       use mod_kepsaut
       use mod_kepsaut_hyper
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= tau
       
       real(TREAL) ,dimension(3):: X, Xp, err, errp, Vb, Xperr
       real(TREAL) ,dimension(3):: DX, DXp 
       integer i, nplan
       real(TREAL) :: tau, mpsbeta, mu_helio
       
       ! write(*,*) "sysplanewtH_pasA", mdt
       tau = mdt
       nplan = this%m_plan_nb
       
      do i= 1,nplan
         X = this%m_plan_coordX(:,i)
         err = this%m_plan_errX(:,i)
         Vb = this%m_plan_coordXP(:,i)
         errp = this%m_plan_errXP(:,i)
         mpsbeta = this%m_plan_mpsbeta(i)
         mu_helio = this%m_plan_mu_helio(i)
         Xp = mpsbeta*Vb
         Xperr = mpsbeta*errp
         
         call keps_new_aut(X,Xp,mu_helio,tau,DX,DXp,this%m_plan_arret)
         !call kepsaut(X,Xp,mu_helio,tau,DX,DXp,this%m_plan_arret)
         
         call  kahansum3(X,err, DX)
         call  kahansum3(Vb,errp, DXp/mpsbeta)

         this%m_plan_coordX(:,i) = X
         this%m_plan_coordXP(:,i) = Vb
         this%m_plan_errX(:,i) = err
         this%m_plan_errXP(:,i) = errp
         
      end do
     
      end  subroutine sysplanewtH_pasA  

!***********************************************************************
!> @brief execute le pas B2 de l'integrateur
!***********************************************************************
      subroutine sysplanewtH_pasB2(this, tau)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: tau !< pas de temps
       
      real(TREAL),dimension(3,this%m_plan_nb,this%m_plan_nb) :: TD
      real(TREAL),dimension(3,this%m_plan_nb)  :: Acc
      real(TREAL),dimension(3)                 :: Vect
      real(TREAL), dimension(3):: errp, Vb
      real(TREAL) :: mplj
      integer i, j, nplan
      
      !write(*,*) "sysplanewtH_pasB2", tau
      nplan = this%m_plan_nb
      
!----- Calcul du tableaux TD 
      do i=1,nplan
         do j=i+1,nplan
            Vect= this%m_plan_coordX(:,i)-this%m_plan_coordX(:,j)
            TD(:,i,j) = Vect/sqrt(dot_product(Vect,Vect))**3
            TD(:,j,i) = -TD(:,i,j)
         end do
      end do
!----- Calcul des acceleration Acc
      Acc = REALC(0.E0)
      do i=1,nplan
         do j=1,i-1
            mplj = this%m_plan_mass(j)
            Acc(:,i) = Acc(:,i) + mplj*TD(:,i,j)
         enddo
         do j=i+1,nplan
            mplj = this%m_plan_mass(j)
            Acc(:,i) = Acc(:,i) + mplj*TD(:,i,j)
         enddo
      enddo 
      Acc = this%m_cG*Acc

! -- step with the compensated somation !! 
!!      XPplan = XPplan - tau * Acc
      do i=1,nplan
         Vb = this%m_plan_coordXP(:,i)
         errp = this%m_plan_errXP(:,i)
         call  kahansum3(Vb,errp, - tau*Acc(:,i))
         this%m_plan_coordXP(:,i) = Vb
         this%m_plan_errXP(:,i) = errp
      end do
       
      end  subroutine sysplanewtH_pasB2  

!***********************************************************************
!> @brief execute le pas B1 de l'integrateur
!***********************************************************************
      subroutine sysplanewtH_pasB1(this, tau)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: tau !< pas de temps
       
       real(TREAL), dimension(3,this%m_plan_nb):: aux
       real(TREAL), dimension(3):: delt,err,Xh
       integer iv, ip, nplan
       
       !write(*,*) "sysplanewtH_pasB1", tau

!        X(:,ip) = X(:,ip) 
!     &             + tau*(aux(:,ip)-mpsm0(ip)*XP(:,ip))

      nplan = this%m_plan_nb
      do iv = 1,3
         aux(iv,:)=dot_product(this%m_plan_mpsm0,                        &
     &                         this%m_plan_coordXP(iv,:))
      enddo 

      do ip=1,nplan
         delt = tau*(aux(:,ip)-this%m_plan_mpsm0(ip)                     &
     &     *this%m_plan_coordXP(:,ip))

         Xh = this%m_plan_coordX(:,ip)
         err = this%m_plan_errX(:,ip)
         call  kahansum3(Xh,err, delt)
         this%m_plan_coordX(:,ip) = Xh
         this%m_plan_errX(:,ip) = err

      end do

      end  subroutine sysplanewtH_pasB1 

!***********************************************************************
!> @brief execute le pas B de l'integrateur = B1(tau)B2(tau)B1(tau) 
!***********************************************************************
      subroutine sysplanewtH_pasB(this, mdt)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
       ! write(*,*) "sysplanewtH_pasB", mdt
       call sysplanewtH_pasB1(this,mdt/REALC(2.E0))
       call sysplanewtH_pasB2(this,mdt)
       call sysplanewtH_pasB1(this,mdt/REALC(2.E0))

      end  subroutine sysplanewtH_pasB  

!***********************************************************************
!> @brief execute le pas C de l'integrateur
!***********************************************************************
      subroutine sysplanewtH_pasC(this, mdt)
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
       write(*,*) "sysplanewtH_pasC", mdt
       write(*,*) "pas de correcteur en heliocentrique"
       stop
       
      end  subroutine sysplanewtH_pasC  


!***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine sysplanewtH_energie(this, H)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), intent(out) :: H  !< energie calculee
       type(compensatedsum) :: Hcomp  !< erreur sur l' energie calculee
       
       call sysplanewtH_energie_with_err(this, Hcomp) 
       H = Hcomp%v

      end  subroutine sysplanewtH_energie  


!***********************************************************************
!> @brief calcule l'energie de this avec l'erreur
!***********************************************************************
      subroutine sysplanewtH_energie_with_err(this, H)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       type(compensatedsum), intent(out) :: H  !< energie calculee
       call EnergieH(this%m_plan_nb,this%m_cG,this%m_star_mass,          &
     &    this%m_plan_mass, this%m_plan_mpsbeta, this%m_plan_coordX,     &
     &    this%m_plan_coordXP,H) 
      ! write(*,*) "sysplanewtH_energie ",H
       
       end  subroutine sysplanewtH_energie_with_err


     
!***********************************************************************
!> @brief calcule l'energie de this : partie keplerienne
!***********************************************************************
      subroutine sysplanewtH_keplerian_energy(this, H)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       type(compensatedsum), intent(out) :: H  !< energie calculee
       call Keplerian_EnergyH(this%m_plan_nb,this%m_cG,this%m_star_mass, &
     &    this%m_plan_mass, this%m_plan_mpsbeta, this%m_plan_coordX,     &
     &    this%m_plan_coordXP,H) 
      ! write(*,*) "sysplanewtH_keplerian_energy ",H
       
       end  subroutine sysplanewtH_keplerian_energy

!***********************************************************************
!> @brief calcule l'energie de this : partie perturbation
!***********************************************************************
      subroutine sysplanewtH_perturbation_energy(this, H)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       type(compensatedsum), intent(out) :: H  !< energie calculee
       call Perturb_EnergyH(this%m_plan_nb,this%m_cG,this%m_star_mass,   &
     &    this%m_plan_mass, this%m_plan_mpsbeta, this%m_plan_coordX,     &
     &    this%m_plan_coordXP,H) 
      ! write(*,*) "sysplanewtH_perturbation_energy ",H
       
      end  subroutine sysplanewtH_perturbation_energy

!***********************************************************************
!> @brief calcule le moment cinetique de this
!***********************************************************************
      subroutine sysplanewtH_mom_cin(this, C)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(1:3), intent(out) :: C  !< moment cinetique
       type(compensatedsum), dimension(1:3) :: Cw

       call sysplanewtH_mom_cin_with_err(this,Cw)
       C(:)=Cw(:)%v
       
      end  subroutine sysplanewtH_mom_cin  

!***********************************************************************
!> @brief calcule le moment cinetique de this
!***********************************************************************
      subroutine sysplanewtH_mom_cin_with_err(this, C)
       use kahansum
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       type(compensatedsum), dimension(1:3), intent(out) :: C  !< moment cinetique
       
       call mon_cinH(this%m_plan_nb,this%m_plan_mass,                    &
     &    this%m_plan_coordX, this%m_plan_coordXP,C) 
       !write(*,*) "sysplanewtH_mom_cin_with_err",C
       
      end  subroutine sysplanewtH_mom_cin_with_err  


!***********************************************************************
!> @brief calcul de l'energie en variables heliocentriques canoniques
!***********************************************************************
      subroutine EnergieH(nplan, cG, m0, mpl, mpsbeta, X,XP,H)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(3,nplan), intent(in) :: X !< position helio
      real(TREAL),dimension(3,nplan), intent(in) :: XP !< vitesse bary
      real(TREAL),dimension(nplan), intent(in) :: mpl !< masse des planetes
      real(TREAL),dimension(nplan), intent(in) :: mpsbeta !< masse/beta des planetes
      real(TREAL), intent(in) :: m0 !< masse du corps central
      real(TREAL), intent(in) :: cG !< constante de gravitation
      type(compensatedsum), intent(out) :: H  !< energie calculee
      type(compensatedsum) :: H0,H1
       
      call Keplerian_EnergyH(nplan, cG, m0, mpl, mpsbeta,               &
     & X, XP, H0)
      call Perturb_EnergyH(nplan,cG, m0, mpl,mpsbeta, X,XP,H1)
      H = H0+H1

      end subroutine EnergieH


!***********************************************************************
!> @brief calcul de l'energie en variables heliocentriques canoniques
!***********************************************************************
      subroutine Keplerian_EnergyH(nplan, cG, m0, mpl, mpsbeta,         &
     &X, XP, H0)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(3,nplan), intent(in) :: X !< position helio
      real(TREAL),dimension(3,nplan), intent(in) :: XP !< vitesse bary
      real(TREAL),dimension(nplan), intent(in) :: mpl !< masse des planetes
      real(TREAL),dimension(nplan), intent(in) :: mpsbeta !< masse/beta des planetes
      real(TREAL), intent(in) :: m0 !< masse du corps central
      real(TREAL), intent(in) :: cG !< constante de gravitation
      type(compensatedsum), intent(out) :: H0 !< energie calculee
      type(compensatedsum) :: T0,U0
      integer :: i
      
      T0 = REALC(0.e0)
      U0 = REALC(0.e0)
      
      do i=1,nplan
          T0=T0+mpl(i)*mpsbeta(i)*dot_product(XP(:,i),XP(:,i))
          U0=U0+mpl(i)/sqrt(dot_product(X(:,i),X(:,i)))
      enddo
       H0 = T0/REALC(2.e0) - cG*m0*U0

      end subroutine Keplerian_EnergyH

!***********************************************************************
!> @brief calcul de l'energie en variables heliocentriques canoniques
!***********************************************************************
      subroutine Perturb_EnergyH(nplan,cG,m0,mpl,mpsbeta,X,XP,H1)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(3,nplan), intent(in) :: X !< position helio
      real(TREAL),dimension(3,nplan), intent(in) :: XP !< vitesse bary
      real(TREAL),dimension(nplan), intent(in) :: mpl !< masse des planetes
      real(TREAL),dimension(nplan), intent(in) :: mpsbeta !< masse/beta des planetes
      real(TREAL), intent(in) :: m0 !< masse du corps central
      real(TREAL), intent(in) :: cG !< constante de gravitation
      type(compensatedsum), intent(out) :: H1  !< energie calculee
      type(compensatedsum) :: T1,U1
      real(TREAL),dimension(3) :: XX
      integer :: i,j
      
      T1 = REALC(0.e0)
      U1 = REALC(0.e0)

      ! write(*,*) 'masses ', mpl
      ! write(*,*) 'mpsbeta ',mpsbeta
      ! write(*,*) X,XP

      do i=1,nplan
         do j=i+1,nplan
          XX = X(:,i) - X(:,j)
           U1 = U1 + mpl(i)*mpl(j)/sqrt(dot_product(XX,XX))
           T1 = T1 + mpl(i)*mpl(j)*dot_product(XP(:,i),XP(:,j))
         end do
      end do
      
        H1 =  - cG*U1 +T1/m0

      end subroutine Perturb_EnergyH 

!***********************************************************************
!> @brief calcul du moment cinetique en variables heliocentriques canoniques
!***********************************************************************
      Subroutine mon_cinH(nplan, mpl,X,XC,C)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(nplan), intent(in) :: mpl !< masse des planetes
      type(compensatedsum),dimension(3), intent(out) :: C !< moment cinetique
      real(TREAL),dimension(3,nplan), intent(in) :: X  !< position helio 
      real(TREAL),dimension(3,nplan), intent(in) :: XC  !< vitesse bary   
      integer :: i
      
      do i=1, 3
            C(i) = REALC(0e0)
      enddo
      do i=1,nplan
         C(1) = C(1) + (X(2,i)*XC(3,i) - X(3,i)*XC(2,i))*mpl(i)
         C(2) = C(2) + (X(3,i)*XC(1,i) - X(1,i)*XC(3,i))*mpl(i)
         C(3) = C(3) + (X(1,i)*XC(2,i) - X(2,i)*XC(1,i))*mpl(i)
      end do         
      end Subroutine mon_cinH

!***********************************************************************
!> @brief dump binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (in) objet a sauvegarder
!***********************************************************************
      subroutine sysplanewtH_dumpbin(this, file) 
       use mod_io_dump
       implicit none
       class(t_sysplanewtH), intent(in):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       call syspla_dumpbin(this, file)

       write(file%m_nf) this%m_plan_mu_helio 
       write(file%m_nf) this%m_plan_mpsm0 
       
      end subroutine  sysplanewtH_dumpbin  

!***********************************************************************
!> @brief restauration binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (inout) objet a sauvegarder
!***********************************************************************
      subroutine sysplanewtH_restorebin(this, file) 
       use mod_io_dump
       implicit none
       class(t_sysplanewtH), intent(inout):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       call syspla_restorebin(this, file)

       if (.not.allocated(this%m_plan_mu_helio)) then
            allocate(this%m_plan_mu_helio(this%m_plan_nb))
            allocate(this%m_plan_mpsm0(this%m_plan_nb))
       endif
       read(file%m_nf) this%m_plan_mu_helio 
       read(file%m_nf) this%m_plan_mpsm0 
 
       end subroutine  sysplanewth_restorebin  

      end module mod_sysplanewtH
