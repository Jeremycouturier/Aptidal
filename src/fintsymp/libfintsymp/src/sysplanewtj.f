!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file sysplanewtj.f 
!!  \brief systeme avec N planetes
!!         pouvant etre integre par un schema symplectique.
!!         les coordonnees des corps sont exprimees en jacobi (position/vitesse).
!!         seules les interactions newtonniennes sont prises en compte.
!!
!!
!
! history : creation 20/07/2012
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme avec N planetes 
! exprime en variables en jacobi
! seules les interactions newtonniennes sont prises en compte
!***********************************************************************
      module mod_sysplanewtJ
       use mod_syspla
       use kahansum
       private EnergieJ,  mon_cinJ
       private Keplerian_EnergyJ, Perturbation_EnergyJ
       
!***********************************************************************
!> @class t_sysplanewtJ
!! classe decrivant 1 systeme planetaires avec N planetes
!! exprimee en variables heliocentriques.
!! seules les interactions newtonniennes sont prises en compte.
!!  
!!  
!!  
!***********************************************************************
      type, extends(t_syspla) :: t_sysplanewtJ 

          real(TREAL), dimension(:), allocatable ::   m_plan_mu_helio !< mu heliocentrique
          real(TREAL), dimension(:), allocatable ::   m_plan_mpsm0    !< masse planete / masse etoile
          real(TREAL), dimension(:), allocatable ::   m_plan_eta    !< eta jacobi
          real(TREAL), dimension(:), allocatable ::   m_plan_mu     !< mu jacobi
          
       contains
          
          procedure :: set_plan_mass => sysplanewtJ_set_plan_mass ! fixe la masse des planetes 

          procedure :: pasA =>sysplanewtJ_pasA
          procedure :: pasB =>sysplanewtJ_pasB
          procedure :: pasC =>sysplanewtJ_pasC
          
          procedure :: energie   =>sysplanewtJ_energie
          procedure :: energie_with_err   =>                            &
     &                         sysplanewtJ_energie_with_err
          procedure :: sysplanewtJ_keplerian_energy
          procedure :: sysplanewtJ_perturbation_energy
          procedure :: mom_cin   =>sysplanewtJ_mom_cin 
          procedure, public :: mom_cin_with_err =>                      &
     &                         sysplanewtJ_mom_cin_with_err 

          procedure :: create_output  => sysplanewtJ_create_output! creation d'une copie pour effectuer les pas de sortie

          procedure :: dumpbin =>  sysplanewtJ_dumpbin ! dump binaire vers un fichier
          procedure :: restorebin => sysplanewtJ_restorebin ! restauration binaire depuis un fichier

      end type t_sysplanewtJ  

      contains
         
!***********************************************************************
!> @brief fixe la masse des planetes 
!!  la fonction set_star_mass doit etre appellee auparavant
!***********************************************************************
      subroutine sysplanewtJ_set_plan_mass(this, m)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(:), intent(in) :: m  !< masse des planetes
       real(TREAL) :: m0
       integer :: i
       
       call syspla_set_plan_mass(this, m)
       call this%set_graphdotname("SYSPLA Jacobi")
       
       allocate(this%m_plan_mu_helio(this%m_plan_nb))
       if( any( shape(this%m_plan_mu_helio) /= shape(m) ) ) then
          stop "Erreur : set_plan_mass : dimensions differentes"
       endif
       allocate(this%m_plan_mpsm0(this%m_plan_nb))
       allocate(this%m_plan_mu(this%m_plan_nb))
       allocate(this%m_plan_eta(0:this%m_plan_nb))

       m0 = this%m_star_mass
       this%m_plan_mu_helio = this%m_cG*(m+m0)
       this%m_plan_mpsm0 = m/m0
       
       this%m_plan_eta(0) = m0 
       do i=1,this%m_plan_nb
         this%m_plan_eta(i) = this%m_plan_eta(i-1) + m(i)
       end do 
       this%m_plan_mu= this%m_cG*this%m_plan_eta(1:)
       write(*,*)'this%m_plan_eta' , this%m_plan_eta
       write(*,*)'this%m_plan_mu' , this%m_plan_mu

      end  subroutine sysplanewtJ_set_plan_mass  

!***********************************************************************
!> @brief creation d'une copie pour effectuer les pas de sortie
!! en general, duplique simplement this
!***********************************************************************
      subroutine sysplanewtJ_create_output(this, sysout)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       class(t_syssymp), allocatable, intent(out) :: sysout  !< systeme cree pour la sortie
       
       call syspla_create_output(this, sysout)
             
      end  subroutine sysplanewtJ_create_output  

!***********************************************************************
!> @brief execute le pas A de l'integrateur
!***********************************************************************
      subroutine sysplanewtJ_pasA(this, mdt)
       use kahansum
       use mod_kepsaut
       use mod_kepsaut_hyper
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= tau
       
       real(TREAL) ,dimension(3):: X, Xp, err, errp
       real(TREAL) ,dimension(3):: DX, DXp 
       integer i, nplan
       real(TREAL) :: tau, mu

       ! write(*,*) "sysplanewtJ_pasA", mdt
       tau = mdt
       nplan = this%m_plan_nb
       
      do i= 1,nplan
         X = this%m_plan_coordX(:,i)
         err = this%m_plan_errX(:,i)
         Xp = this%m_plan_coordXP(:,i)
         errp = this%m_plan_errXP(:,i)
         mu = this%m_plan_mu(i)
!New kepsaut
         ! call kepsaut(X,Xp,mu,tau, DX, DXp, this%m_plan_arret)
!old kepsaut
         call keps_new_aut(X,Xp,mu,tau, DX, DXp, this%m_plan_arret)

         call  kahansum3(X,err, DX)
         call  kahansum3(Xp,errp, DXp)
         this%m_plan_coordX(:,i) = X
         this%m_plan_coordXP(:,i) = Xp
         this%m_plan_errX(:,i) = err
         this%m_plan_errXP(:,i) = errp

      end do
      
      end  subroutine sysplanewtJ_pasA  

!*************************************************************************
!  Computing the acceleration of the perturbation and the corrector !! 
!  Accelera  and Correction  (note correction computes at the same time 
!  accelera)
!*************************************************************************
      subroutine sysplanewtJ_Accelera(this, Acc)
      implicit none
      class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
      real(TREAL),dimension(3,this%m_plan_nb), intent(out)  :: Acc !< acceleration calculee
      
      integer i,j,alp,am1,ap1, nplan
      real(TREAL):: aux, mpli, mplj
      real(TREAL),dimension(3):: Vect,A1,A2,A3,A4
      real(TREAL),dimension(3,this%m_plan_nb):: Ph,Pj,Tv
      real(TREAL),dimension(3,this%m_plan_nb,this%m_plan_nb):: TD
      
 
      nplan = this%m_plan_nb

      Pj = this%m_plan_coordX
      call coord_j2h(nplan,this%m_plan_mass0,this%m_plan_eta,Pj,Ph)
*----- Calcul des tableaux TD et Tv
      do i=1,nplan
         aux = dot_product(Ph(:,i),Ph(:,i))
         aux = sqrt(aux)*aux        
         TD(:,i,i) = Ph(:,i)/aux
         aux = dot_product(Pj(:,i),Pj(:,i))        
         aux = sqrt(aux)*aux        
         Tv(:,i) = Pj(:,i)/aux
         do j=i+1,nplan
            Vect(:) = Ph(:,i)-Ph(:,j)
            aux = dot_product(Vect(:),Vect(:))
            aux = sqrt(aux)*aux        
            TD(:,i,j) = Vect(:)/aux
         end do
      end do

!----- Calcul des acceleration Acc
      Acc = REALC(0.e0)
      do i=2,nplan
         mpli = this%m_plan_mass(i)
         Acc(:,1) = Acc(:,1) + mpli*(TD(:,i,i)+TD(:,1,i))
      end do
 
      do alp = 2, nplan
         am1 = alp-1
         ap1 = alp+1

         Acc(:,alp) = this%m_star_mass*TD(:,alp,alp)
         do i=1,am1
            mpli = this%m_plan_mass(i)
            Acc(:,alp) = Acc(:,alp) - mpli*TD(:,i,alp)
         end do
         Acc(:,alp) = this%m_plan_eta(alp)/this%m_plan_eta(am1)           &
     &                 *Acc(:,alp) - this%m_plan_eta(alp)*Tv(:,alp)

         A1 = REALC(0.e0) 
         A2 = REALC(0.e0) 
         A3 = REALC(0.e0) 
         do i=ap1,nplan
            mpli = this%m_plan_mass(i)
            A1 = A1 + mpli*TD(:,i,i)
            A2 = A2 + mpli*TD(:,alp,i)
            A4 = REALC(0.e0)
            do j=1,am1
               mplj = this%m_plan_mass(j)
               A4 = A4 - mplj*TD(:,j,i)
            end do
            A3 = A3 + mpli*A4
         end do
         A1 = this%m_star_mass/this%m_plan_eta(am1)*A1
         A3 = A3/this%m_plan_eta(am1)
         Acc(:,alp) = Acc(:,alp) + (A1 + A2 + A3) 
      end do
      end subroutine sysplanewtJ_Accelera

!***********************************************************************
!> @brief execute le pas B de l'integrateur
!***********************************************************************
      subroutine sysplanewtJ_pasB(this, mdt)
       use kahansum
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !< pas de temps
       
      real(TREAL),dimension(3,this%m_plan_nb)  :: Acc
      real(TREAL), dimension(3):: errp, Vj
      real(TREAL) :: tau, cG
      integer i, nplan
     
      nplan = this%m_plan_nb
      tau = mdt
      
        
!----- Calcul des acceleration Acc
      call sysplanewtJ_Accelera(this, Acc)
! -- step with the compensated somation !! 
!!      XPplan = XPplan - tau * Acc
      cG = this%m_cG
      do i=1,nplan
         Vj = this%m_plan_coordXP(:,i)
         errp = this%m_plan_errXP(:,i)
         call  kahansum3(Vj,errp, - tau*cG*Acc(:,i))
         this%m_plan_coordXP(:,i) = Vj
         this%m_plan_errXP(:,i) = errp
      end do
      
      end  subroutine sysplanewtJ_pasB 


!***********************************************************************
!> @brief execute le pas C de l'integrateur
!***********************************************************************
      subroutine sysplanewtJ_pasC(this, mdt)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
       write(*,*) "sysplanewtJ_pasC", mdt
       write(*,*) "pas de correcteur en jacobi"
       stop
       
      end  subroutine sysplanewtJ_pasC  


!***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine sysplanewtJ_energie(this, H)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       real(TREAL), intent(out) :: H !< energie calculee
       type(compensatedsum) :: Hcomp
       
       
       call sysplanewtJ_energie_with_err(this, Hcomp)
       H = Hcomp%v
       
      end  subroutine sysplanewtJ_energie  

!***********************************************************************
!> @brief calcule l'energie de this avec l'erreur sur l'energie
!***********************************************************************
      subroutine sysplanewtJ_energie_with_err(this, H)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       type(compensatedsum), intent(out) :: H !< energie calculee
       
       call EnergieJ(this%m_plan_nb,this%m_cG,this%m_star_mass,         &
     &    this%m_plan_mass0, this%m_plan_eta,                           &
     &      this%m_plan_coordX, this%m_plan_coordXP,H) 
       ! write(*,*) "sysplanewtJ_energie_with_err ",H
       
      end  subroutine sysplanewtJ_energie_with_err  

!***********************************************************************
!> @brief calcule l'energie Keplerienne de this
!***********************************************************************
      subroutine sysplanewtJ_keplerian_energy(this, H)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       type(compensatedsum), intent(out) :: H  !< energie calculee
       
       call Keplerian_EnergyJ(this%m_plan_nb,this%m_cG,          
     &    this%m_plan_mass0, this%m_plan_eta,                            &
     &    this%m_plan_coordX, this%m_plan_coordXP,H) 
      ! write(*,*) "sysplanewtJ_keplerian_energy ",H
       
      end  subroutine sysplanewtJ_keplerian_energy


 !***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine sysplanewtJ_perturbation_energy(this, H)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       type(compensatedsum), intent(out) :: H  !< energie calculee
       
       call Perturbation_EnergyJ(this%m_plan_nb,this%m_cG,               &
     &    this%m_star_mass,                                              &
     &    this%m_plan_mass0, this%m_plan_eta,                            &
     &    this%m_plan_coordX, H) 
      ! write(*,*) "sysplanewtJ_perturbation_energy ",H
       
      end  subroutine  sysplanewtJ_perturbation_energy 

!***********************************************************************
!> @brief calcule le moment cinetique de this
!***********************************************************************
      subroutine sysplanewtJ_mom_cin_with_err(this, C)
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       type(compensatedsum), dimension(1:3), intent(out) :: C  !< moment cinetique
       
       integer nplan
       real(TREAL),dimension(3,this%m_plan_nb) :: Ph, Vh, Vb
       nplan = this%m_plan_nb
       
       call coord_j2h(nplan,this%m_plan_mass0,this%m_plan_eta,            &
     &   this%m_plan_coordX,Ph)
       call coord_j2h(nplan,this%m_plan_mass0,this%m_plan_eta,            &
     &   this%m_plan_coordXP,Vh)
       call coord_h2b(nplan, this%m_plan_mass0,Vh,Vb)
       call mon_cinJ(nplan, this%m_plan_mass, Ph,Vb,C) 
       !write(*,*) "sysplanewtJ_mom_cin",C
       
      end  subroutine sysplanewtJ_mom_cin_with_err  


!***********************************************************************
!> @brief calcule le moment cinetique de this
!***********************************************************************
      subroutine sysplanewtJ_mom_cin(this, C)
       use kahansum
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(1:3), intent(out) :: C  !< moment cinetique
       type(compensatedsum), dimension(1:3) :: Cw

       call sysplanewtJ_mom_cin_with_err(this,Cw)
       C(:)=Cw(:)%v
       
      end  subroutine sysplanewtJ_mom_cin  

!***********************************************************************
! calcul de l'energie en variables jacobi 
!***********************************************************************
      subroutine EnergieJ(nplan, cG, m0, mpl, eta, X,XP,H)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(3,nplan), intent(in) :: X !< position jacobi
      real(TREAL),dimension(3,nplan), intent(in) :: XP !< vitesse jacobi
      real(TREAL),dimension(0:nplan), intent(in) :: mpl !< masse des planetes
      real(TREAL),dimension(0:nplan), intent(in) :: eta !< eta des planetes
      real(TREAL), intent(in) :: m0 !< masse du corps central
      real(TREAL), intent(in) :: cG !< constante de gravitation
      type(compensatedsum), intent(out) :: H  !< energie calculee
      type(compensatedsum) :: H0, H1
      
      call Keplerian_EnergyJ(nplan, cG, mpl, eta, X,XP,H0)
      call Perturbation_EnergyJ(nplan, cG, m0, mpl, eta, X,H1)

      H = H0+H1
      
      end subroutine EnergieJ


!***********************************************************************
! calcul de l'energie en variables jacobi : partie keplerienne
!***********************************************************************
      subroutine Keplerian_EnergyJ(nplan, cG, mpl, eta, X,XP,H0)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(3,nplan), intent(in) :: X !< position jacobi
      real(TREAL),dimension(3,nplan), intent(in) :: XP !< vitesse jacobi
      real(TREAL),dimension(0:nplan), intent(in) :: mpl !< masse des planetes
      real(TREAL),dimension(0:nplan), intent(in) :: eta !< eta des planetes
      real(TREAL), intent(in) :: cG !< constante de gravitation
      type(compensatedsum), intent(out) :: H0  !< energie keplerienne calculee
      type(compensatedsum) :: T0,U0
      real(TREAL) ::  delta,aux
     
      integer :: i
      ! integer , save:: nbcall=0
      
      T0 = REALC(0.e0)
      U0 = REALC(0.e0)
      
      do i=1,nplan
         aux = dot_product(XP(:,i),XP(:,i))
         T0 = T0 + aux*mpl(i)*eta(i-1)/eta(i)
         delta = sqrt(dot_product(X(:,i),X(:,i)))
         U0 = U0 + mpl(i)*eta(i-1)/delta
      end do
      
      H0 = REALC(.5e0)*T0 - cG*U0
      
      end subroutine Keplerian_EnergyJ

!***********************************************************************
! calcul de l'energie en variables jacobi  : partie perturbation
!***********************************************************************
      subroutine Perturbation_EnergyJ(nplan, cG, m0, mpl, eta, X,H1)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(3,nplan), intent(in) :: X !< position jacobi
      real(TREAL),dimension(0:nplan), intent(in) :: mpl !< masse des planetes
      real(TREAL),dimension(0:nplan), intent(in) :: eta !< eta des planetes
      real(TREAL), intent(in) :: m0 !< masse du corps central
      real(TREAL), intent(in) :: cG !< constante de gravitation
      type(compensatedsum), intent(out) :: H1  !< energie calculee
      real(TREAL) :: delta,aux
      type(compensatedsum) :: U1
      real(TREAL),dimension(3) :: Xvect
      real(TREAL),dimension(3,nplan) :: Ph 
      integer :: i,j
      ! integer , save:: nbcall=0
      
      U1 = REALC(0.e0)

      
      call coord_j2h(nplan,mpl,eta,X,Ph)
      do i=1,nplan
         if (i.gt.1) then
!--------Attention : ici il y a une simplification
          aux = eta(i-1)/sqrt(dot_product(X(:,i),X(:,i)))
          aux = aux - m0/sqrt(dot_product(Ph(:,i),Ph(:,i)))
          U1 = U1 + mpl(i)*aux   
         end if
         
         do j=i+1,nplan
          Xvect = Ph(:,i) - Ph(:,j)
          delta = sqrt(dot_product(Xvect,Xvect))
          U1 = U1 - mpl(i)*mpl(j)/delta
         end do
      end do
 
      H1 =   cG*U1
            
      end subroutine Perturbation_EnergyJ
      
!***********************************************************************
! calcul du moment cinetique en variables jacobi
!***********************************************************************
      Subroutine mon_cinJ(nplan, mpl,X,XC,C)
      use kahansum
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes
      real(TREAL),dimension(nplan), intent(in) :: mpl !< masse des planetes
      type(compensatedsum),dimension(3), intent(out) :: C !< moment cinetique
      real(TREAL),dimension(3,nplan), intent(in) :: X  !< position jacobi 
      real(TREAL),dimension(3,nplan), intent(in) :: XC  !< vitesse jacobi   
      integer :: i
      
      do i=1,3
      C(i) = REALC(0e0)
      enddo
      do i=1,nplan
         C(1) = C(1) + (X(2,i)*XC(3,i) - X(3,i)*XC(2,i))*mpl(i)
         C(2) = C(2) + (X(3,i)*XC(1,i) - X(1,i)*XC(3,i))*mpl(i)
         C(3) = C(3) + (X(1,i)*XC(2,i) - X(2,i)*XC(1,i))*mpl(i)
      end do         
      end Subroutine mon_cinJ

!***********************************************************************
!> @brief dump binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (in) objet a sauvegarder
!***********************************************************************
      subroutine sysplanewtJ_dumpbin(this, file) 
       use mod_io_dump
       implicit none
       class(t_sysplanewtJ), intent(in):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file

       call syspla_dumpbin(this, file)

       write(file%m_nf) this%m_plan_mu_helio 
       write(file%m_nf) this%m_plan_mpsm0 
       write(file%m_nf) this%m_plan_eta
       write(file%m_nf) this%m_plan_mu
       
      end subroutine  sysplanewtJ_dumpbin  

!***********************************************************************
!> @brief restauration binaire de this vers le fichier num_file
!! file (in) fichier de dump
!! this (inout) objet a sauvegarder
!***********************************************************************
      subroutine sysplanewtJ_restorebin(this, file) 
       use mod_io_dump
       implicit none
       class(t_sysplanewtJ), intent(inout):: this  !< dummy argument
       type(t_io_dump), intent(in) :: file
       
       call syspla_restorebin(this, file)

       if (.not.allocated(this%m_plan_mu_helio)) then
            allocate(this%m_plan_mu_helio(this%m_plan_nb))
            allocate(this%m_plan_mpsm0(this%m_plan_nb))
            allocate(this%m_plan_mu(this%m_plan_nb))
            allocate(this%m_plan_eta(0:this%m_plan_nb))
       endif

       read(file%m_nf) this%m_plan_mu_helio 
       read(file%m_nf) this%m_plan_mpsm0 
       read(file%m_nf) this%m_plan_eta
       read(file%m_nf) this%m_plan_mu
 
       end subroutine  sysplanewtJ_restorebin  

      end module mod_sysplanewtJ
