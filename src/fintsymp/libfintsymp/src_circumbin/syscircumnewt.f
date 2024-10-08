!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file syscircumnewt.f 
!!  \brief systeme planetaires avec 2 etoiles et N planetes
!!         pouvant etre integre par un schema symplectique.
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!!         seules les interactions newtonniennes sont prises en compte.
!!
!!
!
! history : creation 19/11/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
! module pour  1 systeme a 2 etoiles avec N planetes 
! exprime en variables heliocentriques
! seules les interactions newtonniennes sont prises en compte
!***********************************************************************
      module mod_syscircumnewt
       use mod_syspla
       !private mom_cin_circum , Energie_circum
       
!***********************************************************************
!> @class t_syscircumnewt
!! classe decrivant 1 systeme planetaire avec 2 etoiles et N planetes
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!! seules les interactions newtonniennes sont prises en compte.
!!  
!! attention : t_syspla::m_star_mass = masse de la 1ere etoile
!!             t_syspla::m_plan_nb = nombre de planetes +1 (2nde etoile a l'indice 1)
!!             t_syspla::m_plan_mass(1) = masse de la 2nde etoile
!!             t_syspla::m_plan_mass(2:) = masse des planetes
!!             t_syspla::m_plan_mpsbeta(1) = m_etoile1/beta_1 avec beta_1=(m_etoile0*m_etoile1)/(m_etoile0+m_etoile1) pour la 2nde etoile
!!             t_syspla::m_plan_mpsbeta(2:) = m_p/beta_p avec beta_p=((m_etoile+m_etoile1)*m_p)/((m_etoile0+m_etoile1)+m_p) pour les planetes 
!!  
!***********************************************************************
      type, extends(t_syspla) :: t_syscircumnewt 
      
          real(TREAL), dimension(:), allocatable ::   m_mu !< mu(1) = G(somme masses 2 etoiles) , mu(2:) = G*(masses 2 etoiles + masse de la planete)
          real(TREAL) ::m_summass_star !< somme de la masse des 2 etoiles
          
       contains
          
          procedure :: set_plan_nb => syscircumnewt_set_plan_nb ! fixe le nombre d'etoiles et de planetes 
          procedure :: set_mass => syscircumnewt_set_mass ! fixe la masse des etoiles et des planetes
          
          procedure :: pasA =>syscircumnewt_pasA
          procedure :: pasB =>syscircumnewt_pasB
          procedure :: pasC =>syscircumnewt_pasC
          
          procedure :: energie   =>syscircumnewt_energie 
          procedure :: mom_cin   =>syscircumnewt_mom_cin 

      end type t_syscircumnewt  
      
      contains
         
!***********************************************************************
!> @brief fixe le nombre de planetes
!***********************************************************************
      subroutine syscircumnewt_set_plan_nb(this, nplan)
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       integer,  intent(in) :: nplan  !<nombre de planete
       
       call set_plan_nb(this, nplan+1) ! +1 utile pour stocker la 2nde etoile a l'indice 1
      end  subroutine syscircumnewt_set_plan_nb  

!***********************************************************************
!> @brief fixe la masse des 2 etoiles et des planetes
!***********************************************************************
      subroutine syscircumnewt_set_mass(this, cG, m)
       use mod_coordcircumbin
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: cG          !< constante de Gauss AU**3/an**2
       real(TREAL), dimension(0:), intent(in) :: m  !< masse des etoiles et planetes (0:1) = etoiles, (2:) : planetes
       real(TREAL) :: mstar ! somme des masses des etoiles
       integer nplan
       
       nplan = this%get_plan_nb()
       call set_star_mass(this, cG, m(0))
       mstar = m(0)+m(1) 
       this%m_summass_star = mstar
       call this%set_plan_mass(m(1:nplan))
       call this%set_graphdotname("SYS CIRCUMBINAIRE")

       allocate(this%m_mu(1:nplan))
       call circum_calc_mu_mpsbeta(nplan,cG, m(0:nplan),                        &
     &           this%m_mu,this%m_plan_mpsbeta)

      end  subroutine syscircumnewt_set_mass  

!***********************************************************************
!> @brief execute le pas A de l'integrateur
!***********************************************************************
      subroutine syscircumnewt_pasA(this, mdt)
       use kahansum
       use mod_kepsaut
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= tau
       
       real(TREAL) ,dimension(3):: X, Xp, err, errp, Vhat
       real(TREAL) ,dimension(3):: DX, DXp 
       integer i, nplan
       real(TREAL) :: tau, mpsbeta, mu_helio
       
       tau = mdt
       nplan = this%m_plan_nb
       
      do i= 1,nplan
         X = this%m_plan_coordX(:,i)
         err = this%m_plan_errX(:,i)
         Vhat = this%m_plan_coordXP(:,i)
         errp = this%m_plan_errXP(:,i)
         mpsbeta = this%m_plan_mpsbeta(i)
         mu_helio = this%m_mu(i)
         Xp = mpsbeta*Vhat
         
         call keps_new_aut(X,Xp,mu_helio,tau,DX,DXp,this%m_plan_arret)

         call  kahansum3(X,err, DX)
         call  kahansum3(Vhat,errp, DXp/mpsbeta)
         this%m_plan_coordX(:,i) = X
         this%m_plan_coordXP(:,i) = Vhat
         this%m_plan_errX(:,i) = err
         this%m_plan_errXP(:,i) = errp
      end do
     
      end  subroutine syscircumnewt_pasA  

!***********************************************************************
!> @brief execute le pas B2 de l'integrateur
!***********************************************************************
      subroutine syscircumnewt_pasB2(this, tau)
       use kahansum
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: tau !< pas de temps
       
      real(TREAL),dimension(3,this%m_plan_nb,this%m_plan_nb) :: TD
      real(TREAL),dimension(3,this%m_plan_nb)  :: Acc
      real(TREAL),dimension(3)                 :: Vect
      real(TREAL), dimension(3):: errp, Vhat
      real(TREAL) :: mstar ! somme des masses des 2 etoiles
      real(TREAL) :: mplj, factmass, n0, n1, m0,m1, vi3
      real(TREAL),dimension(2:this%m_plan_nb) :: vi_n0_v1,vi_n1_v1
      real(TREAL),dimension(3) :: v1, vi
      integer i, j, nplan

      mstar = this%m_summass_star
      m0 = this%m_star_mass
      m1 = this%m_plan_mass(1)
      n0 = -m1/mstar !=-m1/(m0+m1)
      n1 = m0/mstar !=m0/(m0+m1)
      nplan = this%m_plan_nb
      v1(:) = this%m_plan_coordX(:,1) ! 2nde etoile

      Acc = REALC(0.E0)
      
      do i=2,nplan
        vi = this%m_plan_coordX(:,i)
        vi_n0_v1(i)=REALC(1.E0)/sqrt(dot_product(vi-n0*v1,vi-n0*v1))**3 
        vi_n1_v1(i)=REALC(1.E0)/sqrt(dot_product(vi-n1*v1,vi-n1*v1))**3 
      enddo

!
!----- pour la 2nde etoile----------------
!      
      do i=2,nplan
        factmass = m0*this%m_plan_mass(i)/mstar
        vi = this%m_plan_coordX(:,i)
        Acc(:,1) = Acc(:,1) + factmass*(vi_n0_v1(i)-vi_n1_v1(i))*vi            &
     &         +factmass*(-n0*vi_n0_v1(i)+n1*vi_n1_v1(i))*v1
      enddo
!
!----- pour les planetes----------------
!
      factmass = m0*m1/mstar
      do i=2,nplan
        vi = this%m_plan_coordX(:,i)
        vi3 = REALC(1.E0)/sqrt(dot_product(vi,vi))**3 
        Acc(:,i) = Acc(:,i)                                                   &
     &             + (m0*vi_n0_v1(i)+m1*vi_n1_v1(i)-mstar*vi3)*vi             &
     &             +factmass*(vi_n0_v1(i)-vi_n1_v1(i))*v1
      enddo


!      
!----- Calcul du tableaux TD 
      do i=2,nplan
         do j=i+1,nplan
            Vect= this%m_plan_coordX(:,i)-this%m_plan_coordX(:,j)
            TD(:,i,j) = Vect/sqrt(dot_product(Vect,Vect))**3
            TD(:,j,i) = -TD(:,i,j)
         end do
      end do
!----- Calcul des acceleration Acc
      do i=2,nplan
         do j=2,i-1
            mplj = this%m_plan_mass(j)
            Acc(:,i) = Acc(:,i) + mplj*TD(:,i,j)
         enddo
         do j=i+1,nplan
            mplj = this%m_plan_mass(j)
            Acc(:,i) = Acc(:,i) + mplj*TD(:,i,j)
         enddo
      enddo 
      Acc = this%m_cG*Acc

! -- step with the compensated sommation !! 
!!      XPplan = XPplan - tau * Acc
      do i=1,nplan
         Vhat= this%m_plan_coordXP(:,i)
         errp = this%m_plan_errXP(:,i)
         call  kahansum3(Vhat,errp, - tau*Acc(:,i))
         this%m_plan_coordXP(:,i) = Vhat
         this%m_plan_errXP(:,i) = errp
      end do
       
      end  subroutine syscircumnewt_pasB2  

!***********************************************************************
!> @brief execute le pas B1 de l'integrateur
!***********************************************************************
      subroutine syscircumnewt_pasB1(this, tau)
       use kahansum
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: tau !< pas de temps
       
       real(TREAL), dimension(3,this%m_plan_nb):: aux
       real(TREAL), dimension(3):: delt,err,V
       integer iv, ip, nplan
       real(TREAL) :: mstar ! somme des masses des 2 etoiles
       
!        X(:,ip) = X(:,ip) 
!     &             + tau*(aux(:,ip)-mpsmstar(ip)*XP(:,ip))


      mstar = this%m_summass_star
      nplan = this%m_plan_nb
      do iv = 1,3
         aux(iv,:)=dot_product(this%m_plan_mass(2:nplan)/mstar,         &
     &                         this%m_plan_coordXP(iv,2:nplan))
      enddo 

      do ip=2,nplan
         delt = tau*(aux(:,ip)-this%m_plan_mass(ip)/mstar               &
     &     *this%m_plan_coordXP(:,ip))

         V = this%m_plan_coordX(:,ip)
         err = this%m_plan_errX(:,ip)
         call  kahansum3(V,err, delt)
         this%m_plan_coordX(:,ip) = V
         this%m_plan_errX(:,ip) = err
      end do

      end  subroutine syscircumnewt_pasB1 

!***********************************************************************
!> @brief execute le pas B de l'integrateur = B1(tau)B2(tau)B1(tau) 
!***********************************************************************
      subroutine syscircumnewt_pasB(this, mdt)
       use kahansum
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
       ! write(*,*) "syscircumnewt_pasB", mdt
       call syscircumnewt_pasB1(this,mdt/REALC(2.E0))
       call syscircumnewt_pasB2(this,mdt)
       call syscircumnewt_pasB1(this,mdt/REALC(2.E0))

      end  subroutine syscircumnewt_pasB  

!***********************************************************************
!> @brief execute le pas C de l'integrateur
!***********************************************************************
      subroutine syscircumnewt_pasC(this, mdt)
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(in) :: mdt !<= cm*dt
       
       write(*,*) "syscircumnewt_pasC", mdt
       write(*,*) "pas de correcteur en heliocentrique"
       stop
       
      end  subroutine syscircumnewt_pasC  


!***********************************************************************
!> @brief calcule l'energie de this
!***********************************************************************
      subroutine syscircumnewt_energie(this, H)
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), intent(out) :: H  !< energie calculee
       
       call Energie_circum(this%m_plan_nb,this%m_cG,this%m_summass_star  &
     &    , this%m_star_mass, this%m_plan_mass, this%m_mu,               &
     &    this%m_plan_mpsbeta,this%m_plan_coordX,this%m_plan_coordXP,H) 
      ! write(*,*) "syscircumnewt_energie ",H
       
      end  subroutine syscircumnewt_energie  

!***********************************************************************
!> @brief calcule le moment cinetique de this
!***********************************************************************
      subroutine syscircumnewt_mom_cin(this, C)
       implicit none
       class(t_syscircumnewt), intent(inout):: this  !< dummy argument
       real(TREAL), dimension(1:3), intent(out) :: C  !< moment cinetique
       
       call mom_cin_circum(this%m_plan_nb,this%m_plan_mass,             &
     &    this%m_plan_coordX, this%m_plan_coordXP,C) 
       !write(*,*) "syscircumnewt_energie",C
       
      end  subroutine syscircumnewt_mom_cin  

!***********************************************************************
!> @brief calcul de l'energie en variables heliocentriques canoniques
!***********************************************************************
      subroutine Energie_circum(nplan,cG,mstar,m0,mpl,mu,mpsbeta,X,XP,   &
     &  H)
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes +1 (pour la 2nde etoile en 1)
      real(TREAL),dimension(3,nplan), intent(in) :: X !< position circumbinaire (v)
      real(TREAL),dimension(3,nplan), intent(in) :: XP !< vitesse circumbinaire (v^)
      real(TREAL),dimension(nplan), intent(in) :: mpl !< masse de la 2nde etoile (en 1) et des planetes (en 2:)
      real(TREAL),dimension(nplan), intent(in) :: mu !< G*(m0+m1) de la 2nde etoile (en 1) et G*(m0+m1+mp) des planetes (en 2:)
      real(TREAL),dimension(nplan), intent(in) :: mpsbeta !< masse/beta de la 2nde etoile (en 1) et des planetes (en 2:)
      real(TREAL), intent(in) :: mstar !< somme des masses des 2 etoiles
      real(TREAL), intent(in) :: m0 !< masse de la 1ere etoile
      real(TREAL), intent(in) :: cG !< constante de gravitation
      real(TREAL), intent(out) :: H  !< energie calculee
      real(TREAL) :: H0, T0,U0,T1,U1,U2, factm
      real(TREAL),dimension(3) :: XX
      
      real(TREAL) :: n0, n1, m1, usvi
      real(TREAL) :: vi_n0_v1,vi_n1_v1
      real(TREAL),dimension(3) :: v1, vi
      integer i, j

      m1 = mpl(1)
      n0 = -m1/mstar !=-m1/(m0+m1)
      n1 = m0/mstar !=m0/(m0+m1)
      v1(:) = X(:,1) ! 2nde etoile

      
      T1 = REALC(0.e0)
      U1 = REALC(0.e0)
      T0 = REALC(0.e0)
      U0 = REALC(0.e0)
      U2 = REALC(0.E0)
      ! write(*,*) 'masses ', mpl
      ! write(*,*) 'mpsbeta ',mpsbeta
      ! write(*,*) X,XP
      
      do i=1,nplan
         T0=T0+mpl(i)*mpsbeta(i)*dot_product(XP(:,i),XP(:,i))
         factm = mu(i)*(mpl(i)/mpsbeta(i))
         U0=U0+factm/sqrt(dot_product(X(:,i),X(:,i)))
      enddo
      H0 = T0/REALC(2.e0) - U0


      do i=2,nplan
         do j=i+1,nplan
          XX = X(:,i) - X(:,j)
          U1 = U1 + mpl(i)*mpl(j)/sqrt(dot_product(XX,XX))
          T1 = T1 + mpl(i)*mpl(j)*dot_product(XP(:,i),XP(:,j))
         end do
      end do
      
      ! seconde partie planetes-etoile
      do i=2,nplan
        vi = X(:,i)
        usvi = REALC(1.E0)/sqrt(dot_product(vi,vi))
        vi_n0_v1=REALC(1.E0)/sqrt(dot_product(vi-n0*v1,vi-n0*v1))
        vi_n1_v1=REALC(1.E0)/sqrt(dot_product(vi-n1*v1,vi-n1*v1)) 
        U2=U2+mpl(i)*m0*(vi_n0_v1-usvi)+mpl(i)*m1*(vi_n1_v1-usvi)
      end do
      
       H = H0 - cG*U1 + T1/mstar - cG*U2
      end subroutine Energie_circum

!***********************************************************************
!> @brief calcul du moment cinetique en variables circum-binaires
!***********************************************************************
      Subroutine mom_cin_circum(nplan, mpl,X,XC,C)
      implicit none
      integer,intent(in) :: nplan !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),dimension(nplan), intent(in) :: mpl !< masse des planetes et de la 2nde etoile
      real(TREAL),dimension(3), intent(out) :: C !< moment cinetique
      real(TREAL),dimension(3,nplan), intent(in) :: X  !< position circum-binaires
      real(TREAL),dimension(3,nplan), intent(in) :: XC  !< vitesse circum-binaires   
      integer :: i
      
      C = REALC(0e0)
      do i=1,nplan
         C(1) = C(1) + (X(2,i)*XC(3,i) - X(3,i)*XC(2,i))*mpl(i)
         C(2) = C(2) + (X(3,i)*XC(1,i) - X(1,i)*XC(3,i))*mpl(i)
         C(3) = C(3) + (X(1,i)*XC(2,i) - X(2,i)*XC(1,i))*mpl(i)
      end do         
      end Subroutine mom_cin_circum

      end module mod_syscircumnewt
