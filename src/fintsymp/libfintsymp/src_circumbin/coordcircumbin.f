!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file coordcircumbin.f 
!!  \brief routine de conversion pour les coordonnees circum-binaires
!!
! history : creation 19/11/2015
!***********************************************************************
#include "realc.f"

       module mod_coordcircumbin
       
       contains

      subroutine circum_calc_mu_mpsbeta(npl,cG,mpl,mu,mpsbeta)
!*******************************************************
!!> calcule mu et m_p/beta_p pour les quantites circum-binaires
!! calcul particulier pour la 2nde etoile (voir le code)
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),intent(in) :: cG
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes (en 2:)
      real(TREAL),dimension(1: npl),intent(out) :: mu !< G*(m+m_p) de la 2nde etoile (en 1) et des planetes (en 2:)
      real(TREAL),dimension(1: npl),intent(out) :: mpsbeta !< m_p/beta_p de la 2nde etoile (en 1) et des planetes (en 2:)
      real(TREAL) :: mstar ! somme des masses des 2 etoiles

         mstar = mpl(0)+mpl(1)
         ! pour la 2nde etoile
         mu(1) = cG*mstar
         mpsbeta(1) = mstar/mpl(0) 

         ! pour les planetes
         mu(2:npl) = cG*(mstar + mpl(2:npl))
         mpsbeta(2:npl) = (mpl(2:)+mstar)/mstar 

      end subroutine circum_calc_mu_mpsbeta

      subroutine circum_calc_mu_betasmp(npl,cG,mpl,mu,betasmp)
!*******************************************************
!!> calcule mu et beta_p/m_p pour les quantites circum-binaires
!! calcul particulier pour la 2nde etoile (voir le code)
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),intent(in) :: cG
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes (en 2:)
      real(TREAL),dimension(1:npl),intent(out) :: mu !< G*(m+m_p) de la 2nde etoile (en 1) et des planetes (en 2:)
      real(TREAL),dimension(1:npl),intent(out) :: betasmp !< /beta_p/m_p de la 2nde etoile (en 1) et des planetes (en 2:)
      real(TREAL) :: mstar ! somme des masses des 2 etoiles

         mstar = mpl(0)+mpl(1)
         ! pour la 2nde etoile
         mu(1) = cG*mstar
         betasmp(1) = mpl(0)/mstar 

         ! pour les planetes
         mu(2:npl) = cG*(mstar + mpl(2:npl))
         betasmp(2:npl) = mstar/(mpl(2:npl)+mstar)

      end subroutine circum_calc_mu_betasmp



      subroutine circum_ell_hcan2vvhat(npl,mpl,cG,ell_kh_c,V,Vhat)
!*******************************************************
!!> passage des elements elliptiques (a,la,k,h,q,p) canoniques 
!! aux positions/vitesses (V,Vhat) circum-binaires 
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!! les elements canoniques sont les elements associes a (V, Vhat)
!*******************************************************
       use mod_elliptid
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),intent(in) :: cG
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes en (2:)
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh_c !< elements elliptiques canoniques de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(out) :: V,Vhat !< position-vitesse de la 2nde etoile et des planetes
      real(TREAL),dimension(3,npl) :: Vc
      real(TREAL),dimension(1:npl) :: cmu ! G*(m+m_p) 
      real(TREAL),dimension(1:npl) :: betasmp ! beta_p/m_p
      integer i

      call circum_calc_mu_betasmp(npl, cG, mpl, cmu, betasmp)
      do i=1,npl
         call ellipx1(ell_kh_c(:,i),cmu(i),V(:,i),Vc(:,i))
!------ facteur beta(i)/m(i)
         Vhat(:,i) = betasmp(i)*Vc(:,i)
      end do
      
      end subroutine circum_ell_hcan2vvhat

      subroutine circum_vvhat2ell_hcan(npl,mu,mpsbeta,V,Vhat,ell_kh_c,       &
     &                 arret)
!*******************************************************
!!> passage positions/vitesses (V,Vhat) circum-binaires a 
!! des elements elliptiques (a,la,k,h,q,p) canoniques  
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!! les elements canoniques sont les elements associes a (V, Vhat)
!*******************************************************
       use mod_elliptid
       use mod_arret
       implicit none 
       integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
       real(TREAL),dimension(1:npl),intent(in) :: mu ! G*(m+m_p)  de la 2nde etoile (en 1) et des planetes en (2:)
       real(TREAL),dimension(1:npl),intent(in) :: mpsbeta ! m_p/beta_p de la 2nde etoile (en 1) et des planetes en (2:)
       real(TREAL),dimension(6,npl),intent(out) :: ell_kh_c !< elements elliptiques canoniques de la 2nde etoile (en 1) et des planetes en (2:)
       real(TREAL),dimension(3,npl),intent(in) :: V,Vhat !< position-vitesse de la 2nde etoile et des planetes
       type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
       real(TREAL),dimension(3,npl) :: Vc
       integer i
       
       i = 0 
       do while ((i.lt.npl).and.(arret%m_stop.eq.0))
          i = i + 1
          Vc(:,i) = mpsbeta(i)*Vhat(:,i)
          call xyzkhqp1(V(:,i),Vc(:,i),mu(i), ell_kh_c(:,i), arret)
       end do
       if (arret%m_stop.ne.0) then
          call arret%set_body(i) 
       end if     
             
      end subroutine circum_vvhat2ell_hcan


      subroutine circum_PhVh2vvhat(npl,mpl,Ph,Vh,V,Vhat)
!*******************************************************
!!> passage des positions et vitesses heliocentriques (Ph,Vh) 
!! aux positions/vitesses (V,Vhat) circum-binaires 
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(in) :: Ph,Vh !< position-vitesse heliocentrique de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(out) :: V,Vhat !< position-vitesse circumbinaire de la 2nde etoile et des planetes
      real(TREAL) :: massstar ! masse totale des etoiles
      real(TREAL) :: mtotale ! masse totale du systeme
      real(TREAL), dimension(3) :: ndr0 ! =dr_0/dt
      integer i, j

      massstar = mpl(0)+mpl(1)
      mtotale = sum(mpl(0:npl))
      
      !2nde etoile
      V(:,1)=Ph(:,1)
      Vhat(:,1)=mpl(0)/massstar*Vh(:,1)
      
      ! les planetes
      do j=1, 3
        ndr0(j)=dot_product(mpl(1:npl),Vh(j,1:npl))/mtotale
      enddo  
      do i=2,npl
         V(:,i)=Ph(:,i)-mpl(1)/massstar*Ph(:,1)
         Vhat(:,i)=Vh(:,i)-ndr0
      end do
      
      end subroutine circum_PhVh2vvhat

      subroutine circum_vvhat2PhVh(npl,mpl,V,Vhat,Ph,Vh)
!*******************************************************
!!> passage des positions/vitesses (V,Vhat) circum-binaires 
!! aux positions et vitesses heliocentriques (Ph,Vh) 
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!*******************************************************
       use mod_elliptid
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(in) :: V,Vhat !< position-vitesse circumbinaire de la 2nde etoile et des planetes
      real(TREAL),dimension(3,npl),intent(out) :: Ph,Vh !< position-vitesse heliocentrique de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL) :: massstar ! masse totale des etoiles
      real(TREAL), dimension(3) :: dv0
      integer i, j

      massstar = mpl(0)+mpl(1)
      
      !2nde etoile
      Ph(:,1)=V(:,1)
      Vh(:,1)=massstar/mpl(0)*Vhat(:,1)
      
      ! les planetes
      do j=1, 3
        dv0(j)=dot_product(mpl(2:npl),Vhat(j,2:npl))/massstar
      enddo  
      do i=2,npl
         Ph(:,i)=V(:,i)+mpl(1)/massstar*V(:,1)
         Vh(:,i)=Vhat(:,i)+mpl(1)/mpl(0)*Vhat(:,1)+dv0
      end do
      
      end subroutine circum_vvhat2PhVh

      subroutine circum_vdot2vhat(npl,mpl, betasmp, Vdot,Vhat)
!*******************************************************
!!> passage des vitesses dV/dt  (derivee de V) 
!! aux vitesses Vhat (=Vtitde/m_i) circum-binaires 
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes en (2:)
      real(TREAL),dimension(1:npl),intent(in) :: betasmp !< beta_p/m_p de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(in) :: Vdot !< vitesse (derivee de V) de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(out) :: Vhat !< vitesse circumbinaire de la 2nde etoile et des planetes
      real(TREAL) :: massstar ! masse totale des etoiles
      real(TREAL) :: mtotale ! masse totale du systeme
      real(TREAL), dimension(3) :: dv0 
      integer i, j

      massstar = mpl(0)+mpl(1)
      mtotale = sum(mpl(0:npl))
      
      !2nde etoile
      Vhat(:,1)=betasmp(1)*Vdot(:,1)
      
      ! les planetes
      do j=1, 3
        dv0(j)=dot_product(mpl(2:npl),Vdot(j,2:npl))/mtotale
      enddo  
      do i=2,npl
         Vhat(:,i)=Vdot(:,i)-dv0
      end do
      
      end subroutine circum_vdot2vhat

      subroutine circum_vhat2vdot(npl,mpl, mpsbeta, Vhat, Vdot)
!*******************************************************
!!> passage des vitesses Vhat (=Vtitde/m_i) 
!! aux  dV/dt  (derivee de V)  circum-binaires 
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!*******************************************************
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes en (2:)
       real(TREAL),dimension(1:npl),intent(in) :: mpsbeta ! m_p/beta_p de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(in) :: Vhat !< vitesse circumbinaire de la 2nde etoile et des planetes
      real(TREAL),dimension(3,npl),intent(out) :: Vdot !< vitesse (derivee de V) de la 2nde etoile (en 1) et des planetes en (2:)
      
      real(TREAL) :: massstar ! masse totale des etoiles
      real(TREAL), dimension(3) :: dv0
      integer i, j

      massstar = mpl(0)+mpl(1)

      !2nde etoile
      Vdot(:,1)=mpsbeta(1)*Vhat(:,1)
      
      ! les planetes
      do j=1, 3
        dv0(j)=dot_product(mpl(2:npl),Vhat(j,2:npl))/massstar
      enddo  
      do i=2,npl
         Vdot(:,i)=Vhat(:,i)+dv0
      end do
      
      end subroutine circum_vhat2vdot

      subroutine circum_ell_h2vvhat(npl,mpl,cG,ell_kh,V,Vhat)
!*******************************************************
!!> passage des elements elliptiques (a,la,k,h,q,p) non canoniques 
!! aux positions/vitesses (V,Vhat) circum-binaires 
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!! les elements non canoniques sont les elements associes a (V, Vdot)
!*******************************************************
       use mod_elliptid
      implicit none 
      integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
      real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes en (2:)
      real(TREAL), intent(in) :: cG !< constante de gauss en UA, an
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh !< elements elliptiques non canoniques de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL),dimension(3,npl),intent(out) :: V,Vhat !< position-vitesse de la 2nde etoile et des planetes
      real(TREAL),dimension(3,npl) :: Vdot
      integer i
      real(TREAL),dimension(1:npl) :: mu ! G*(m+m_p)  de la 2nde etoile (en 1) et des planetes en (2:)
      real(TREAL),dimension(1:npl) :: betasmp ! beta_p/m_p de la 2nde etoile (en 1) et des planetes en (2:)


      call circum_calc_mu_betasmp(npl, cG, mpl, mu, betasmp)
      
      do i=1,npl
         call ellipx1(ell_kh(:,i),mu(i),V(:,i),Vdot(:,i))
      end do
      call circum_vdot2vhat(npl,mpl,betasmp,Vdot,Vhat)
      
      end subroutine circum_ell_h2vvhat

      subroutine circum_vvhat2ell_h(npl,mpl,mu,mpsbeta,V,Vhat,          &
     &                 ell_kh, arret)
!*******************************************************
!!> passage positions/vitesses (V,Vhat) circum-binaires a 
!! des elements elliptiques (a,la,k,h,q,p) canoniques  
!!         les coordonnees des corps sont exprimees en :
!!              * 2eme etoile  reperee par rapport a la premiere etoile
!!              * les planetes sont reperees par rapport au barycentre des 2 etoiles
!! les elements non canoniques sont les elements associes a (V, Vdot)
!*******************************************************
       use mod_elliptid
       use mod_arret
       implicit none 
       integer, intent(in) :: npl !< nombre de planetes +1 (pour la 2nde etoile)
       real(TREAL),dimension(0:npl),intent(in) :: mpl !< masses des 2 etoiles (en 0:1) et des planetes en (2:)
       real(TREAL),dimension(1:npl),intent(in) :: mu ! G*(m+m_p)  de la 2nde etoile (en 1) et des planetes en (2:)
       real(TREAL),dimension(1:npl),intent(in) :: mpsbeta ! m_p/beta_p de la 2nde etoile (en 1) et des planetes en (2:)
       real(TREAL),dimension(6,npl),intent(out) :: ell_kh !< elements elliptiques non canoniques de la 2nde etoile (en 1) et des planetes en (2:)
       real(TREAL),dimension(3,npl),intent(in) :: V,Vhat !< position-vitesse de la 2nde etoile et des planetes
       type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
       real(TREAL),dimension(3,npl) :: Vdot
       integer i
       
       call circum_vhat2vdot(npl, mpl, mpsbeta, Vhat, Vdot)
       i = 0 
       do while ((i.lt.npl).and.(arret%m_stop.eq.0))
          i = i + 1
          call xyzkhqp1(V(:,i),Vdot(:,i),mu(i), ell_kh(:,i), arret)
       end do
       if (arret%m_stop.ne.0) then
          call arret%set_body(i) 
       end if     
             
      end subroutine circum_vvhat2ell_h

      end  module mod_coordcircumbin
