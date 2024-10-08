!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file elliptiDpart.f 
!!  \brief conversion de coordonnees elliptiques pour les particules
!!
! history : creation 26/08/2014
!***********************************************************************
#include "realc.f"
**********************************************************
!
!   Routines de transformation d'elements elliptique 
!   heliocentriques (canoniques et non canoniques)
!   en positions vitesses pour les particules
!**********************************************************
!   
!   npl  : nombre de planetes.
!   mpl(0:npl) : masse du Soleil (0) puis des planetes (1:npl).
!   cG : contante de la gravitation.  
!   cmu(1:npl) : coefficient mu du probleme de kepler.
!   ell_kh_nc(6:npl) : a,la,k,h,q,p pour chaque planete.
!                     calcule a partir des elements cononiques.
!   ell_kh_c(6:npl)  : a,la,k,h,q,p pour chaque planete calcule
!                      a partir des elements non cononiques.
!   ell_pi(6:npl) : a,e,i,la,pi,Omega pour chaque planete. 
!   ell_om(6:npl) : a,e,i,l,omega,Omega pour chaque planete.
!   Ph(3:npl) : positions heliocentriques des planetes.
!   Vh(3:npl) : vitesses heliocentriques des planetes. 
!   Vb(3:npl) : vitesses barycentriques des planetes. 
!   
*****************************************************************
       module mod_elliptidpart
       use mod_elliptid
       
       contains


      subroutine ell_hcan2PhVh_part(npl,mpl,cG,ell_kh_c, npart,          &
     &                  ell_kh_c_part, Phpart,Vhpart)
!*******************************************************
! passage des elements elliptiques (a,la,k,h,q,p) 
! aux positions et vitesse heliocentriques 
! canoniques pour les particules.
! base sur ell_hcan2PhVh.
!*******************************************************
       use mod_coordpart
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),intent(in) :: cG
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh_c
      integer, intent(in) :: npart !< nombre de particules
      real(TREAL),dimension(6,npart),intent(in) :: ell_kh_c_part !< elements elliptiques (a,la,k,h,q,p) des particules
      real(TREAL),dimension(3,npart),intent(out) :: Phpart !< positions heliocentriques canoniques pour les particules
      real(TREAL),dimension(3,npart),intent(out) :: Vhpart !< vitessses heliocentriques canoniques pour les particules

      real(TREAL),dimension(3,npl) :: Vb
      real(TREAL),dimension(3,npl) :: Ph
      real(TREAL),dimension(npl) :: cmu
      real(TREAL),dimension(3,npart) :: Vbpart
      real(TREAL) :: cmupart
      integer i

!------ Positions et vitesses heliocentriques pour les planetes
      do i=1,npl
         cmu(i) = cG*(mpl(0) + mpl(i))
         call ellipx1(ell_kh_c(:,i),cmu(i),Ph(:,i),Vb(:,i))
!------ facteur beta((i)/m(i)
         Vb(:,i) = mpl(0)/(mpl(0) + mpl(i))*Vb(:,i)
      end do
      
!------ Positions et vitesses heliocentriques pour les particules
         cmupart = cG*mpl(0)
      do i=1,npart
         call ellipx1(ell_kh_c_part(:,i),cmu(i),Phpart(:,i),Vbpart(:,i))
      end do
      
!------ vitesses heliocentriques
      call coord_b2h_part(npl,mpl,Vb, npart, Vbpart, Vhpart)
      end subroutine

      end module
