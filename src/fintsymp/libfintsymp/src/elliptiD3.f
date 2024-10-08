#include "realc.f"
**********************************************************
!
!   Routines de transformation d'elements elliptique 
!   heliocentriques (canoniques et non canoniques)
!   en positions vitesses.
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
!   Controle de l'ellipticite des orbites : 
!      pl_stop(1) = 0 : tout est elliptique,
!                     ==> pl_stop(2) = 0
!      pl_stop(1) = i :  la ieme planete n'a pas une
!                        orbite elliptique.
!                     ==> pl_stop(2) = -4
!
!
!  liste des routines :
!  keplkh2(l,k,h,a)
!  ellipx1(ell_kh,cmu,x,xp)
!  xyzkhqp1(x,xd,cmu,ell_kh)
!  PhVh2ell_hcan(npl,mpl,cG,Ph,Vh,ell_kh_c)
!  ell_hcan2PhVh(npl,mpl,cG,ell_kh_c,Ph,Vh)
!  ell_h2hcan(npl,mpl,cG,ell_kh_nc,ell_kh_c)
!  ell_hcan2h(npl,mpl,cG,ell_kh_c,ell_kh_nc)
!  ell_pi2ell_kh(npl,ell_pi,ell_kh)
!  ell_kh2ell_pi(npl,ell_kh,ell_pi)
!  ell_kh2ell_om(npl,ell_kh,ell_om)
!  ell_om2ell_kh(npl,ell_om,ell_kh)
!
!
! 20/03/2014 M. Gastineau : modification por la gestion des arrets
*****************************************************************
       module mod_elliptid
       use mod_arret
       
       contains
       
       subroutine keplkh2(l,k,h,a)
!---------------------------------------------------
! Resolution de l'equation de Kepler 
! pour les variables L=M+ pi et A=E+pi
! precision 2*(eps-mach) pour des excentricites
! inferieures a 0.3
! temps de calculs :
!     pour la premiere iteration (recherche du depart) 
!                     2 (cos sin)
!                    11 (* /)
!                     9 (+ -)
!    pour la deuxieme iteration
!                     2 (cos sin)
!                    25 (* /)
!                    14 (+ -)
! F. Joutel 26/4/99
! 26/1/01 appel eventuel a la methode de Newton
! 14/05/2001 P. Robutel : passage en fortran90
!---------------------------------------------------
       implicit none
       real(TREAL),intent(in) :: l,k,h
       real(TREAL),intent(out) :: a
       real(TREAL) sa,ca,se,ce,fa,f1a,f2a,f3a,f4a,f5a
       real(TREAL) d1,d2,d3,d4,d5
       real(TREAL) :: eps = REALC(2.e0)*epsilon(REALC(1.e0))
       integer ::i
       integer :: imax = 20
! depart methode d'ordre 3
		 a=l       
       ca=cos(a)
       sa=sin(a)
       se=k*sa-h*ca
       ce=k*ca+h*sa
       fa=a-se-l
       f1a=REALC(1.e0)-ce
       f2a=se/REALC(2.e0)
       f3a=ce/REALC(6.e0)
       d1=-fa/f1a
       d2=-fa/(f1a-d1*f2a)
       d3 =-fa/(f1a+d2*(f2a+d2*f3a))
       a=a+d3
! methode d'ordre 6
       ca=cos(a)
       sa=sin(a)
       se=k*sa-h*ca
       ce=k*ca+h*sa
       fa=a-se-l
       f1a=REALC(1.e0)-ce
       f2a=se/REALC(2.e0)
       f3a=ce/REALC(6.e0)
       f4a=-se/REALC(24.e0)
       f5a=-ce/REALC(120.e0)
       d1=-fa/f1a
       d2=-fa/(f1a-d1*f2a)
       d3=-fa/(f1a+d2*(f2a+d2*f3a))
       d4=-fa/(f1a+d3*(f2a+d3*(f3a+d3*f4a)))
       d5=-fa/( f1a+d4*(f2a+d4*(f3a+d4*(f4a+d4*f5a))))
       a=a+d5
!       return
! on teste la precision obtenue
       i=0
100    i=i+1
       ca=cos(a)
       sa=sin(a)
       se=k*sa-h*ca
       fa=a-se-l
       ce=k*ca+h*sa
       f1a=REALC(1.e0)-ce
       d1=-fa/f1a
! si la precison n'est pas bonne, on continue les calculs
! en iterant la methode d'ordre 1       
       if (ABSR(d1)/max1(REALC(1.e0),ABSR(a)).gt.eps) then
         if (i.eq.imax) then
!          write(*,*) 'erreur fatale dans elliptid:keplkh2'
!          write(*,*) 'erreur :',ABSR(d1)/dmax1(REAL(1.e0),ABSR(a))
           return
         endif
         a=a+d1
         goto 100
       else
         return 
       endif 
       end subroutine





      subroutine ellipx1(ell_kh,cmu,x,xp)
!--------------------------------------------------------
!   Conversion elements elliptiques en position-vitesse 
! J. Laskar
! F. Joutel 26/4/99
! 14/05/2001 P. Robutel : passage en fortran90
!---------------------------------------------------------     
      implicit none
      integer i,j
      real(TREAL),intent(in) :: ell_kh(6),cmu
      real(TREAL),intent(out) :: x(3),xp(3)
      real(TREAL) :: a,l,k,h,q,p,n,phi,ki,psi,sf,cf
      real(TREAL) :: rlmf,umrsa,asr,rsa,rna2sr,f
      real(TREAL),dimension(3,2) :: rot
      real(TREAL),dimension(2) :: tx1(2),tx1t(2)
      a=ell_kh(1) 
      l=ell_kh(2) 
      k=ell_kh(3) 
      h=ell_kh(4) 
      q=ell_kh(5) 
      p=ell_kh(6) 
      n=sqrt(cmu/(a**3)) 
      phi=sqrt(1-k**2-h**2) 
      ki=sqrt(1-q**2-p**2) 
      psi=1/(1+phi) 
 
!---- matrice de rotation ----------------------------------------------

      rot(1,1)=1-2*p**2 
      rot(1,2)=2*p*q 
      rot(2,1)=2*p*q 
      rot(2,2)=1-2*q**2
      rot(3,1)=-2*p*ki 
      rot(3,2)= 2*q*ki 
 
!---- calcul de la longitude excentrique f -----------------------------
!---- f = anomalie excentrique e + longitude du periapse omegapi -------
 

      call keplkh2(l,k,h,f) 
      sf    =sin(f) 
      cf    =cos(f) 
      rlmf  =-k*sf+h*cf 
      umrsa =k*cf+h*sf 
      asr   =1/(1-umrsa) 
      rsa   =1/asr 
      rna2sr=n*a*asr 
 
!---- calcul de tx1 et tx1t --------------------------------------------
 
      tx1(1) =a*(cf-psi*h*rlmf-k) 
      tx1(2) =a*(sf+psi*k*rlmf-h) 
      tx1t(1)=rna2sr*(-sf+psi*h*umrsa) 
      tx1t(2)=rna2sr*( cf-psi*k*umrsa) 

!---- calcul de xyz ----------------------------------------------------
 
      do i=1,3 
         x(i) =0 
         xp(i)=0 
         do j=1 ,2 
            x(i)  =x(i)  +rot(i,j)*tx1(j) 
            xp(i) =xp(i) +rot(i,j)*tx1t(j)
          end do
       end do
       end subroutine


       SUBROUTINE XYZKHQP1(X,XD,CMU,ELL, arret) 
!--------------------------------------------------------
!   Conversion position-vitesse en elements elliptiques
!   convient pour les mouvements elliptiques non
!   rectilignes ( C <> 0) et les inclinaisons differentes
!   de Pi ( cos(i) <> -1)
!   gere le cas : cos(i)==-1 et e<>0
! J. Laskar (1986-2004)  (voir notes de cours)
! version revisee par rapport a la version originale de 1986 
! amelioration dans toutes les precisions 
! par rapport a la version XYZKHQP1 modifiee par Fred.
! BUG corrige en longitude quand e est nul
! precision moins bonne (mais encore OK)
!  pour tres grandes inclinaisons ( i > 1 radian)
!  v. 1.00  (jxl 11/06/2004)
!  ELL : a, lambda, k, h, q, p
!  rajout des test arret  18/8/2004
! modification des arrets M. Gastineau 20/03/2014
!---------------------------------------------------------     
      implicit none
      REAL(TREAL), intent (in) ::  X,XD,CMU
      REAL(TREAL), intent (out) ::  ELL
      type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
      DIMENSION X(3),XD(3),ELL(6) 
      REAL(TREAL)  R,V2,RV,C1,C2,C3,DC,CC,AA,aux0
      REAL(TREAL)  a11,a12,a21,a22,c11,c12,c21,c22,k1,h1,k2,h2,K,H
      REAL(TREAL)  V,SMU,FAC1,FAC2,b12,b22,sinF,cosF,F,USQA,aux1
      REAL(TREAL)  AUX,AUXP, EPS
      DIMENSION V(3)
      
!------------ normalisation des vitesses
      SMU=SQRT(CMU)
      V(1)=XD(1)/SMU
      V(2)=XD(2)/SMU
      V(3)=XD(3)/SMU
!------------ quantitees utilies 
      R=SQRT(X(1)**2+X(2)**2+X(3)**2) 
      V2=V(1)**2+V(2)**2+V(3)**2 
      RV=(X(1)*V(1)+X(2)*V(2)+X(3)*V(3))
      C1=X(2)*V(3)-X(3)*V(2) 
      C2=X(3)*V(1)-X(1)*V(3) 
      C3=X(1)*V(2)-X(2)*V(1)
      CC=C1**2+C2**2+C3**2
      DC= SQRT(CC)
!------------ demi grand axe a
      AA = R/(REALC(2.e0)-R*V2)
      ELL(1)=AA
      if (AA.le.(REALC(0e0))) then
!------  cas non elliptique
         call arret%set_error(-4)
         goto 9999
      endif
      USQA= sqrt(REALC(2.e0)/R-V2)
!------------ q,p
!      aux0 = sqrt(2*(CC+DC*C3))
      aux0 = 2*(CC+DC*C3)
      if (aux0.lt.(REALC(0e0))) then
!------  cas non elliptique
         call arret%set_error(-4)
         goto 9999
      endif
      if (aux0.eq.REALC(0e0)) then
!------  cas inclinaison =pi
         ELL(5)=REALC(1.e0)
         ELL(6)=REALC(0.e0)
!------------ k,h
         a11= V2-REALC(1.e0)/R
         a21=-RV 
         K = a11*X(1)+a21*V(1)
         H = -(a11*X(2)+a21*V(2))
         ELL(3)=K
         ELL(4)=H
!------------ calcul de la longitude
         EPS = epsilon(REALC(1.E0))
         IF (ABS(K)+ABS(H) .GT. EPS) THEN     
           AUX=RV/SQRT(AA)   
           AUXP=REALC(1.e0)-R/AA     
           F= ATAN2(K*AUX+H*AUXP,K*AUXP-H*AUX)
           ELL(2)=F-AUX 
          else
           ELL(2)=1D99
           call arret%set_error(-4)
           goto 9999
          endif 
       else
!------  cas inclinaison <>pi
       aux0 = sqrt(aux0)
       ELL(5)=-C2/aux0
       ELL(6)= C1/aux0 
!------------ k,h
!------------ coefs de  e * M^-1
      a11= V2-REALC(1.e0)/R
      a12= RV/(R*DC) 
      a21=-RV 
      a22= DC-R/DC 
!------------ calcul de (h,k)  
      c11=X(1)*a11+V(1)*a21 
      c12=X(1)*a12+V(1)*a22 
      c21=X(2)*a11+V(2)*a21 
      c22=X(2)*a12+V(2)*a22 
      FAC1=C1/(DC+C3) 
      FAC2=C2/(DC+C3) 
      k1= c11-FAC1*(X(3)*a11+V(3)*a21) 
      h1=-c12+FAC1*(X(3)*a12+V(3)*a22) 
      h2= c21-FAC2*(X(3)*a11+V(3)*a21) 
      k2= c22-FAC2*(X(3)*a12+V(3)*a22) 
      H =(h1+h2)/REALC(2.e0)
      K =(k1+k2)/REALC(2.e0) 
      ELL(3)=K
      ELL(4)=H
!------------ calcul de la longitude
      b12=V(1)-FAC1*V(3)
      b22=V(2)-FAC2*V(3)
      aux1 = (R*V2-REALC(1.e0))/(REALC(1.e0)+DC*USQA)
      sinF=-b12*R*USQA +H*aux1
      cosF= b22*R*USQA +K*aux1 
      F=atan2(sinF,cosF) 
      ELL(2)=F-RV*USQA
      endif   
9999  continue
      END subroutine


        
                
      subroutine PhVb2ell_hcan(npl,mpl,cG,Ph,Vb,ell_kh_c, arret)
!*******************************************************
! passage positions heliocentriques
! et vitesse barycentriques
! aux  elements elliptiques heliocentriques canoniques
! (a,la,k,h,q,p) 
!
!   En sortie, arret contient un erreur si probleme non elliptique
!   indice de la premiere planete non elliptique est stocke dans arret.
!
! 4/9/2004 J. Laskar  
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),intent(in) :: cG
      type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(3,npl) :: Ph 
      real(TREAL),dimension(6,npl),intent(out) :: ell_kh_c    
      real(TREAL),dimension(3,npl) :: Vb
      real(TREAL),dimension(npl) :: cmu
      integer i

      i = 0
!------ elements heliocentriques canoniques
      do while ((i.lt.npl).and.(arret%m_stop.eq.0))
        i = i + 1
        cmu(i) = cG*(mpl(0) + mpl(i))
!------ facteur m(i)/beta(i)
         Vb(:,i) = (REALC(1e0) + mpl(i)/mpl(0))*Vb(:,i)
         call xyzkhqp1(Ph(:,i),Vb(:,i),cmu(i),ell_kh_c(:,i), arret)
      end do
      if (arret%m_stop.ne.0) then
         call arret%set_body(i) 
      end if     
      end subroutine


      subroutine PhVh2ell_hcan(npl,mpl,cG,Ph,Vh,ell_kh_c, arret)
!*******************************************************
! passage positions et vitesse heliocentriques 
! aux  elements elliptiques heliocentriques canoniques
! (a,la,k,h,q,p) 
!
!   En sortie, arret contient un erreur si probleme non elliptique
!   indice de la premiere planete non elliptique est stocke dans arret.
!
! 17/05/2001 : P. Robutel
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),intent(in) :: cG
      type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(3,npl) :: Ph,Vh
      real(TREAL),dimension(6,npl),intent(out) :: ell_kh_c    
      real(TREAL),dimension(3,npl) :: Vb
      real(TREAL),dimension(npl) :: cmu
      integer i

      i = 0
!------ vitesses barycentriques
      call coord_h2b(npl,mpl,Vh,Vb)
!------ elements heliocentriques canoniques
      do while ((i.lt.npl).and.(arret%m_stop.eq.0))
        i = i + 1
        cmu(i) = cG*(mpl(0) + mpl(i))
!------ facteur m(i)/beta(i)
         Vb(:,i) = (REALC(1e0) + mpl(i)/mpl(0))*Vb(:,i)
         call xyzkhqp1(Ph(:,i),Vb(:,i),cmu(i),ell_kh_c(:,i), arret)
      end do
      if (arret%m_stop.ne.0) then
         call arret%set_body(i) 
      end if     
      end subroutine
   


      subroutine ell_hcan2PhVh(npl,mpl,cG,ell_kh_c,Ph,Vh)
!*******************************************************
! passage des elements elliptiques (a,la,k,h,q,p) 
! aux positions et vitesse heliocentriques 
! heliocentriques canoniques 
! correction Jxl 2004
!*******************************************************
! passage des positions et vitesse heliocentriques 
! aux elements elliptiques (a,la,k,h,q,p) 
! heliocentriques canoniques 
!
! 17/05/2001 : P. Robutel
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),intent(in) :: cG
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh_c
      real(TREAL),dimension(3,npl),intent(out) :: Ph,Vh
      real(TREAL),dimension(3,npl) :: Vb
      real(TREAL),dimension(npl) :: cmu
      integer i

!------ Positions et vitesses heliocentriques
      do i=1,npl
         cmu(i) = cG*(mpl(0) + mpl(i))
         call ellipx1(ell_kh_c(:,i),cmu(i),Ph(:,i),Vb(:,i))
!------ facteur beta((i)/m(i)
         Vb(:,i) = mpl(0)/(mpl(0) + mpl(i))*Vb(:,i)
      end do
!------ vitesses heliocentriques
      call coord_b2h(npl,mpl,Vb,Vh)
      end subroutine

 
      subroutine ell_h2hcan(npl,mpl,cG,ell_kh_nc,ell_kh_c, arret)
!*******************************************************
! passage des elements elliptiques (a,la,k,h,q,p) 
! heliocentriques (non canoniques) 
! aux elements heliocentriques canoniques
!
!   En sortie, arret contient un erreur si probleme non elliptique
!   indice de la premiere planete non elliptique est stocke dans arret.
!
! 16/05/2001 : P. Robutel
! 3/9/04  correction erreur jxl
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),intent(in) :: cG
      type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh_nc
      real(TREAL),dimension(6,npl),intent(out) :: ell_kh_c    
      real(TREAL),dimension(3,npl) :: Ph,Vh,Vb
      real(TREAL),dimension(npl) :: cmu
      integer i

!------ Positions et vitesses heliocentriques
      do i=1,npl
         cmu(i) = cG*(mpl(0) + mpl(i))
         call ellipx1(ell_kh_nc(:,i),cmu(i),Ph(:,i),Vh(:,i))
      end do
!------ vitesses barycentriques
      call coord_h2b(npl,mpl,Vh,Vb)
!------ elements heliocentriques canoniques
      do i=1,npl
!------ facteur m(i)/beta(i)
!         Vb(:,i) = ((REALC(1e0) + mpl(i)/mpl(0))*Vb(:,i)
         Vb(:,i) = (mpl(0) + mpl(i)/mpl(0))*Vb(:,i)
         call xyzkhqp1(Ph(:,i),Vb(:,i),cmu(i),ell_kh_c(:,i), arret)
      end do 
      end subroutine
   


      subroutine ell_hcan2h(npl,mpl,cG,ell_kh_c,ell_kh_nc, arret)
!*******************************************************
! passage des elements elliptiques (a,la,k,h,q,p) 
! heliocentriques canoniques 
! aux elements heliocentriques non canoniques 
!
!
!   En sortie, arret contient un erreur si probleme non elliptique
!   indice de la premiere planete non elliptique est stocke dans arret.
!
! 16/05/2001 : P. Robutel
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),intent(in) :: cG
      type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh_c
      real(TREAL),dimension(6,npl),intent(out) :: ell_kh_nc    
      real(TREAL),dimension(3,npl) :: Ph,Vh,Vb
      real(TREAL),dimension(npl) :: cmu
      integer i

!------ Positions et vitesses heliocentriques
      do i=1,npl
         cmu(i) = cG*(mpl(0) + mpl(i))
         call ellipx1(ell_kh_c(:,i),cmu(i),Ph(:,i),Vb(:,i))
!------ facteur beta((i)/m(i)
         Vb(:,i) = mpl(0)/(mpl(0) + mpl(i))*Vb(:,i)
      end do
!------ vitesses heliocentriques
      call coord_b2h(npl,mpl,Vb,Vh)
!------ elements heliocentriques non canoniques
      i = 0
      do while ((i.lt.npl).and.(arret%m_stop.eq.0))
        i = i + 1
         call xyzkhqp1(Ph(:,i),Vh(:,i),cmu(i),ell_kh_nc(:,i), arret)
      end do 
      if (arret%m_stop.ne.0) then
         call arret%set_body(i) 
      end if     
      end subroutine   
   




      subroutine ell_kh2ell_pi(npl,ell_kh,ell_pi)
!*******************************************************
! passage des elements elliptiques (a,la,k,h,q,p)  
! aux (a,e,i,la,pi,Omega) 
!
! 17/05/2001 : P. Robutel
!
!*******************************************************
      implicit none
      integer, intent(in) :: npl
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh
      real(TREAL),dimension(6,npl),intent(out) :: ell_pi 
      real(TREAL),dimension(npl) :: aux
      real(TREAL),parameter :: eps = epsilon(REALC(1.e0))
!      eps = 2.26d-16
!----- a et lambda 
      ell_pi(1,:) = ell_kh(1,:)
      ell_pi(4,:) = ell_kh(2,:)
!----- excentricite et pi
      ell_pi(2,:) = sqrt(ell_kh(3,:)**2 + ell_kh(4,:)**2)
      where (ell_pi(2,:).le.eps)
         ell_pi(2,:) = REALC(0e0)
         ell_pi(5,:) = REALC(0e0)
      elsewhere
         ell_pi(5,:) = atan2(ell_kh(4,:),ell_kh(3,:))
      end where
!----- inclinaison et Omega
      aux(:) = sqrt(ell_kh(5,:)**2 + ell_kh(6,:)**2)
      where (aux(:).le.eps)
         ell_pi(3,:) = REALC(0e0)
         ell_pi(6,:) = REALC(0e0)
      elsewhere
         ell_pi(3,:) = REALC(2e0)*asin(aux(:)) 
         ell_pi(6,:) = atan2(ell_kh(6,:),ell_kh(5,:))
      end where      
      end subroutine


      subroutine ell_pi2ell_kh(npl,ell_pi,ell_kh)
!*******************************************************
! passage des elements elliptiques (a,e,i,la,pi,Omega) 
! aux (a,la,k,h,q,p)  
!
! 17/05/2001 : P. Robutel
!
!*******************************************************
      implicit none
      integer, intent(in) :: npl
      real(TREAL),dimension(6,npl),intent(in) :: ell_pi
      real(TREAL),dimension(6,npl),intent(out) :: ell_kh 
      real(TREAL),dimension(npl) :: aux
      
      ell_kh(1,:) = ell_pi(1,:)
      ell_kh(2,:) = ell_pi(4,:)
      ell_kh(3,:) = ell_pi(2,:)*cos(ell_pi(5,:))
      ell_kh(4,:) = ell_pi(2,:)*sin(ell_pi(5,:))
      aux(:) = sin(ell_pi(3,:)/REALC(2e0))
      ell_kh(5,:) = aux*cos(ell_pi(6,:))
      ell_kh(6,:) = aux*sin(ell_pi(6,:))
      end  subroutine
      
      
      subroutine ell_kh2ell_om(npl,ell_kh,ell_om)
!*******************************************************
! passage des elements elliptiques (a,la,k,h,q,p)  
! aux (a,e,i,M,om,Omega) 
!
! 17/05/2001 : P. Robutel
!
!*******************************************************
      implicit none
      integer, intent(in) :: npl
      real(TREAL),dimension(6,npl),intent(in) :: ell_kh
      real(TREAL),dimension(6,npl),intent(out) :: ell_om 
      real(TREAL),dimension(6,npl) :: ell_pi 
      real(TREAL),parameter :: eps = epsilon(REALC(1.e0))

      call ell_kh2ell_pi(npl,ell_kh,ell_pi)
      ell_om(1:3,:) = ell_pi(1:3,:)
      ell_om(6,:) = ell_pi(6,:)
      ell_om(5,:) = ell_pi(5,:) - ell_pi(6,:)
      ell_om(4,:) = ell_pi(4,:) - ell_pi(5,:)
      end subroutine


      subroutine ell_om2ell_kh(npl,ell_om,ell_kh)
!*******************************************************
! passage des elements elliptiques (a,e,i,M,om,Omega) 
! aux (a,la,k,h,q,p)  
!
! 17/05/2001 : P. Robutel
!
!*******************************************************
      implicit none
      integer, intent(in) :: npl
      real(TREAL),dimension(6,npl),intent(in) :: ell_om
      real(TREAL),dimension(6,npl),intent(out) :: ell_kh 
      real(TREAL),dimension(6,npl) :: ell_pi
      
      ell_pi(1:3,:) = ell_om(1:3,:)
      ell_pi(6,:) = ell_om(6,:)
      ell_pi(5,:) = ell_om(5,:) + ell_om(6,:)
      ell_pi(4,:) = ell_om(4,:) + ell_pi(5,:)
      call ell_pi2ell_kh(npl,ell_pi,ell_kh)
      end subroutine
 
      end module
