!*********************************************************
!
!   Routines de transformation des positions vitesses 
!   entre les systemes de coordonnees heliocentrique,
!   barycentrique et de Jacobi. 
!**********************************************************
!   
!   npl  : nombre de planetes.
!   mpl(0:npl) : masse du Soleil (0) puis des planetes (1:npl).
!   eta(0:npl) : eta(i) = sum(mpl(0:i))
!   Xh(3:npl) : positions ou vitesses heliocentriques des planetes.
!   Xb(3:npl) : positions ou vitesses barycentriques des planetes. 
!   Xj(3:npl) : positions ou vitesses de Jacobi des planetes. 
!
!
!  liste des routines :
!  coord_b2h(npl,mpl,Xb,Xh)
!  coord_h2b(npl,mpl,Xh,Xb)
!  coord_h2j(npl,mpl,eta,Xh,Xj)
!  coord_j2h(npl,mpl,eta,Xj,Xh)
!
!   17/05/2001  P. Robutel
!
!  version 1.1  : (jxl) on a rajoute les routines de plan_invariant.f
!  
!  pass_invariant
!  coord_ang
!  rot_3(teta,u,ru)
!  rot_1(teta,u,ru)
!  coor_inv(X,X_inv,teta1,teta3)
!
!*****************************************************************
#include "realc.f"
	 

      subroutine coord_vb2vh(npl,mpl,Vb,Vh)
!*******************************************************
! Transformation des vitesses barycentriques (Vb)
! en vitesses heliocentriques (Vh)
! npl: nombre de planetes 
! mpl: masses soleil+planetes
!
! 4/9/2004  J. Laskar
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(3,npl),intent(in):: Vb
      real(TREAL),dimension(3,npl),intent(out):: Vh
      real(TREAL),dimension(3) :: vb0
      real(TREAL) :: masstot
      integer i
!----------- calcul de la vitesse barycentrique du soleil
      masstot = sum(mpl)
      do i=1,3
         Vb0(i) = -dot_product(mpl(1:),Vb(i,:))/mpl(0)
      enddo
!------ 
      do i=1,3
        Vh(i,1:) = Vb(i,:)-Vb0(i)
      end do
      end     
   

      subroutine coord_vh2vb(npl,mpl,Vh,Vb)
!*******************************************************
! Transformation des vitesses heliocentriques (Vh)
! en vitesses barycentriques (Vb)
! npl: nombre de planetes 
! mpl: masses soleil+planetes
!
! 4/9/2004  J. Laskar
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(3,npl),intent(in):: Vh
      real(TREAL),dimension(3,npl),intent(out):: Vb
      real(TREAL),dimension(3) :: vb0
      real(TREAL) :: masstot
      integer i
!----------- calcul de la vitesse barycentrique du soleil
      masstot = sum(mpl)
      do i=1,3
         Vb0(i) = -dot_product(mpl(1:),Vh(i,:))/masstot
      enddo
!------ 
      do i=1,3
        Vb(i,1:) = Vh(i,:)+Vb0(i)
      end do
      end     
   


      subroutine coord_b2h(npl,mpl,Xb,Xh)
!*******************************************************
! Transformation d'une quantite (position ou vitesse) 
! barycentrique en heliocentrique
!
! 15/05/2001 : Philippe Robutel
!*******************************************************
      implicit none 
      integer, intent(in) :: npl
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(3,npl),intent(in):: Xb
      real(TREAL),dimension(3,npl),intent(out):: Xh
      integer i

      do i=1,3
!------ Jocelyn Couetdic 30/04/04 Ajout de '/mpl(0)'
        Xh(i,:) =  dot_product(mpl(1:),Xb(i,:))/mpl(0) + Xb(i,:)
      end do
      end     
   
   
      subroutine coord_h2b(npl,mpl,Xh,Xb)
!*******************************************************
! Transformation d'une quantite (position ou vitesse) 
! heliocentrique en barycentrique
!
! 15/05/2001 : Philippe Robutel
!*******************************************************
      implicit none 
      integer,intent(in) :: npl
      real(TREAL),dimension(0:npl),intent(in) :: mpl
      real(TREAL),dimension(3,npl),intent(in) :: Xh
      real(TREAL),dimension(3,npl),intent(out) :: Xb
      integer i
     
      do i=1,3
        Xb(i,:) = Xh(i,:) - dot_product(mpl(1:),Xh(i,:))/sum(mpl)
      end do
      end     
        
	 
      subroutine coord_h2j(npl,mpl,eta,Xh,Xj)
!*******************************************************
! Transformation d'une quantite (position ou vitesse) 
! heliocentrique en Jacobi
!
! 15/05/2001 : Philippe Robutel
!*******************************************************
      implicit none 
      integer,intent(in) :: npl
      real(TREAL),dimension(0:npl),intent(in) :: mpl,eta
      real(TREAL),dimension(3,npl),intent(in) :: Xh
      real(TREAL),dimension(3,npl),intent(out) :: Xj
      integer i
      real(TREAL),dimension(3):: SX,SXX
          
      Xj(:,1) = Xh(:,1)
      SX = mpl(1)*Xh(:,1)
      SXX = SX/eta(1)
      do i=2,npl
        Xj(:,i) = Xh(:,i) - SXX
        if (i.lt.npl) then
           SX = SX + mpl(i)*Xh(:,i)
           SXX = SX/eta(i)
        end if
      enddo
      end   
      
        
  
      subroutine coord_j2h(npl,mpl,eta,Xj,Xh)
!*******************************************************
! Transformation d'une quantite (position ou vitesse) 
! Jacobi en heliocentrique 
!
! 15/05/2001 : Philippe Robutel
!*******************************************************
      implicit none 
      integer,intent(in) :: npl
      real(TREAL),dimension(0:npl),intent(in) :: mpl,eta
      real(TREAL),dimension(3,npl),intent(in) :: Xj
      real(TREAL),dimension(3,npl),intent(out) :: Xh
      integer i
      real(TREAL),dimension(3):: SX
      
      Xh(:,1) = Xj(:,1)
      SX = mpl(1)*Xj(:,1)/eta(1)     
      do i=2,npl
        Xh(:,i) = Xj(:,i) + SX
        if (i.lt.npl) then
          SX = SX + mpl(i)*Xj(:,i)/eta(i)
        end if
      enddo
      end     

!----------------------------------------------------------------------*
!       routines de plan_invariant    *
!----------------------------------------------------------------------*

      subroutine pass_invar(nplan,X,Xp,C,X_inv,Xp_inv)
      implicit none
      integer,intent(in) :: nplan
      real(TREAL),dimension(3,nplan),intent(in) :: X,Xp 
      real(TREAL),intent(in) :: C(3) 
      real(TREAL),dimension(3,nplan),intent(out) :: X_inv,Xp_inv 
      integer :: i
      real(TREAL) :: teta,phi,teta1,teta3,CN(3)
      real(TREAL) :: pi = REALC(3.14159265358979e0)
      
!--------- Determination des angle de rotation

      CN(:) = C(:)/sqrt(dot_product(C(:),C(:))) 
      write(*,*)'pass_invar'
      write(*,*) C
      write(*,*) CN
      call coord_ang(CN,teta,phi)
      if (sin(teta).ge.REALC(0e0)) then
         teta1 = -teta
         teta3 = mod(-pi/REALC(2.e0)-phi,pi)
      else
         teta1 = teta
         teta3 = mod(pi/REALC(2.e0)-phi,pi)
      end if
      write(*,*)'teta,phi',teta,phi
      write(*,*)'teta1,teta3',teta1,teta3
!------- Passage en coord. du plan invariant
      do i=1,nplan
         call coor_inv(X(:,i),X_inv(:,i),teta1,teta3)
	 call coor_inv(Xp(:,i),Xp_inv(:,i),teta1,teta3)
      end do  
      write(*,*)'pass_invar fin'    
      end
      
      
      subroutine coord_ang(X,teta,phi)
!----------------------------------------------------------------------*
!        Determination des angles  teta et phi pour X donne par :      *
!        X = (sin(teta)*cos(phi),sin(teta)*sin(phi),cos(teta))         *
!----------------------------------------------------------------------*
!
      implicit none 
      real(TREAL),intent(in) :: X(1:3)
      real(TREAL),intent(out) :: teta,phi
      real(TREAL) :: test
      real(TREAL) :: eps = REALC(2.e0)*epsilon(REALC(1.e0))
      
      test = ABSR(X(3) - REALC(1e0))
      if (test.le.eps) then
         phi = atan2(REALC(-1e0),REALC(0e0))
         teta = REALC(0e0)
      else
         phi = atan2(X(2),X(1))
         if (ABSR(sin(phi)).lt.REALC(1e-1)) then
            teta = atan2(X(1)/cos(phi),X(3))
         else
            teta = atan2(X(2)/sin(phi),X(3))
         end if
      end if
      end

      subroutine rot_3(teta,u,ru)
!----------------------------------------------------------------------*
!        ru = R_3(teta)u                                               *
!----------------------------------------------------------------------*
      implicit none
      real(TREAL),intent(in) :: teta,u(3)
      real(TREAL),intent(out):: ru(3)
      real(TREAL) :: c,s
      
      c = cos(teta)
      s = sin(teta)
      ru(1) = c*u(1) - s*u(2)
      ru(2) = s*u(1) + c*u(2)
      ru(3) = u(3)
      end


      subroutine rot_1(teta,u,ru)
!----------------------------------------------------------------------*
!        ru = R_1(teta)u                                               *
!----------------------------------------------------------------------*
      implicit none
      real(TREAL),intent(in) :: teta,u(3)
      real(TREAL),intent(out) :: ru(3)
      real(TREAL) :: c,s

      c = cos(teta)
      s = sin(teta)
      ru(1) = u(1)
      ru(2) = c*u(2) - s*u(3)
      ru(3) = s*u(2) + c*u(3)
      end



      subroutine coor_inv(X,X_inv,teta1,teta3)
!----------------------------------------------------------------------*
!        Transformation des coordonnees initiales en plan invariant    *
!----------------------------------------------------------------------*
!
      implicit none
      real(TREAL),intent(in) :: X(3),teta1,teta3
      real(TREAL),intent(out) :: X_inv(3)
      real(TREAL) :: XX(3)     
      
      call rot_3(teta3,X,XX)
      call rot_1(teta1,XX,X_inv)
      end


