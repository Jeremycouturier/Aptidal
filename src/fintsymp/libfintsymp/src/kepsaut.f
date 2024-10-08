!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!  \file kepsaut.f 
!  \brief Implementation of :
!    kepsaut
!
!    ce fichier ne doit rien contenir de specifique a un projet 
!
!
!
! history : creation 2014/03/19
!
!   17/05/2001  P. Robutel, J. Laskar, F. Joutel
!   revise le 3/9/2004 pour inclure kepsaut et annexes (jxl)
!
!
! Modification: by A. Farres (10/01/2011)
!    1 -> change real(8) by real(TREAL)
!    2 -> change dabs by  abs to be agreement with (TREAL)
!    3 -> EPS = 2.23d-16 changed for EPS = epsilon(REALC(1.e0)) 
! NOTE: line 509:69... on peut pas changer avec REALC .. il me done ++ 
! des errors dans la compilation !
!
! Modification: by A. Farres (17/01/2011)
!   two new routines
!    1)  subroutine keps_new_aut(X,XP,cmu,dt,Xn, XDn) 
!			 like kepsaut but it return Xn and XDn the increment 
!			 in position and velocities for a given step TAU = dt. 
!		 i.e. X = X + Xn*tau // XP = XP + XDn*tau  (in pasA and pasAH)
!    2)  subroutine prods_new(a,c,X,b,d,Y,R,S) 
!	       like prods but modified to be used by { keps_new_aut } 
!			 it now returns R = aX+bY   and   S = c*X+d*Y
!
! Modification: by A. Farres (20/10/2011)
!    1) We no longer use the approximation of sin(x)-x by the polynomial 
!     as it gives bad results when we consider extended arithmetics. 
!     (It might be interesting to get a better approximation)
!    2) In subroutine keps_new_aut we have changed the way of evaluating 
!      a1, a2, a3, a4 to get less rounding errors (changes proposed by Mickael)
!      BEFORE: 
!         a1= cm1/(REALC(2.e0)-R*V2)
!			a2=dt +(smx)*sqa*a/smu
!			a3=-s/((REALC(2.e0)-R*V2)*(REALC(1.e0)-ce*c+se*s)*(sqa*a))*smu
!			a4=(a2*a3-a1)/(REALC(1.e0)+a1)
!		NOW: 
!			a1= cm1/rv2 
!			a3=-s/(rv2*e*(a*sqa))*smu
!			a4= cm1/e
!			a2 = (a1+a4+a1*a4)/a3;
!
! Modification: by A. Farres (25/10/2011)
!   1 -> change ABS by ABSR that is defined in realc.f (this way we 
!        can deal with TREAL = 8,10 and 16 )
!   2 -> change DMAX1 by DMAX1 so that it soportes REAL(16) 
!
! Modification: by M. gastineau (29/11/2018)
!        appel de newt2018 au lieu de newt2006
!        meilleure condition d'arret de l'algo de newton
!***********************************************************************
#include "realc.f"

      module mod_kepsaut
         use mod_arret
         
         private NEWT
        
      contains
        
      SUBROUTINE NEWT(DM,A,B,X,C,S,EXP,CM1,SMX, arret)
!*---------------------------------------------------
!* RESOLUTION DE L'ÈQUATION F(X)=X-A*SIN(X)-B*COS(X)+B-DM=0
!* SORTIE DE X, C=COS(X),S=SIN(X),EXP=1-A*C-B*S,CM1,SMX
!   rajout du test d'erreur de P.R (jxl) 
!-----------------------------------------------------
      IMPLICIT NONE
      type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
       REAL(TREAL),INTENT(IN):: DM,A,B
       REAL(TREAL),INTENT(OUT):: X,C,S,EXP,CM1,SMX
       REAL(TREAL) B0,B1,F,F1
       REAL(TREAL) D1,EPS,DIFF
       INTEGER I,NMAX
       !EPS=2.22D-16
       EPS = epsilon(REALC(1.e0))
!       NMAX=10
       NMAX=15

!* DEPART METHODE D'ORDRE 2
       F=DM
       F1=REALC(1.e0)-A
       D1=F/F1
       X=D1
!* METHODE DE NEWTON
       DO I=1,NMAX
           C=COS(X)
           S=SIN(X)
           B0=A*S+B*C
           B1=A*C-B*S
           F=B-DM+X-B0
           F1=REALC(1.e0)-B1
           DIFF=F/F1
           X=X-DIFF
!*  NORMALEMEMT X ET DIFF SONT PETITS
!          IF (ABS(DIFF)/MAX1(1.D0,ABS(X)).LT.1.5D0*EPS) GOTO 100
!        IF ((ABS(DIFF)/MAX1((REALC(1.e0)),ABS(X))).LT. 
        IF ((ABS(DIFF)/MAX1((REALC(1.e0)),ABS(X))).LT. 
     & ((REALC(1.5e0))*EPS)) GOTO 100
        END DO
        WRITE(*,*) 'ERREUR FATALE DANS KEPSAUT:NEWT'
        WRITE(*,*) 'X,DIFF',X,DIFF
!------------------- gestion de l'erreur 
        call arret%set_error(-3)
!       STOP     
100     C=COS(X)
        S=SIN(X)
! CALCUL DE L'EXPRESSION EXP    
        EXP=F1-DIFF*B0
! CALCUL DE COS(X)-1    
        CM1=-REALC(2.e0)*SIN(X*REALC(0.5e0))**2        
! CALCUL DE SIN(X)-X
!        IF (ABS(X).LT.REALC(0.2e0)) THEN               
        IF (0.EQ.1) THEN               
         SMX = (REALC(2623.e0)/REALC(4929724800.e0)*X**9-
     &          REALC(97.e0)/REALC(1053360.e0)*X**7
     &   +REALC(59.e0)/REALC(9880.e0)*X**5-X**3/(REALC(6.e0)))
     &     /(1+REALC(7.e0)/REALC(494.e0)*X**2+
     &    REALC(23.e0)/REALC(326040.e0)*X**4)
         ELSE
          SMX=S-X
         ENDIF
        END subroutine
        
      SUBROUTINE NEWT2006(DM,A,B,X,C,S,EXP,CM1,SMX, arret)
*---------------------------------------------------
* RESOLUTION DE L'ÈQUATION F(X)=X-A*SIN(X)-B*COS(X)+B-DM=0
* SORTIE DE X, C=COS(X),S=SIN(X),EXP=1-A*C-B*S,CM1,SMX
!   rajout du test d'erreur de P.R (jxl) 
!   10/04/06 : erreur releve a 4eps
!   05/05/06 : F. Joutel : dans le cas de fortes excentricites,
!       les calculs se font mal lorsque A est trËs proche de L
!       on passe alors en quadruple prÈcision
!   20/10/06 : J. Couetdic : ajout d un compteur pour les
!       passages en quadruple
!-----------------------------------------------------
      IMPLICIT NONE
       type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
       REAL(TREAL),INTENT(IN):: DM,A,B
       REAL(TREAL),INTENT(OUT):: X,C,S,EXP,CM1,SMX
       REAL(TREAL) B0,B1,F,F1
       REAL(TREAL) D1,EPS,DIFF
       real(16) :: CQ,SQ,B0Q,B1Q,FQ,F1Q,DIFFQ,XQ
       real(16) :: AQ,BQ,DMQ
       INTEGER I,NMAX
       EPS = epsilon(REALC(1.e0))
       NMAX=15
!* DEPART METHODE D'ORDRE 2
       F=DM
       F1=REALC(1.e0)-A
       D1=F/F1
       X=D1
!* METHODE DE NEWTON
       DO I=1,NMAX
           C=COS(X)
           S=SIN(X)
           B0=A*S+B*C
           B1=A*C-B*S
           F=B-DM+X-B0
           F1=REALC(1.e0)-B1
           DIFF=F/F1
           X=X-DIFF
!*  NORMALEMEMT X ET DIFF SONT PETITS
!          IF (ABS(DIFF)/MAX1(1.D0,ABS(X)).LT.1.5D0*EPS) GOTO 100
!        IF ((ABS(DIFF)/MAX1((REALC(1.e0)),ABS(X))).LT. 
        IF ((ABS(DIFF)/MAX1((REALC(1.e0)),ABS(X))).LT. 
     & ((REALC(1.5e0))*EPS)) GOTO 100
        END DO

        AQ=A
        BQ=B
        DMQ=DM
        XQ=X
91     DO I=1,NMAX
           CQ=COS(XQ)
           SQ=SIN(XQ)
           B0Q=AQ*SQ+BQ*CQ
           B1Q=AQ*CQ-BQ*SQ
           FQ=BQ-DMQ+XQ-B0Q
           F1Q=1.Q0-B1Q
           DIFFQ=FQ/F1Q
           XQ=XQ-DIFFQ
*  NORMALEMEMT X ET DIFF SONT PETITS
           IF (ABS(DIFFQ)/MAX(1.Q0,ABS(XQ)).LT.1.5*EPS) THEN
             X=XQ
             GOTO 100
           ENDIF
           X=XQ
        END DO        
        WRITE(*,*) 'ERREUR FATALE DANS KEPSAUT:NEWT'
        WRITE(*,*) 'X,DIFF',X,DIFF
!------------------- gestion de l'erreur 
        call arret%set_error(-3)
!       STOP     
100     C=COS(X)
        S=SIN(X)
! CALCUL DE L'EXPRESSION EXP        
        EXP=F1-DIFF*B0
! CALCUL DE COS(X)-1        
        CM1=-2.D0*SIN(X*0.5D0)**2         
! CALCUL DE SIN(X)-X
!        IF (ABS(X).LT.0.2D0) THEN               
        IF (0.EQ.1) THEN               
         SMX = (2623.D0/4929724800.D0*X**9-
     &          97.D0/1053360.D0*X**7
     &   +59.D0/9880.D0*X**5-X**3/6.D0)
     &     /(1+7.D0/494.D0*X**2+
     &    23.D0/326040.D0*X**4)
         ELSE
          SMX=S-X
         ENDIF
        END SUBROUTINE NEWT2006

!-----------------------------------------------------
!> prend une position-vitesse, calcule les elements
!! elliptiques, ajoute N*dt a la longitude moyenne
!! puis revient au une position-vitesse 
!!   Les variables X et XD sont modifies par la routine
!!
!! si une erreur s'est produite, arret met a jour le code d'erreur 
!! mais pas le temps    
!!
!! -> 17/1/2011 A. FARRES
!!   modification de KEPSAUT ! 
!!    now the ruttine returns DX an DPX that is what needs to 
!!    be advanced outside to do the compensated somation !
!------------------------------------------------------    
      subroutine keps_new_aut(X,XP,cmu,dt,Xn, XDn, arret)
      implicit none     
      type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
      real(TREAL),intent(in) :: cmu,dt
      real(TREAL),intent(in) :: X(3),XP(3)
      real(TREAL),intent(out) :: Xn(3),XDn(3)
      real(TREAL) smu
      real(TREAL):: r,v2,xv,rv2,a,ce,se,dm,c,s,e
      real(TREAL):: a1,a2,a3,a4,delta,sqa   
      real(TREAL) ::cm1,smx
      smu=sqrt(cmu)
      R=sqrt(X(1)**2+X(2)**2+X(3)**2)
      v2=(XP(1)**2+XP(2)**2+XP(3)**2)/cmu
      XV=(X(1)*XP(1)+X(2)*XP(2)+X(3)*XP(3))/smu
      rv2=(REALC(2.e0)-R*V2)
      a=R/rv2
      sqa=sqrt(a)
      ce=R*V2-REALC(1.e0)   
      se=XV/sqa
      dm=smu*(dt/(sqa*a))
!      call newt(dm,ce,se,delta,c,s)
!      call newt(dm,ce,se,delta,c,s,e,cm1,smx, arret)
!      call newt2006(dm,ce,se,delta,c,s,e,cm1,smx, arret)
      call newt2018(dm,ce,se,delta,c,s,e,cm1,smx, arret)
! on enleve l'identite 
!!      a1= cm1/(REALC(2.e0)-R*V2)  
      a1= cm1/rv2  
!!      a2=dt +(smx)*sqa*a/smu
!      a3=-s/((REALC(2.e0)-R*V2)*(REALC(1.e0)-ce*c+se*s)*(sqa*a))*smu
!      a3=-s/(rv2*e*(a*sqa))*smu
      a3=-s/(R*e*sqa)*smu
!      a4= cm1/(REALC(1.e0)-ce*c+se*s)
      a4= cm1/e
!      a4=(a2*a3-a1)/(REALC(1.e0)+a1)      
      a2 = (a1+a4+a1*a4)/a3;
!      call prods(a1,a3,X,a2,a4,Xp,Xn,XDn)
      Xn = a1*X + a2*Xp
      XDn = a3*X + a4*Xp
      end subroutine keps_new_aut

      SUBROUTINE NEWT2018(DM,A,B,X,C,S,EXP,CM1,SMX, arret)
*---------------------------------------------------
* RESOLUTION DE L'ÈQUATION F(X)=X-A*SIN(X)-B*COS(X)+B-DM=0
* SORTIE DE X, C=COS(X),S=SIN(X),EXP=1-A*C-B*S,CM1,SMX
!   rajout du test d'erreur de P.R (jxl) 
!   10/04/06 : erreur releve a 4eps
!   05/05/06 : F. Joutel : dans le cas de fortes excentricites,
!       les calculs se font mal lorsque A est trËs proche de L
!       on passe alors en quadruple prÈcision
!   20/10/06 : J. Couetdic : ajout d un compteur pour les
!       passages en quadruple
!-----------------------------------------------------
      IMPLICIT NONE
       type(t_arret), intent(inout) :: arret !< cause de l'arret (non modifie si aucune erreur)
       REAL(TREAL),INTENT(IN):: DM,A,B
       REAL(TREAL),INTENT(OUT):: X,C,S,EXP,CM1,SMX
       REAL(TREAL) B0,B1,F,F1, F2
       REAL(TREAL) D1,EPS,DIFF
       real(16) :: CQ,SQ,B0Q,B1Q,FQ,F1Q,DIFFQ,XQ
       real(16) :: AQ,BQ,DMQ
       INTEGER I,NMAX
       EPS = epsilon(REALC(1.e0))
       NMAX=15
!* DEPART METHODE D'ORDRE 2
       F=DM
       F1=REALC(1.e0)-A
       D1=F/F1
       X=D1
!* METHODE DE NEWTON
       DO I=1,NMAX
           C=COS(X)
           S=SIN(X)
           B0=A*S+B*C
           B1=A*C-B*S
           F=B-DM+X-B0
           F1=REALC(1.e0)-B1
           F2=B0
           DIFF=F/F1
           X=X-DIFF
!*  NORMALEMEMT X ET DIFF SONT PETITS
           IF (ABS(DIFF).LT.(REALC(1.2e0)*EPS*ABS(X))) GOTO 100
        END DO

        AQ=A
        BQ=B
        DMQ=DM
        XQ=X
91     DO I=1,NMAX
           CQ=COS(XQ)
           SQ=SIN(XQ)
           B0Q=AQ*SQ+BQ*CQ
           B1Q=AQ*CQ-BQ*SQ
           FQ=BQ-DMQ+XQ-B0Q
           F1Q=1.Q0-B1Q
           DIFFQ=FQ/F1Q
           XQ=XQ-DIFFQ
*  NORMALEMEMT X ET DIFF SONT PETITS
           IF (ABS(DIFFQ).LT.1.5Q0*EPS*ABS(XQ)) THEN
             X=XQ
             GOTO 100
           ENDIF
           X=XQ
        END DO        
        WRITE(*,*) 'ERREUR FATALE DANS KEPSAUT:NEWT'
        WRITE(*,*) 'X,DIFF',X,DIFF
!------------------- gestion de l'erreur 
        call arret%set_error(-3)
!       STOP     
100     C=COS(X)
        S=SIN(X)
! CALCUL DE L'EXPRESSION EXP        
        EXP=F1-DIFF*B0
! CALCUL DE COS(X)-1        
        CM1=-REALC(2.E0)*SIN(X*REALC(0.5E0))**2         
! CALCUL DE SIN(X)-X
        SMX=S-X
        END SUBROUTINE NEWT2018
      


      end  module mod_kepsaut

