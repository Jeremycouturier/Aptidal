!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file dopri8t.f 
!!  \brief version modernisee de dopri8.f mais requiert ode.f
!!
!
! history : creation 07/12/2015
!***********************************************************************
#include "realc.f"

!***********************************************************************
!> module pour l'integration avec DOPRI 
!***********************************************************************
      module mod_dopri
      use mod_ode

#if TREAL==8    
!***********************************************************************
!> @class t_dopri
!! integration avec dopri :
!! il faut appeller :
!! 1. start : pour le remplissage du vecteur d'etat  et l'initialisation de dopri
!! 2. run : pour le remplissage du vecteur d'etat 
!!
!***********************************************************************
      type, extends(t_schemaode) :: t_dopri
          
          real(TREAL), dimension(:), allocatable :: m_y !< vecteur d'etat
          
       contains         
          procedure :: start => dopri_start  !< fonction appellee lors du demarrage
          procedure :: run => dopri_run  !< fonction appellee lors d'un pas normal
          procedure :: check_error => dopri_check_error  ! controle si une erreur s'est produite
       !   procedure, private :: COEFST !< coefficient de dopri8
      end type t_dopri  
#endif
     
      contains

#if TREAL==8    
!***********************************************************************
! fonction pour demarrer une integration  au temps tstart  avec l'integrateur
!***********************************************************************
       subroutine dopri_start(this, system, tstart, dt)
       implicit none
       class(t_dopri), intent(inout) :: this
       class(t_sysode), intent(inout) :: system ! systeme qui est integre
       real(TREAL), intent(in) :: dt ! pas d'integration
       real(TREAL), intent(in) :: tstart ! temps initial
        COMMON/ASD/IPTR,NFPTR,NERROR
        integer IPTR,NFPTR,NERROR
       
       ! remplissage du vecteur d'etat
       call system%allocate_y(this%m_y)
       call system%fill_to_y(this%m_y, tstart)
       ! initialise les variables de dopri
       IPTR = 0
       NFPTR = 6
       NERROR = 0
       
       end subroutine dopri_start

!***********************************************************************
! fonction pour executer un pas normal d'integration avec DOPRI
!***********************************************************************
      subroutine dopri_run(this, system, tstart, tend, dt)
      implicit none
      class(t_dopri), intent(inout) :: this
      class(t_sysode), intent(inout) :: system ! systeme qui est integre
      real(TREAL), intent(in) :: dt ! pas d'integration
      real(TREAL), intent(in) :: tstart ! temps initial
      real(TREAL), intent(in) :: tend ! temps final = tstart+dt
      
      call DOPRI8(size(this%m_y),system, tstart, this%m_y,tend,           &
     &    this%m_eps, this%m_hmax, this%m_h)
         call system%fill_from_y(this%m_y)
      
      end subroutine dopri_run

!***********************************************************************
!> @brief controle si une erreur s'est produite. 
!! retourne true en cas d'erreur
!***********************************************************************
       function dopri_check_error(this) result(ret)
        implicit none
        class(t_dopri), intent(inout) :: this  !< dummy argument
        logical ret ! code de retour
        COMMON/ASD/IPTR,NFPTR,NERROR
        integer IPTR,NFPTR,NERROR
        ret = .false.
        if (NERROR.ne.0) then
          ret = .true.        
        endif         
       end function dopri_check_error


        SUBROUTINE DOPRI8(N,system,X,Y,XEND,EPS,HMAX,H)
C---------------------------------------------------------
C       NUMERICAL SOLUTION OF A SYSTEM OF FIRST ORDER
C       ORDINARY DIFFERRENTIAL EQUATIONS Y'=F(X,Y).
C       THIS IS AN EMBEDDED RUNGE-KUTTA METHOD OF ORDER (7)8
C       DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL).
C       C.F. SECTION II.6
C
C       PAGE 435 HAIRER, NORSETT & WANNER
C
C       INPUT PARAMETERS
C       ----------------
C       N            DIMENSION OF THE SYSTEM ( N.LE.51)
C       FCN          NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
C                    FIRST DERIVATIVE F(X,Y):
C                      SUBROUTINE FCN(N,X,Y,F)
C                      REAL*8 X,Y(N),F(N)
C                      F(1)=....  ETC.
C       X            INITIAL X-VALUE
C       XEND         FINAL X-VALUE (XEND-X POSITIVE OR NEGATIVE)
C       Y(N)         INITIAL VALUES FOR Y
C       EPS          LOCAL TOLERANCE
C       HMAX         MAXIMAL STEPSIZE
C       H            INITIAL STEPSIZE GUESS
C       OUTPUT PARAMETERS
C       -----------------
C       Y(N) SOLUTION AT XEND
C    
C       EXTERNAL SUBROUTINE (TO BE SUPPLIED BY THE USER)
C       -------------------
C       SOLDOPRI       THIS SUBROUTINE IS CALLED AFTER EVERY
C                    STEP
C                       SUBROUTINE SOLDOPRI(NR,X,Y,N)
C                       REAL*8 X,Y(N)
C                    FURNISHES THE SOLUTION Y AT THE NR-TH
C                    GRID-POINT X (THE INITIAL VALUE IS CON-
C                    SIDERED AS THE FIRST GRID-POINT).
C                    SUPPLIED A DUMMY SUBROUTINE, IF THE SOLUTION
C                    IS NOT DESIRED AT THE INTERMEDIATE POINTS.
C--------------------------------------------------------------------
C   Routine reentrante (Utilisation de la variable FIRST)
C   Modifie le 9/5/95 pour une recuperation de l'echec (NERROR)
C   et la prise en compte du pas XEND-X (variable LAST)
C   Modif le 24/5/95 pour un bug sur l'emploi de LAST
C   Modif le 31/5/95 pour plus de modularite et suppression des 
C   boucles numerotees
C   Modif le 26/6/95 pour prendre en compte l'u.d.p. local dans la recherche 
C   du dernier pas
C   (*m/4)
C   ASD version 1.04 (26/6/95)
C--------------------------------------------------------------------
        IMPLICIT REAL(TREAL) (A-H,O-Z)
        class(t_sysode), intent(inout) :: system ! systeme qui est integre
        REAL(TREAL) K1(51),K2(51),K3(51),K4(51),K5(51),K6(51),K7(51),
     1   Y(N),Y1(51)
        LOGICAL REJECT,FIRST,LAST
        DATA FIRST /.TRUE./
        COMMON/DOP8ST/NFCN,NSTEP,NACCPT,NREJCT
        COMMON/ASD/IPTR,NFPTR,NERROR
        SAVE FIRST  
C-------COMMON DOP8ST CAN BE USED FOR STATISTICS
C-------   NFCN     NUMBER OF FUNCTION EVALUATIONS
C-------   NSTEP    NUMBER OF COMPUTEF STEPS 
C-------   NACCPT   NUMBER OF ACCEPTED STEPS
C-------   NREJCT   NUMBER OF REJECTED STEPS
C*******   NERROR   ERROR NUMBER
C                      + 0 : OK
C                      + 1 : MORE THAN NMAX STEPS
C                      + 2 : STEP IS TOO SMALL
C
        COMMON/COEFRK/C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,
     &A21,A31,A32,A41,A43,A51,A53,A54,A61,A64,A65,A71,A74,A75,A76,
     &A81,A84,A85,A86,A87,A91,A94,A95,A96,A97,A98,A101,A104,A105,A106,
     &A107,A108,A109,A111,A114,A115,A116,A117,A118,A119,A1110,A121,
     &A124,A125,A126,A127,A128,A129,A1210,A1211,A131,A134,A135,A136,
     &A137,A138,A139,A1310,A1311,B1,B6,B7,B8,B9,B10,B11,B12,B13,
     &BH1,BH6,BH7,BH8,BH9,BH10,BH11,BH12
         DATA NMAX_REFMIN/2000/, UROUND/2.23d-16/,FACSEC/0.03D0/
         DATA NMAX_REFMAX/20000000/
C-------  NMAX    MAXIMAL NUMBER OF STEPS
C-------  UROUND  SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.0
C-------          (TO BE ADAPTED BY THE USER)
**--- Modif *m/4 (31/5/95) pour adapter le programme a des pas plus petits
C-------  FACSEC  Facteur de securite pour ne pas prendre un pas
C-------          plus petit que (u.d.p. local)/FACSEC (FACSEC<=1)  
      HINPUT=H
      HFINAL=H
      XINPUT=X
      IF (FIRST) THEN
          CALL COEFST
          FIRST=.FALSE.
      ENDIF
      POSNEG=DSIGN(1.D0,XEND-X)
C------- INITIAL PREPARATIONS

      HMAXLOC=DABS(HMAX)
      H=DMIN1(DMAX1(1.D-10,DABS(H)),HMAXLOC)
      H=DSIGN(H,POSNEG)
      EPSLOC=DMAX1(EPS,13.D0*UROUND)
      REJECT=.FALSE.
      LAST=.FALSE.
      NACCPT=0
      NREJCT=0       
      NFCN=0
      NSTEP=0
      NERRTEMP = 0
!---- controle du nombre de pas max
      NMAX = ABS((XEND-X)/(H/NMAX_REFMIN));
      if (NMAX>NMAX_REFMAX) NMAX = NMAX_REFMAX
      if (NMAX<NMAX_REFMIN) NMAX = NMAX_REFMIN

      CALL SOLDOPRI(NACCPT,X,Y,N)
C-------BASIC INTEGRATION STEP
1     CONTINUE
      IF (NSTEP.GT.NMAX.OR.X+FACSEC*H.EQ.X) GOTO 79
**------modif *m/4 (21/6/95) pour respecter l'u.d.p. local
      GLOC=DMAX1(DMIN1(DABS(XEND),DABS(X)),1.D0)
      IF( DABS(XEND-X).LE. UROUND*GLOC) GO TO 100
**      IF((X-XEND)*POSNEG+UROUND.GT.0D0) GO TO 100
**------end modif
**------modif *m/4 (9/5/95) pour prendre en compte le dernier pas
      IF((X+H-XEND)*POSNEG.GE.0.D0) then
         H=XEND-X
         LAST=.TRUE.
      ENDIF
**- -----end modif
      CALL system%fcn(N,X,Y,K1)
**------Modif *m/4 (9/5/95) Le test est rejete au cas de rejet du pas
**       pour eviter de passer dessus lors du dernier pas
** 2       if (nstep.gt.nmax.or.x+0.03D0*h.eq.x) goto 79
**------End Modif
2     CONTINUE
      NSTEP=NSTEP+1
C-----THE FIRST 9 STAGES
      DO I=1,N
         Y1(I)=Y(I)+H*A21*K1(I)
      END DO
      CALL system%fcn(N,X+C2*H,Y1,K2)
      DO I=1,N
         Y1(I)=Y(I)+H*(A31*K1(I)+A32*K2(I))
      END DO        
      CALL system%fcn(N,X+C3*H,Y1,K3)
      DO I=1,N
         Y1(I)=Y(I)+H*(A41*K1(I)+A43*K3(I))
      END DO
      CALL system%fcn(N,X+C4*H,Y1,K4)
      DO I=1,N
         Y1(I)=Y(I)+H*(A51*K1(I)+A53*K3(I)+A54*K4(I))
      END DO        
      CALL system%fcn(N,X+C5*H,Y1,K5)
      DO I=1,N
         Y1(I)=Y(I)+H*(A61*K1(I)+A64*K4(I)+A65*K5(I))
      END DO
      CALL system%fcn(N,X+C6*H,Y1,K6)
      DO I=1,N
         Y1(I)=Y(I)+H*(A71*K1(I)+A74*K4(I)+A75*K5(I)+A76*K6(I))
      END DO
      CALL system%fcn(N,X+C7*H,Y1,K7)
      DO I=1,N
         Y1(I)=Y(I)+H*(A81*K1(I)+A84*K4(I)+A85*K5(I)+A86*K6(I)
     &         +A87*K7(I))
      END DO       
      CALL system%fcn(N,X+C8*H,Y1,K2)
      DO I=1,N
         Y1(I)=Y(I)+H*(A91*K1(I)+A94*K4(I)+A95*K5(I)+A96*K6(I)
     &         +A97*K7(I)+A98*K2(I))
      END DO
      CALL system%fcn(N,X+C9*H,Y1,K3)
      DO I=1,N
         Y1(I)=Y(I)+H*(A101*K1(I)+A104*K4(I)+A105*K5(I)+A106*K6(I)
     &         +A107*K7(I)+A108*K2(I)+A109*K3(I))
      END DO
C------- COMPUTE INTERMEDIATE SUMS TO SAVE MEMORY
      DO I=1,N
         Y11S=A111*K1(I)+A114*K4(I)+A115*K5(I)+A116*K6(I)+A117*K7(I)
     &         +A118*K2(I)+A119*K3(I)
         Y12S=A121*K1(I)+A124*K4(I)+A125*K5(I)+A126*K6(I)+A127*K7(I)
     &         +A128*K2(I)+A129*K3(I)
         K4(I)=A131*K1(I)+A134*K4(I)+A135*K5(I)+A136*K6(I)+A137*K7(I)
     &         +A138*K2(I)+A139*K3(I)
         K5(I)=B1*K1(I)+B6*K6(I)+B7*K7(I)+B8*K2(I)+B9*K3(I)
         K6(I)=BH1*K1(I)+BH6*K6(I)+BH7*K7(I)+BH8*K2(I)+BH9*K3(I)
         K2(I)=Y11S
         K3(I)=Y12S
      END DO
C-----THE LAST 4 STAGES
      CALL system%fcn(N,X+C10*H,Y1,K7)
      DO I=1,N
         Y1(I)=Y(I)+H*(K2(I)+A1110*K7(I))
      END DO
      CALL system%fcn(N,X+C11*H,Y1,K2)
      XPH=X+H
      DO I=1,N
         Y1(I)=Y(I)+H*(K3(I)+A1210*K7(I)+A1211*K2(I))
      END DO     
      CALL system%fcn(N,XPH,Y1,K3)
      DO I=1,N
         Y1(I)=Y(I)+H*(K4(I)+A1310*K7(I)+A1311*K2(I))
      END DO          
      CALL system%fcn(N,XPH,Y1,K4)
      NFCN=NFCN+13
      DO I=1,N
         K5(I)=Y(I)+H*(K5(I)+B10*K7(I)+B11*K2(I)+B12*K3(I)+B13*K4(I))
         K6(I)=Y(I)+H*(K6(I)+BH10*K7(I)+BH11*K2(I)+BH12*K3(I))       
      END DO
C-----ERROR ESTIMATION
      ERR=0.D0
      DO I=1,N
         DENOM=DMAX1(1.D-6,DABS(K5(I)),DABS(Y(I)),2.D0*UROUND/EPSLOC)
         ERR=ERR+((K5(I)-K6(I))/DENOM)**2
      END DO
      ERR=DSQRT(ERR/DFLOAT(N))
C-----COMPUTATION OF HNEW
C-----WE REQUIRE .333 <=HNEW/W<=6.
      FAC=DMAX1((1.D0/6.D0),DMIN1(3.D0,(ERR/EPSLOC)**(1.D0/8.D0)/.9D0))
      HNEW=H/FAC
      IF(ERR.GT.EPSLOC) GOTO 51
C-----STEP IS ACCEPTED
      NACCPT=NACCPT+1
      DO I=1,N
         Y(I)=K5(I)
      END DO
      X=XPH
      CALL SOLDOPRI(NACCPT+1,X,Y,N)
      IF (DABS(HNEW).GT.HMAXLOC) HNEW=POSNEG*HMAXLOC
      IF (REJECT) HNEW=POSNEG*DMIN1(DABS(HNEW),DABS(H))
      REJECT=.FALSE.
      H=HNEW
      IF (LAST) GO TO 100 
      ! MG 2016/01/05 : evite de memoriser le dernier pas (car en general, il est incomplet)
      HFINAL = H
      GOTO 1
C-----STEP IS REJECT
51    CONTINUE
      REJECT=.TRUE.
**----Modif *m/4 (24/5/95) pour corriger l'erreur sur le dernier pas
**----  Meme si c'etait le dernier, cela ne l'est plus
      LAST=.FALSE.
**----endmodif
      H=HNEW
      IF(NACCPT.GE.1) NREJCT=NREJCT+1
**---- Modif *m/4 (9/5/95) pour prendre en compte le pas XEND-X
      IF ((NSTEP.GT.NMAX).OR.(X+FACSEC*H.EQ.X)) GOTO 79
**---- End Modif
      NFCN=NFCN-1
      GOTO 2
C-----FAIL EXIT
********* a partir d'ici sortie avec erreur.
*           On peut eventuellement
*          recuperer nerror par le common pour traiter 
*          l'erreur dans le programme appelant.
79    CONTINUE
**--- Modif *m/4 (9/5/95) pour recuperer l'erreur
      IF (NSTEP.GT.NMAX) NERRTEMP=1
      IF (X+FACSEC*H.EQ.X) NERRTEMP=2
100   NERROR=NERROR+NERRTEMP
      IF (IPTR.GE.0) THEN
         IF (NERRTEMP.GT.0) THEN
            WRITE(NFPTR,*) ' probleme dans dopri8.'
            WRITE(NFPTR,*) 'SORTIE POUR X=',X
             IF (NERRTEMP.EQ.1) THEN
                WRITE(NFPTR,*) 'TROP D''ETAPES'
             ENDIF
             IF (NERRTEMP.EQ.2) THEN
                WRITE(NFPTR,*) 'LE PAS EST DEVENU TROP PETIT, PAS :',h
             ENDIF
                WRITE(NFPTR,*) 'APPEL : DOPRI8(',N,',',XINPUT,',',       &
     &           XEND,',',EPS,',',HMAX,',',HINPUT,')' 
         ENDIF
         IF (IPTR.GE.1) THEN
            WRITE(NFPTR,*) 'Routine DOPRI8'
            IF (IPTR.GE.2) THEN
               write(NFPTR,*) 'RESULTATS'
               WRITE(NFPTR,*) 'X=',X
               WRITE(NFPTR,*) 'VALEURS DE Y'
               DO I=1,N
                  WRITE(NFPTR,*) Y(I)
               END DO
             ENDIF
         ENDIF
      ENDIF
      ! MG 2016/01/05 : evite de memoriser le dernier pas (car en general, il est incomplet)
      IF (NERRORTEMP.EQ.0.AND.ABS(H/HFINAL).LT.(1.D0/6.D0)) THEN
         H = HFINAL
      ENDIF 
**--- End Modif
      RETURN
      END  SUBROUTINE
C
        SUBROUTINE COEFST
C------ THIS ROUTINE SETS THE COEFFICIENTS FOR THE DORMAND-PRINCE
C------ METHOD OF ORDER 8 WITH ERROR ESTIMATOR OF ORDER 7 AND 13 STAGES
        IMPLICIT REAL(TREAL) (A-H,O-Z)
        COMMON/COEFRK/C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,
     &  A21,A31,A32,A41,A43,A51,A53,A54,A61,A64,A65,A71,A74,A75,A76,
     &  A81,A84,A85,A86,A87,A91,A94,A95,A96,A97,A98,A101,A104,A105,A106,
     &  A107,A108,A109,A111,A114,A115,A116,A117,A118,A119,A1110,A121,
     &  A124,A125,A126,A127,A128,A129,A1210,A1211,A131,A134,A135,A136,
     &  A137,A138,A139,A1310,A1311,B1,B6,B7,B8,B9,B10,B11,B12,B13,
     &  BH1,BH6,BH7,BH8,BH9,BH10,BH11,BH12
        C2=1.D0/18.D0
        C3=1.D0/12.D0
        C4=1.D0/8.D0
        C5=5.D0/16.D0
        C6=3.D0/8.D0
        C7=59.D0/400.D0
        C8=93.D0/200.D0
        C9=5490023248.D0/9719169821.D0
        C10=13.D0/20.D0
        C11=1201146811.D0/1299019798.D0
        C12=1.D0
        C13=1.D0
        A21=C2
        A31=1.D0/48.D0
        A32=1.D0/16.D0
        A41=1.D0/32.D0
        A43=3.D0/32.D0
        A51=5.D0/16.D0
        A53=-75.D0/64.D0
        A54=-A53
        A61=3.D0/80.D0
        A64=3.D0/16.D0
        A65=3.D0/20.D0
        A71=29443841.D0/614563906.D0
        A74=77736538.D0/692538347.D0
        A75=-28693883.D0/1125.D6
        A76=23124283.D0/18.D8
        A81=16016141.D0/946692911.D0
        A84=61564180.D0/158732637.D0
        A85=22789713.D0/633445777.D0
        A86=545815736.D0/2771057229.D0
        A87=-180193667.D0/1043307555.D0
        A91=39632708.D0/573591083.D0
        A94=-433636366.D0/683701615.D0
        A95=-421739975.D0/2616292301.D0
        A96=100302831.D0/723423059.D0
        A97=790204164.D0/839813087.D0
        A98=800635310.D0/3783071287.D0
        A101=246121993.D0/1340847787.D0
        A104=-37695042795.D0/15268766246.D0
        A105=-309121744.D0/1061227803.D0
        A106=-12992083.D0/490766935.D0
        A107=6005943493.D0/2108947869.D0
        A108=393006217.D0/1396673457.D0
        A109=123872331.D0/1001029789.D0
        A111=-1028468189.D0/846180014.D0
        A114=8478235783.D0/508512852.D0
        A115=1311729495.D0/1432422823.D0
        A116=-10304129995.D0/1701304382.D0
        A117=-48777925059.D0/3047939560.D0
        A118=15336726248.D0/1032824649.D0
        A119=-45442868181.D0/3398467696.D0
        A1110=3065993473.D0/597172653.D0
        A121=185892177.D0/718116043.D0
        A124=-3185094517.D0/667107341.D0
        A125=-477755414.D0/1098053517.D0
        A126=-703635378.D0/230739211.D0
        A127=5731566787.D0/1027545527.D0
        A128=5232866602.D0/850066563.D0
        A129=-4093664535.D0/808688257.D0
        A1210=3962137247.D0/1805957418.D0
        A1211=65686358.D0/487910083.D0
        A131=403863854.D0/491063109.D0
        A134=-5068492393.D0/434740067.D0
        A135=-411421997.D0/543043805.D0
        A136=652783627.D0/914296604.D0
        A137=11173962825.D0/925320556.D0
        A138=-13158990841.D0/6184727034.D0
        A139=3936647629.D0/1978049680.D0
        A1310=-160528059.D0/685178525.D0
        A1311=248638103.D0/1413531060.D0
        B1=14005451.D0/335480064.D0
        B6=-59238493.D0/1068277825.D0
        B7=181606767.D0/758867731.D0
        B8=561292985.D0/797845732.D0
        B9=-1041891430.D0/1371343529.D0
        B10=760417239.D0/1151165299.D0
        B11=118820643.D0/751138087.D0
        B12=-528747749.D0/2220607170.D0
        B13=1.D0/4.D0
        BH1=13451932.D0/455176623.D0
        BH6=-808719846.D0/976000145.D0
        BH7=1757004468.D0/5645159321.D0
        BH8=656045339.D0/265891186.D0
        BH9=-3867574721.D0/1518517206.D0
        BH10=465885868.D0/322736535.D0
        BH11=53011238.D0/667516719.D0
        BH12=2.D0/45.D0
        RETURN
        END SUBROUTINE
C
        SUBROUTINE SOLDOPRI(K,X,Y,N)   
        implicit none
         integer, intent(in) :: K,N
         real(TREAL) :: X      
         real(TREAL), dimension(:) :: Y      
        RETURN
        END SUBROUTINE
#endif        
      end module mod_dopri
