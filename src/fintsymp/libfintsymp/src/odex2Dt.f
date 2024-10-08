!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file odex2Dt.f 
!!  \brief version modernisee de odex2D.f mais requiert ode.f
!!   integrateur du second ordre
!!
!
! history : creation 05/01/2016
!***********************************************************************
#include "realc.f"

#if TREAL==8
!***********************************************************************
!> module pour l'integration avec odex2
!***********************************************************************
      module mod_odex2
      use mod_ode

!***********************************************************************
!> @class t_odex2
!! integration avec odex2 :
!! il faut appeller :
!! 1. start : pour le remplissage du vecteur d'etat  et l'initialisation
!de odex2
!! 2. run : pour le remplissage du vecteur d'etat 
!!
!***********************************************************************
      type, extends(t_schemaode) :: t_odex2
          
          real(TREAL), dimension(:), allocatable :: m_y !< vecteur d'etat
          
       contains         
          procedure :: start => odex2_start  !< fonction appellee lors du demarrage
          procedure :: run => odex2_run  !< fonction appellee lors d'un pas normal
          procedure :: check_error => odex2_check_error  ! controle si une erreur s'est produite
      end type t_odex2  
      
      contains
     
!***********************************************************************
! fonction pour demarrer une integration  au temps tstart  avec
! l'integrateur
!***********************************************************************
       subroutine odex2_start(this, system, tstart, dt)
       implicit none
       class(t_odex2), intent(inout) :: this
       class(t_sysode), intent(inout) :: system ! systeme qui est integre
       real(TREAL), intent(in) :: dt ! pas d'integration
       real(TREAL), intent(in) :: tstart ! temps initial
!        COMMON/ASD/IPTR,NFPTR,NERROR
!        integer IPTR,NFPTR,NERROR
       
       ! remplissage du vecteur d'etat
       call system%allocate_y(this%m_y)
       call system%fill_to_y(this%m_y, tstart)
       ! initialise les variables de odex1
!       IPTR = 0
!       NFPTR = 6
!       NERROR = 0
       
       end subroutine odex2_start

!***********************************************************************
! fonction pour executer un pas normal d'integration avec odex2
!***********************************************************************
      subroutine odex2_run(this, system, tstart, tend, dt)
      implicit none
      class(t_odex2), intent(inout) :: this
      class(t_sysode), intent(inout) :: system ! systeme qui est integre
      real(TREAL), intent(in) :: dt ! pas d'integration
      real(TREAL), intent(in) :: tstart ! temps initial
      real(TREAL), intent(in) :: tend ! temps final = tstart+dt
      integer n
      n = size(this%m_y)
      call ODEX2D(n/2,system, tstart, this%m_y, this%m_y(n/2+1), tend,         &
     &    this%m_eps, this%m_hmax, this%m_h)
         call system%fill_from_y(this%m_y)
      
      end subroutine odex2_run

!***********************************************************************
!> @brief controle si une erreur s'est produite. 
!! retourne true en cas d'erreur
!***********************************************************************
       function odex2_check_error(this) result(ret)
        implicit none
        class(t_odex2), intent(inout) :: this  !< dummy argument
        logical ret ! code de retour
        !COMMON/ASD/IPTR,NFPTR,NERROR
        !integer IPTR,NFPTR,NERROR
        ret = .false.
        ! if (NERROR.ne.0) then
        !  ret = .true.        
        !endif         
       end function odex2_check_error

       SUBROUTINE ODEX2D(N,system,X,Y,YP,XEND,EPS,HMAX,H)
C-----------------------------------------------------------
C      NUMERICAL SOLUTION OF A SYSTEM OF SECOND ORDER
C      ORDINARY DIFFERENTIAL EQUATIONS Y"=F(X,Y).
C      THIS IS AN EXTRAPOLATION-ALGORITHM, BASED ON THE
C      EXPLICIT MIDPOINT RULE (WITH STEP SIZE CONTROL
C      AND ORDER SELECTION).
C      C.F. SECTION II.13
C-- HAIRER, NORSET & WANNER p433
C    
C      INPUT PARAMETERS
C      ----------------
C      N             DIMENSION OF THE SYSTEM (N.LE.NDMAX)
C      FCN           NAME (EXTERNAL) OF SUBROUTINE COMPUTING
C                    THE SECOND DERIVATIVE F(X,Y)
C                      SUBROUTINE FCN(N,X,Y,F)
C                      REAL*8 X,Y(N),F(N)
C                      F(1)= ... ETC.
C      X             INITIAL X-VALUE
C      XEND          FINAL X-VALUE (XEND .GT. X)
C      Y(N)          INITIAL VALUED FOR Y
C      YP(N)         INITIAL VALUES FOR Y'
C      EPS           LOCAL TOLERANCE
C      HMAX          MAXIMAL STEPSIZE
C      H             INITIAL STEPSIZE GUESS
C
C      OUTPUT PARAMETERS
C      -----------------
C      Y(N)          SOLUTION AT XEND
C      YP(N)         DERIVATIVE OF SOLUTION AT XEND
C
C      EXTERNAL SUBROUTINE (TO BE SUPPLIED BY THE USER)
C      ------------------------------------------------
C      SOLOUT2       THIS SUBROUTINE IS CALLED AFTER
C                    EVERY SUCCESSFULLY COMPUTED STEP
C                    (AND THE INITIAL VALUE):
C                      SUBROUTINE SOLOUT2(NR,X,Y,YP,N)
C                      REAL*8 X,Y(N),YP(N)
C                    FURNISHES THE SOLUTION (Y,YP) AT THE
C                    NR-TH GRID-POINT X (THE INITIAL VALUE
C                    IS CONSIDERED AS THE FIRST GRID-POINT).
C
C                    SUPPLY A DUMMY ROUTINE IF THE SOLUTION IS 
C                    NOT DESIRED AT THE INTERMEDIATE POINTS.
C--------------------------------------------------------------- 
C modif 11/12/98 51 -> NDMAX
C----------------------------------------------------------------
      IMPLICIT double precision (A-H,O-Z)
      PARAMETER(NDMAX=54) 
      LOGICAL REJECT,LAST
      dimension Y(N),YP(N)
      class(t_sysode), intent(inout) :: system ! systeme qui est integre
      COMMON /STATODEX/ NFCN,NSTEP,NACCPT,NREJCT
C---- COMMON STATODEX CAN BE USED FOR STATISTICS
C---       NFCN      NUMBER OF FUNCTIONS EVALATION
C---       NSTEP     NUMBER OF COMPUTED STEPS
C---       NACCPT    NUMBER OF ACCEPTED STEPS
C---       NREJCT    NUMBER OF REJECTED STEPS
      COMMON /EXTABL2/ DZ(NDMAX),T(9,NDMAX),TP(9,NDMAX),
     1     A(9),HH(9),W(9),NJ(9)
      COMMON/EXTVAR/ERR,FAC,EPSD4,UROUND,FAC1,FAC2,SAFE2,posneg
      DATA NJ/2,4,6,8,10,12,14,16,18/
      DATA A/2.D0,4.D0,7.D0,11.D0,16.D0,22.D0,29.D0,37.D0,46.D0/
      DATA NMAX/80000/,KM/9/,UROUND/2.26D-32/
C---       NMAX      MAXIMAL NUMBER OF STEPS
C---       UROUND    SMALLEST NUMBER SATISFYING 1.D0+UROUND > 1.D0
C---                 (TO BE ADAPTED BY THE USER)
      DATA FAC1/2.D-2/, FAC2/4.D0/, FAC3/.9D0/, FAC4/.8D0/
      DATA SAFE1/.65D0/, SAFE2/.94D0/
C--- INITIAL PREPARATIONS
      IF (N.GT.NDMAX) THEN
        Write(*,*) ' Dans ODEX, dimension trop grande'
        stop
      ENDIF
      NMIN=3
      EPSD4=EPS*SAFE1
      NSTEP=0
      NREJCT=0
      NACCPT=0
      NFCN=0
      K=MAX0(NMIN,MIN0(8,INT(-LOG10(EPS)*.6D0+1.5D0)))
      posneg = sign(1.D0,xend-x)
      H=MIN(abs(H),abs(HMAX),abs(XEND-X)/2.D0)
      h=sign(h,posneg)
      CALL SOLOUT2(NACCPT+1,X,Y,YP,N,k,h)
      ERR=0.D0
      W(1)=0.D0
      REJECT=.FALSE.
      LAST=.FALSE.
C--- IS XEND REACHED IN THE NEXT STEP ?
10    H1=abs(XEND-X)
      IF (H1.LE.UROUND) GOTO 110
      H=MIN(abs(H),H1,abs(HMAX))
      H=sign(h,posneg)
C     IF (H.GE.H1-UROUND) LAST=.TRUE.
      IF((X+H-XEND)*POSNEG.GE.0.D0) then
         H=XEND-X
         LAST=.TRUE.
      ENDIF
      CALL system%FCN(N,X,Y,DZ)
      NFCN=NFCN+1
C--- THE FIRST AND THE LAST STEP
      IF (NSTEP.EQ.0 .OR. LAST) THEN
         NSTEP=NSTEP+1
         DO 20 J=1,K
         KC=J
         CALL STOERM(J,X,Y,YP,H,HMAX,N,system)
20       IF (J.GT.1 .AND. ERR.LE.EPS) GOTO 60
         GOTO 55
      ENDIF
C--- BASIC INTEGRATION STEP
30    CONTINUE
      NSTEP=NSTEP+1
      IF (NSTEP.GE.NMAX) GOTO 120
      KC=K-1
      DO 40 J=1,KC
40    CALL STOERM(J,X,Y,YP,H,HMAX,N,system)
C--- CONVERGENCE MONITOR
      IF (K.EQ.NMIN-1 .OR. REJECT) GO TO 50
      IF (ERR.LE.EPS) GO TO 60
      IF (ERR/EPS .GT. (FLOAT(NJ(K+1)*NJ(K))/4.D0)**2) GO TO 100
50    CALL STOERM(K,X,Y,YP,H,HMAX,N,system)
      KC=K
      IF (ERR.LE.EPS) GO TO 60
C--- HOPE IN CONVERGENCE IN LINE K+1
55    IF (ERR/EPS .GT. (FLOAT(NJ(K+1))/2.D0)**2) GO TO 100
      KC=K+1
      CALL STOERM(KC,X,Y,YP,H,HMAX,N,system)
      IF (ERR.GT.EPS) GO TO 100
C--- STEP IS ACCEPTED
60    X=X+H
      DO 70 I=1,N
      YP(I)=TP(1,I)
70    Y(I)=T(1,I)
      NACCPT=NACCPT+1
      CALL SOLOUT2(NACCPT+1,X,Y,YP,N,K,h)
C--- COMPUTE OPTIMAL ORDER
      IF (KC.EQ.NMIN-1) THEN
         KOPT=NMIN
         IF (REJECT) KOPT=NMIN-1
         GO TO 80
      ENDIF
      IF (KC.LE.K) THEN
         KOPT=KC
         IF (W(KC-1).LT.W(KC)*FAC3) KOPT=KC-1
         IF (W(KC).LT.W(KC-1)*FAC3) KOPT=MIN0(KC+1,KM-1)
      ELSE
         KOPT=KC-1
         IF (KC.GT.NMIN .AND. W(KC-2).LT.W(KC-1)*FAC3) KOPT=KC-2
         IF (W(KC).LT.W(KOPT)*FAC3) KOPT=MIN0(KC,KM-1)
      ENDIF
C--- AFTER A REJECTED STEP
80    IF (REJECT) THEN
         K=MIN0(KOPT,KC)
         H=MIN(abs(H),abs(HH(K)))
         H=sign(H,posneg)
         REJECT=.FALSE.
         GO TO 10
      ENDIF
C--- COMPUTE STEPSIZE FOR NEXT STEP
      IF (KOPT.LE.KC) THEN
         H=HH(KOPT)
      ELSE
         IF (KC.LT.K .AND. W(KC).LT.W(KC-1)*FAC4) THEN
            H=HH(KC)*A(KOPT+1)/A(KC)
         ELSE
            H=HH(KC)*A(KOPT)/A(KC)
         ENDIF
      ENDIF
      K=KOPT
      GO TO 10
C--- STEP IS REJECTED
100   K=MIN0(K,KC)
      IF (K.GT.NMIN-1 .AND. W(K-1).LT.W(K)*FAC3) K=K-1
      NREJCT=NREJCT+1
      H=HH(K)
      REJECT=.TRUE.
      GO TO 30
C--- SOLUTION EXIT
110   CONTINUE
      RETURN
C--- FAIL EXIT
120   WRITE(6,*) ' MORE THAN  ',NMAX,' STEPS '
      WRITE(6,*) ' EXIT OF ODEX2 AT X = ',X
      RETURN
      END SUBROUTINE
C
      SUBROUTINE STOERM(J,X,Y,YP,H,HMAX,N,system)
C--- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
C--- EXTRAPOLATION TABLE AND PROVIDES AN ESTIAMTION
C--- OF THE OPTIMAL STEPSIZE
      IMPLICIT double precision (A-H,O-Z)
      PARAMETER(NDMAX=54)
      class(t_sysode), intent(inout) :: system ! systeme qui est integre
      dimension Y(N),YP(N),DY(NDMAX),YH1(NDMAX),ZH1(NDMAX) 
      COMMON/STATODEX/NFCN,NSTEP,NACCPT,NREJCT
      COMMON /EXTABL2/ DZ(NDMAX),T(9,NDMAX),TP(9,NDMAX),
     1     A(9),HH(9),W(9),NJ(9)
      COMMON/EXTVAR/ERR,FAC,EPSD4,UROUND,FAC1,FAC2,SAFE2,posneg
      HJ=H/FLOAT(NJ(J))
      HJ2=HJ*2.D0
C--- EULER STARTING STEP
      DO 30 I=1,N
      YH1(I)=Y(I)
30    ZH1(I)=YP(I)+HJ*DZ(I)
C--- EXPLICIT MIDPOINT (STOERMER) RULE
      M=NJ(J)/2-1
      IF (J.EQ.1) GOTO 37
      DO 35 MM=1,M
      DO 33 I=1,N
33    YH1(I)=YH1(I)+HJ2*ZH1(I)
      CALL system%FCN(N,X+HJ2*FLOAT(MM),YH1,DY)
      DO 35 I=1,N
35    ZH1(I)=ZH1(I)+HJ2*DY(I)
C--- FINAL STEP
37    CONTINUE
      DO 40 I=1,N
40    YH1(I)=YH1(I)+HJ2*ZH1(I)
      CALL system%FCN(N,X+H,YH1,DY)
      DO 43 I=1,N
      T(J,I)=YH1(I)
43    TP(J,I)=ZH1(I)+HJ*DY(I)
      NFCN=NFCN+M+1
C--- POLYNOMIAL EXTRAPOLATION
      IF (J.EQ.1) RETURN
      DO 60 L=J,2,-1
      FAC=(FLOAT(NJ(J))/FLOAT(NJ(L-1)))**2-1.D0
      DO 60 I=1,N
      T(L-1,I)=T(L,I)+(T(L,I)-T(L-1,I))/FAC
      TP(L-1,I)=TP(L,I)+(TP(L,I)-TP(L-1,I))/FAC
60    CONTINUE
      ERR=0.D0
      DO 65 I=1,N
C--- SCALING
      SCAL=MAX(ABS(Y(I)),ABS(T(1,I)),1.D-6,UROUND/EPSD4)
      SCALP=MAX(ABS(YP(I)),ABS(TP(1,I)),1.D-6,UROUND/EPSD4)
65    ERR=ERR+((T(1,I)-T(2,I))/SCAL)**2+((TP(1,I)-TP(2,I))/SCALP)**2  
      ERR=SQRT(ERR/FLOAT(N*2))
C--- COMPUTE OPTIMAL STEPSIZES
      EXPO=1.D0/FLOAT(2*J-1)
      FACMIN=FAC1**EXPO
      FAC=MIN(FAC2/FACMIN,MAX(FACMIN,(ERR/EPSD4)**EXPO/SAFE2))
      FAC=1.D0/FAC
      HH(J)=MIN(abs(H*FAC),abs(HMAX))
      W(J)=A(J)/HH(J)
      HH(J)=sign(HH(J),posneg)
      RETURN
      END SUBROUTINE
C
      SUBROUTINE SOLOUT2(NRPNTS,X,Y,YP,N,K,h)
      IMPLICIT double precision (A-H,O-Z)
      dimension Y(N),YP(N)
      RETURN
      END SUBROUTINE

      end module mod_odex2
#endif

          

