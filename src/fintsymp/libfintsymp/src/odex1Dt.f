!***********************************************************************
!* 
!  (c) ASD/IMCCE
!  programmer : Gastineau
! 
!>  \file odex1Dt.f 
!!  \brief version modernisee de odex1D.f mais requiert ode.f
!!
!
! history : creation 05/01/2016
!***********************************************************************
#if TREAL==8
      MODULE CONSTODEX1
      INTEGER,PARAMETER :: NDMAX=66,KMAX=9,NMAX=800
      REAL(8), SAVE :: FAC3=.9D0,FAC4=.8D0
      REAL(8), SAVE :: FAC1=2.D-2,FAC2=4.D0
      REAL(8), SAVE :: SAFE2=.94D0,SAFE1=.65D0
      REAL(8), SAVE :: UROUND=2.26D-16
      REAL(8),SAVE  :: POSNEG
      END MODULE CONSTODEX1
      
      
      
!***********************************************************************
!> module pour l'integration avec odex1 
!***********************************************************************
      module mod_odex1
      use mod_ode

!***********************************************************************
!> @class t_odex1
!! integration avec odex1 :
!! il faut appeller :
!! 1. start : pour le remplissage du vecteur d'etat  et l'initialisation
!de odex1
!! 2. run : pour le remplissage du vecteur d'etat 
!!
!***********************************************************************
      type, extends(t_schemaode) :: t_odex1
          
          real(TREAL), dimension(:), allocatable :: m_y !< vecteur d'etat
          
       contains         
          procedure :: start => odex1_start  !< fonction appellee lors du demarrage
          procedure :: run => odex1_run  !< fonction appellee lors d'un pas normal
          procedure :: check_error => odex1_check_error  ! controle si une erreur s'est produite
      end type t_odex1  
      
      contains
     
!***********************************************************************
! fonction pour demarrer une integration  au temps tstart  avec
! l'integrateur
!***********************************************************************
       subroutine odex1_start(this, system, tstart, dt)
       implicit none
       class(t_odex1), intent(inout) :: this
       class(t_sysode), intent(inout) :: system ! systeme qui est integre
       real(TREAL), intent(in) :: dt ! pas d'integration
       real(TREAL), intent(in) :: tstart ! temps initial
        COMMON/ASD/IPTR,NFPTR,NERROR
        integer IPTR,NFPTR,NERROR
       
       ! remplissage du vecteur d'etat
       call system%allocate_y(this%m_y)
       call system%fill_to_y(this%m_y, tstart)
       ! initialise les variables de odex1
       IPTR = 0
       NFPTR = 6
       NERROR = 0
       
       end subroutine odex1_start

!***********************************************************************
! fonction pour executer un pas normal d'integration avec odex1
!***********************************************************************
      subroutine odex1_run(this, system, tstart, tend, dt)
      implicit none
      class(t_odex1), intent(inout) :: this
      class(t_sysode), intent(inout) :: system ! systeme qui est integre
      real(TREAL), intent(in) :: dt ! pas d'integration
      real(TREAL), intent(in) :: tstart ! temps initial
      real(TREAL), intent(in) :: tend ! temps final = tstart+dt
      
      call ODEX1D(size(this%m_y),system, tstart, this%m_y,tend,         &
     &    this%m_eps, this%m_hmax, this%m_h)
         call system%fill_from_y(this%m_y)
      
      end subroutine odex1_run

!***********************************************************************
!> @brief controle si une erreur s'est produite. 
!! retourne true en cas d'erreur
!***********************************************************************
       function odex1_check_error(this) result(ret)
        implicit none
        class(t_odex1), intent(inout) :: this  !< dummy argument
        logical ret ! code de retour
        COMMON/ASD/IPTR,NFPTR,NERROR
        integer IPTR,NFPTR,NERROR
        ret = .false.
        if (NERROR.ne.0) then
          ret = .true.        
        endif         
       end function odex1_check_error

      SUBROUTINE ODEX1D (N, system, X, Y, XEND, EPS,HMAX,H)
! ----------------------------------------------------------
!     NUMERICAL SOLUTION OF A SYSTEM OF FIRST ORDER 
!     DIFFERENTIAL EQUATIONS Y'=F(X,Y)
!     THIS IS AN EXTRAPOLATION-ALGORITHM,BASED ON THE
!     EXPLICIT MIDPOINT RULE (WITH STEPSIZE CONTROL
!     AND ORDER SELECTION).
!     C.F. SECTION II.9
!     
!     INPUT PARAMETERS
!     ----------------
!     N           DIMENSION OF THE SYSTEM (N.LE.66)
!     system%FCN  NAME (EXTERNAL ) OF SUBROUTINE COMPUTING THE
!                 FIRST DERIVATIVE F(X,Y)
!                 SUBROUTINE FCN(N,X,Y,F)
!                 REAL*8 X,Y(N),F(N)
!     X           INITIAL X-VALUES
!     XEND        FINAL X-VALUES (XEND.GT.X)
!     Y(N)        INITIAL VALUES
!     EPS         LOCAL TOLERANCE
!     HMAX        MAXIMAL STEPSIZE
!     H           INITIAL STEPSIZE GUESS
!
!     OUTPUT PARAMETERS  
!     -----------------
!     Y(N)        SOLUTION AT XEND
!     
!     EXTERNAL SUBROUTINE (TO BE SUPPLIED BY THE USER)
!     -------------------
!     SOLOUT      THIS SUBROUTINE IS COMPUTED AFTER EVERY
!                 SUCCESSFULLY COMPUTED STEP (AND THE 
!                 INITIAL VALUE)
!                       SUBROUTINE SOLOUT(NR,X,Y,N)
!                       REAL*8 X,Y(N)
!                 FURNISHES THE SOLUTION Y AT THE NR-TH
!                 GRID-POINT X (THE INITIAL VALUE IS CON-
!                 SIDERED AS THE FIRST GRID-POINT).
!                 SUPPLIED A DUMMY SUBROUTINE, IF THE SOLUTION
!                 IS NOT DESIRED AT THE INTERMEDIATE POINTS.
!     
!------------------------------------------------------------
!  VERSION 1.02   (*M/4) ASD 
!   3/9/98
!   PRISE EN COMPTE DE L'ERREUR RELATIVE 3/9/98 
!   51 -> NDMAX, 9-> KMAX, KM=KMAX
!  7/1/00 CORRECTION DU RISQUE DE DEPASSEMEMT DE XEND LORS
!         D'UN REJET D'UN PAS
!  4/2/00 CORRECTION DES APPELS MIDEX
!------------------------------------------------------------
      USE CONSTODEX1
      IMPLICIT REAL(8) (A-H,O-Z)
! --- NMAX     MAXIMAL NUMBER OF STEPS
!     NDMAX DIMENSION MAXIMALE DU SYSTEME
!     KMAX ORDRE MAXIMAL DE LA METHODE 
      LOGICAL REJECT, LAST
      LOGICAL,SAVE :: FIRST=.TRUE.
      class(t_sysode), intent(inout) :: system ! systeme qui est integre
      REAL(8) Y(N)
      COMMON /ODEXST/NFCN,NSTEP,NACCPT,NREJCT,NORDRE
      COMMON /ASD/IPTR,NFPTR,NERROR
C --- COMMON ODSTAT CAN BE USED FOR STATISTICS
C ---    NFCN       NUMBER OF FUNCTION EVALUATIONS
C ---    NSTEP      NUMBER OF COMPUTED STEPS
C ---    NACCPT     NUMBER OF ACCEPTED STEPS
C ---    NREJCT     NUMBER OF REJECTED STEPS
      COMMON /EXTABL/ DZ(NDMAX),T(KMAX,NDMAX),HH(KMAX),
     &  W(KMAX),ERR,FAC,EPSD4
      COMMON /EXTABL1/ A(KMAX),NJ(KMAX)
C --- INITIAL PREPARATIONS
      IF (FIRST) THEN
        DO I=1,KMAX
          NJ(I) = 2*I
          A(I) = 1.D0+ I*(I+1)*1.D0
        END DO
        FIRST=.FALSE.
      ENDIF
      KM=KMAX
      EPSD4=EPS*SAFE1
      NSTEP=0
      NREJCT=0
      NACCPT=0
      NFCN=0
      NERRTEMP=0
      K=MAX0(3,MIN0(8,INT(-DLOG10(EPS)*.6D0+1.5D0)))
      POSNEG = DSIGN(1.D0,XEND-X)
      H=DMIN1(DABS(H),DABS(HMAX),DABS((XEND-X))/2.D0)
      H=DSIGN(H,POSNEG)
      CALL SOLODEX (NACCPT+1,X,Y,N,H)
      ERR =0.D0
      W(1)=0.D0
      REJECT=.FALSE.
      LAST=.FALSE.
C --- IS XEND REACHED IN THE NEXT STEP ?
 10   CONTINUE
      H1=DABS(XEND-X)
      NORDRE=K
! MODIF 3/9/98
!     IF ( H1.LE.UROUND) GOTO 110
      IF (H1.LE.UROUND*DMAX1(1.D0,DABS(XEND),DABS(X))) GOTO 110
      H=DMIN1(DABS(H),H1,DABS(HMAX))
! MODIF 3/9/98
!      IF (H.GE.H1-UROUND) LAST=.TRUE.
       IF (DABS(H1-H) .LE. UROUND*DMAX1(1.D0,H1,H)) LAST=.TRUE.
       H=DSIGN(H,POSNEG)
      CALL system%FCN(N,X,Y,DZ)
      NFCN=NFCN+1
C --- THE FIRST AND LAST STEP
      IF (NSTEP.EQ.0.OR.LAST) THEN
         NSTEP = NSTEP+1
         DO  J=1,K
           KC=J
           CALL MIDEX(J,X,Y,H,HMAX,N,system)
         IF (J.GT.1.AND.ERR.LE.EPS) GOTO 60
         ENDDO
         GOTO 55
      END IF
C --- BASIC INTEGRATION STEP
 30   CONTINUE
      NSTEP=NSTEP+1
      IF (NSTEP.GE.NMAX) GOTO 120
      KC=K-1
      DO J=1,KC
        CALL MIDEX(J,X,Y,H,HMAX,N,system)
      ENDDO
C --- CONVERGENCE MONITOR
      IF (K.EQ.2.OR.REJECT) GOTO 50
      IF (ERR.LE.EPS) GOTO 60
      IF (ERR/EPS.GT.(DFLOAT(NJ(K+1)*NJ(K))/4.D0)**2) GOTO 100
! correction 4/2/00      
50    CALL MIDEX(K,X,Y,H,HMAX,N,system)
! 50   CALL MIDEX(j,X,Y,H,HMAX,N,FCN)
! fin correction 4/2/00
      KC=K
      IF (ERR.LE.EPS) GOTO 60
C --- HOPE FOR CONVERGENCE IN LINE K+1
 55   IF (ERR/EPS.GT.(DFLOAT(NJ(K+1))/2.D0)**2) GOTO 100
      KC = K+1
! correction 4/2/00
      CALL MIDEX(KC,X,Y,H,HMAX,N,system)
!      CALL MIDEX(J,X,Y,H,HMAX,N,FCN)
! fin correction 4/2/00
      IF (ERR.GT.EPS) GOTO 100
C --- STEP IS ACCEPTED
 60   X=X+H
      DO I =1,N
        Y(I)=T(1,I)
      ENDDO
      NACCPT =NACCPT+1
      CALL SOLODEX (NACCPT+1,X,Y,N,H)
C --- COMPUTE OPTIMAL ORDER
      IF (KC.EQ.2) THEN
         KOPT=3
         IF (REJECT) KOPT=2
         GOTO  80
      END IF
      IF (KC.LE.K) THEN
         KOPT=KC
         IF (W(KC-1).LT.W(KC)*FAC3) KOPT=KC-1
         IF (W(KC).LT.W(KC-1)*FAC3) KOPT=MIN0(KC+1,KM-1)
      ELSE
         KOPT = KC-1
         IF (KC.GT.3.AND.W(KC-2).LT.W(KC-1)*FAC3) KOPT=KC-2
         IF (W(KC).LT.W(KOPT)*FAC3) KOPT=MIN0(KC,KM-1)
      END IF
C --- AFTER A REJECTED STEP
 80    CONTINUE
       IF (REJECT) THEN
         K=MIN0(KOPT,KC)
         H=DMIN1(DABS(H),DABS(HH(K)))
         H=DSIGN(H,POSNEG)
         REJECT=.FALSE.
         GOTO 10
      END IF
C --- COMPUTER STEPSIZE FOR NEXT STEP
      IF (KOPT.LE.KC) THEN
         H=HH(KOPT)
      ELSE
         IF (KC.LT.K.AND.W(KC).LT.W(KC-1)*FAC4) THEN
            H=HH(KC)*A(KOPT+1)/A(KC)
         ELSE
            H=HH(KC)*A(KOPT)/A(KC)
         END IF
      END IF
      K=KOPT
      GOTO 10
C --- STEP IS REJECTED
 100  K=MIN0(K,KC)
      IF (K.GT.2.AND.W(K-1).LT.W(K)*FAC3) K=K-1
      NREJCT=NREJCT+1
! MODIF 7/1/00
! POUR EVITER LE CAS D'UN REJET ALORS
! QU'UNE LIGNE DU TABLEAU D'INTERPOLATION DONNE UN
! RESULTAT ACCEPTABLE ET QUE L'ON SE SITUE PRES DE
! XEND, ON FORCE LE CALCUL AVEC UNE NOUVELLE VALEUR
! DE K CORRESPONDANT A LA VALEUR EVENTUELLE DONNANT
! UN AGGRANDISSEMENT DU PAS    
      H=DMIN1(DABS(HH(K)),H1)
      H=DSIGN(H,POSNEG)
!      H=HH(K)
! FIN MODIF 7/1/00
      REJECT=.TRUE.
      GOTO 30
C ---  SORTIE
C******  ECHEC
120   NERRTEMP=1
      IF (IPTR.GE.0) THEN
         WRITE(NFPTR,*) 'ERREUR DANS ODEX1'
         WRITE(NFPTR,*) 'TROP D''ETAPES'
      ENDIF
      NERROR=NERROR+NERRTEMP
C*****  REUSSITE
110   CONTINUE
      IF (IPTR.GE.1) THEN
         WRITE(NFPTR,*) ' ROUTINE ODEX1'
         IF (IPTR.GE.2) THEN
            WRITE(NFPTR,*) ' X =',X
            WRITE(NFPTR,*) 'VALEURS DE Y:'
            DO I=1,N
               WRITE(NFPTR,*) I,' : ',Y(I)
            END DO
          ENDIF
       ENDIF           
      RETURN
      END SUBROUTINE ODEX1D 
C
      SUBROUTINE MIDEX(J,X,Y,H,HMAX,N,system)
C --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
C --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATION
C --- OF THE OPTIMAL STEPSIZE
      USE CONSTODEX1
      IMPLICIT REAL(8) (A-H,O-Z)
      class(t_sysode), intent(inout) :: system ! systeme qui est integre
      DIMENSION Y(N),DY(NDMAX),YH1(NDMAX),YH2(NDMAX)
      COMMON /ODEXST/NFCN,NSTEP,NACCPT,NREJCT,NORDRE
      COMMON /EXTABL/ DZ(NDMAX),T(KMAX,NDMAX),HH(KMAX),
     & W(KMAX),ERR,FAC,EPSD4
      COMMON /EXTABL1/ A(KMAX),NJ(KMAX)
      HJ=H/DFLOAT(NJ(J))
C --- EULER STARTING STEP
      DO I=1,N
        YH1(I)=Y(I)
        YH2(I)=Y(I)+HJ*DZ(I)
      END DO
C --- EXPLICIT MIDPOINT RULE
      M=NJ(J)-1
      DO  MM=1,M
        CALL system%FCN(N,X+HJ*DFLOAT(MM),YH2,DY)
        DO I=1,N
          YS=YH1(I)
          YH1(I)=YH2(I)
          YH2(I)=YS+2.D0*HJ*DY(I)
        END DO
      END DO
C --- FINAL SMOOTHING STEP
      CALL system%FCN(N,X+H,YH2,DY)
      DO I=1,N
        T(J,I)=(YH1(I)+YH2(I)+HJ*DY(I))/2.D0
      END DO
      NFCN=NFCN+NJ(J)
C --- POLYNOMIAL EXTRAPOLATION
      IF (J.EQ.1) RETURN
      DO L=J,2,-1
        FAC=(DFLOAT(NJ(J))/DFLOAT(NJ(L-1)))**2-1.D0
        DO  I=1,N
          T(L-1,I)=T(L,I)+(T(L,I)-T(L-1,I))/FAC
        END DO
      END DO
      ERR = 0.D0
      DO I=1,N
C --- SCALING
        SCAL=DMAX1(DABS(Y(I)),DABS(T(1,I)),1.D-6,UROUND/EPSD4)
        ERR=ERR+((T(1,I)-T(2,I))/SCAL)**2
      END DO
      ERR=DSQRT(ERR/DFLOAT(N))
C --- COMPUTE OPTIMAL STEPSIZES
      EXPO=1.D0/DFLOAT(2*J-1)
      FACMIN=FAC1**EXPO
      FAC=DMIN1(FAC2/FACMIN,DMAX1(FACMIN,(ERR/EPSD4)**EXPO/SAFE2))
      FAC=1.D0/FAC
      HH(J)=DMIN1(DABS(H*FAC),DABS(HMAX))
      W(J)=A(J)/HH(J)
      HH(J)=DSIGN(HH(J),POSNEG)
      RETURN
      END SUBROUTINE MIDEX



********************************
*  PARTIE UTILISABLE PAR L'UTILISATEUR
***********************************
C
      SUBROUTINE SOLODEX (NRPNTS,X,Y,N,H)
      IMPLICIT REAL(8) (A-H,O-Z)
      COMMON /ODEXST/NFCN,NSTEP,NACCPT,NREJCT,NORDRE
      COMMON /ODEXER/NERROR
      COMMON/TOT/PASMIN,NFTOT,NRTOT
      DIMENSION Y(N)
      RETURN
      END SUBROUTINE SOLODEX
      
      end module mod_odex1
#endif
      
