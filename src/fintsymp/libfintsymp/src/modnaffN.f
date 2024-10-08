!-------------------------------------------------------------------------
!   NAFF.0.84 NUMERICAL ANALYSIS OF FUNDAMENTAL FREQUENCIES
!            (27 JANVIER 1996)
!  
!   (C) JACQUES LASKAR 
!       ASTRONOMIE ET SYSTEMES DYNAMIQUES
!       BUREAU DES LONGITUDES
!       75006 PARIS
!       EMAIL : LASKAR@BDL.FR
!
!    MAIN REFERENCES : 
!    LASKAR, J.: THE CHAOTIC MOTION OF THE SOLAR SYSTEM. A NUMERICAL
!    ESTIMATE OF THE SIZE OF THE CHAOTIC ZONES,ICARUS,88,(1990),266--291
!
!   (NAFFMOD - 07/09/01)
!
!**************************************************************************
!  THIS PROGRAMM  CANNOT BE COPIED, 
!  DISTRIBUTED NOR MODIFIED WITHOUT THE AGREEMENT OF THE AUTHOR.
!**************************************************************************
! 
!                        PROGRAMME D'ANALYSE DE FOURIER AUTOMATIQUE
!   -- NAFF --           CALCULE UN TERME A LA FOIS ET LE RETRANCHE
!                        A TABS
!                        NOUVELLE VERSION 26/9/87
!                        MODIFIEE POUR NFS LE 18/10/87
!                        MODIFIEE LE 10/4/88 POUR TRAITER
!                        A LA FOIS UN TABLEAU REEL OU COMPLEXE
!                        CALCUL BATCH SUR CIRCE
! 8/10/90 MODIF POUR REMETTRE SUR VP (FORTRAN STANDARD)
! 13/12/90   MODIF POUR RENDRE PLUS MODULAIRE
! 18/12/90   MODIF IL N'Y A PLUS DE DIMENSIONS A REGLER DANS
!            LES SOUS PROGRAMMES. TABS(2,*) AU LIEU DE TABS(2,0:KTABS)
!            NFR EST ENCORE EN PARAMETER, MAIS EST RAREMENT A CHANGER.
!
!            IL FAUT REMPLACER LE * DES TABLEAUX
!            DIMENSIONES T(*) OU T(N,*) PAR T(1) OU T(N,1)
! 19/12/90
!  16/4/91   COMPILATION SEPAREE DE NAFF EN UNE BIBLIOTHEQUE
!  30/4/91   TRAITEMENT DU CAS OU ON RETROUVE LA MEME FREQUENCE
!  27/1/96   PLUSIEURS FENETRES
!  27/2/97   MODIF POUR DIMENSIONNER TOUT EXPLICITEMENT AVEC 
!            INCLUDE NAFF.INC POUR LES PARAMETRES
!  23/45/1997    CORRECTION DE MAXIQUA.F (F. JOUTEL)
!  22/4/1998 MODIF POUR QUADRUPLE PRECISION (J. LASKAR)
!            CHANGEMENT POUR UN EPSILON MACHINE DONNE PAR NAFF.INC
!  27/5/98   AMELIORATION DE LA RECHERCHE DU MAXIMUM PAR
!            LA METHODE DES SECANTES
!            ROUTINE DE CORRECTION DE LA FREQUENCE PRINCIPALE PAR
!            DEVELOPPEMENT ASYMTOTIQUE (F. JOUTEL)
!  31/5/00   Version 1.91
!  07/09/01  TOUTES LES VARIABLES PUBLIQUES SONT RENOMMEES EN NAF_xxx
!            AJOUT DE TOL DANS LA PARTIE PUBLIQUE (=>NAF_TOL)
!            Version 2.00
!  24/02/03  ...
!***********************************************************************
!          MODULE  NAFF
!----------------------------------------------------------------------
!  PROCEDURE D'UTILISATION:
!  - PHASE D'INITIALISATION
!      1)  INITIALISER LES VARIABLES :
!              NAF_DTOUR,NAF_KTABS,NAF_XH,NAF_NTERM,NAF_IW,NAF_T0,
!              NAF_ICPLX,NAF_ISEC,NAF_NFPRT,NAF_IPRT,NAF_TOL
!      2)  CALL INITNAF
!             CALCULE NAF_UNIANG,NAF_FREFON, NAF_EPSM, NAF_PI 
!             CREE LES TABLEAUX DYNAMIQUES ((NAF_ZTABS)), ((NAF_TFS)), 
!              ((NAF_ZAMP)), ((NAF_ZALP)) et((TWIN))
!             APPELLE INIWIN
!  - PHASE D'UTILISATION
!      1)  REMPLIR LE TABLEAU NAF_ZTABS
!      2)  CALL MFTNAF(NBTERM,EPS)
!          ANALYSE ((NAF_NAF_ZTABS)) ET ESSAIE DE RECHERCHER (NBTERM)  
!          FREQUENCES, (EPS) ETANT L'ERREUR TOLEREE
!          --! MODIFIE ((NAF_ZTABS))
!  - PHASE D'EXPLOITATION OU DE MODIFICATION
!     * CALL PRTABS(KTABS,ZTABS,IPAS)
!         IMPRIME LES ELEMENTS DU TABLEAU COMPLEXE (ZTABS(0:KTABS)) TOUS 
!         LES (IPAS), ((NAF_IPRT)) DOIT ETRE MIS A 1 AVANT LA PROCEDURE
!         POUR IMPRIMER QUELQUE CHOSE 
!     * CALL INIFRE
!         REMET A ZERO ((NAF_TFS)), ((NAF_ZAMP)),((NAF_ZALP)) et ((NAF_NFS))
!         UTILE QUAND ON BOUCLE SUR PLUSIEURS CAS
!  - PHASE DE SORTIE (POUR RECUPERER DE LA PLACE OU POUR REINITIALISER)
!      1) CALL CLEANNAF
!         DESALLOUE LES TABLEAUX DYNAMIQUES ((NAF_ZTABS)), ((NAF_TFS)),
!           ((NAF_ZAMP)), ((ZALP)) et ((TWIN))
!-------------------------------------------------------------------------
!  VARIABLES INDISPENSABLES
!  NAF_NTERM            : LE NOMBRE MAXIMAL DE FREQUENCES A RECHERCHER
!  NAF_KTABS            : LE NOMBRE D'INTERVALLES ENTRE LES DONNEES
!                         (IL Y A DONC NAF_KTABS+1 DONNEES)
!  NAF_XH               : PAS ENTRE DEUX DONNEES 
!                         DANS UNE UNITE DE TEMPS QUELCONQUE
!                         LA DUREE DE L'INTERVALLE EST ALORS NAF_XH*KTABS
!  NAF_T0               : DATE DE DEPART DES DONNEES
!  NAF_DTOUR            : LONGUEUR D'UN TOUR DE CADRAN
!  NAF_ICPLX            :  0 SI LA FONCTION EST REELLE, 1 SINON
!  NAF_IW               : PARAMETRE DE LA FENETRE D'INTEGRATION
!                            NAF_IW=0  PAS DE FENETRE
!                            NAF_IW=1  FENETRE DE HANNING  etc
!  NAF_ISEC             : DRAPEAU DE LA METHODE DES SECANTES
!                           NAF_ISEC=1 oui,  NAF_ISEC=0 non
!  NAF_NFPRT            : NUMERO DE LA CONNECTION PUR L'IMPRESSION
!  NAF_IPRT             : INDIQUE LE TYPE D'IMPRESSION, EN GENERAL :
!                         -1 PAS D'IMPRESSION
!                         O  MESSAGES D'ERREURS EVENTUELLES
!                         1  RESULTATS
!                         2  (DEBUG)
!  NAF_TOL              : TOLERANCE POUR DETERMINER SI 2 FREQUENCES
!                         SONT IDENTIQUES
!  TABLEAU OU RANGER LE SIGNAL
!  NAF_ZTABS(0:NAF_KTABS)   : TABLEAU DE NAF_KTABS+1 COMPLEXES
!
!  VARIABLES CALCULEES PAR NAFF
!  NAF_UNIANG           : UNITE D'ANGLE = 1RD=UNIANG UNITE D'ANGLE
!                          (SI L'UNITE EST EN SECONDES D'ARC
!                              NAF_UNIANG=206264.80624709D0 )
!                        CALCULE EN FONCTION DE DTOUR
!  NAF_FREFON           : FREQUENCE FONDAMENTALE
!                           NAF_FREFON=2*PI/(KTABS*XH) OU EN SECONDES
!                           NAF_FREFON=360.D0*3600.D0/(KTABS*XH)  
!  NAF_NFS              : NOMBRE DE FREQUENCES ENTIEREMENT DETERMINEES
!                         AINSI QUE LEUR AMPLITUDE
!  NAF_TFS(NBTERM)      : TABLEAU DES FREQUENCES
!  NAF_ZAMP(NBTERM)     : AMPLITUDE COMPLEXE DES TERMES
!  NAF_ZALP(NBTERM,NBTERM) : TABLEAU DE CHANGEMENT DE BASE
!
!-----------------------------------------------------------------------
       MODULE NAFF
          IMPLICIT NONE
          INTEGER,SAVE :: NAF_NERROR,NAF_NFPRT,NAF_IPRT
          INTEGER,SAVE :: NAF_NTERM,NAF_KTABS,NAF_NFS
          INTEGER,SAVE :: NAF_ICPLX,NAF_IW,NAF_ISEC
          REAL (8),SAVE  ::NAF_EPSM,NAF_PI
          REAL (8),SAVE :: NAF_DTOUR,NAF_UNIANG,NAF_FREFON
          REAL (8),SAVE :: NAF_XH,NAF_T0,NAF_TOL
          REAL (8), DIMENSION (:),ALLOCATABLE,SAVE :: NAF_TFS
          COMPLEX (8), DIMENSION (:),ALLOCATABLE,SAVE :: NAF_ZAMP
          COMPLEX (8), DIMENSION (:,:),ALLOCATABLE,SAVE ::NAF_ZALP
          COMPLEX (8), DIMENSION(:),ALLOCATABLE,SAVE :: NAF_ZTABS
!
!-----------------------------------------------------------------------          
! VARIABLES A INITIALISER PAR L'UTILISATEUR AVANT DE LANCER INITNAF
          PUBLIC :: NAF_DTOUR,NAF_KTABS,NAF_XH,NAF_NTERM,NAF_IW
          PUBLIC :: NAF_NFPRT,NAF_IPRT,NAF_TOL
! VARIABLES PUBLIQUES INITIALISEES PAR INITNAF
          PUBLIC :: NAF_EPSM,NAF_PI,NAF_UNIANG,NAF_FREFON
! TABLEAUX PUBLICS CREES PAR INITNAF ET DETRUITES PAR CLEANNAFF
          PUBLIC ::  NAF_TFS,NAF_ZAMP,NAF_ZTABS,NAF_ZALP
!-----------------------------------------------------------------------          
! VARIABLES  A INITIALISER PAR L'UTILISATEUR AVANT DE LANCER NAFF
          PUBLIC :: NAF_T0,NAF_ICPLX,NAF_ISEC
! VARIABLES INITIALISES ET UTILISEES PAR NAFF
          PUBLIC :: NAF_NERROR,NAF_NFS
!-----------------------------------------------------------------------
! ROUTINES PUBLIQUES
          PUBLIC :: INITNAF,CLEANNAF,MFTNAF,PRTABS
          PUBLIC :: INIFRE, INIWIN
          PUBLIC :: CORRECTION
!
!-----------------------------------------------------------------------          
! VARIABLES PRIVEES
          REAL (8),DIMENSION(:),ALLOCATABLE,SAVE :: TWIN, NAF_VT
          REAL (8),SAVE :: AF,BF, NAF_PHISEC0
          COMPLEX (8),DIMENSION(:),ALLOCATABLE,SAVE :: NAF_ZWORK1
          COMPLEX (8),DIMENSION(:),ALLOCATABLE,SAVE :: NAF_ZWORK2
          PRIVATE :: TWIN,AF,BF,NAF_ZWORK1,NAF_ZWORK2,NAF_VT,NAF_PHISEC0                   
! ROUTINES PRIVEES
          PRIVATE :: FFTMAX,FOUR1,PUISS2       
          PRIVATE :: SECANTES,MAXIQUA
          PRIVATE :: FUNC,FUNCP,FREFIN,ZARDYD
          PRIVATE :: DERMODFFNU, PRODERFFNU,PRONUNU, PROFFNU
          PRIVATE :: PHISEC0
!------------------------------------------------------------------------          
          CONTAINS
 
 
!

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!           ROUTINES DE NAFF
!-----------------------------------------------------------------------
      SUBROUTINE INITNAF
      IMPLICIT NONE
!-----------------------------------------------------------------------
!        INITNAFF
!  - effectue les initialisations necessaires :  PI,UNIANG,FREFON,EPSM
!  - alloue les tableaux dynamiques NAF_TFS,NAF_ZAMP,NAF_ZALP,
!    NAF_ZTABS,TWIN, NAF_ZWORK1, NAF_ZWORK2, NAF_VT
!  - initialise NAF_VT
!  - appelle INIWIN
!-----------------------------------------------------------------------
      INTEGER it
!----------------- PREMIERES INITIALISATIONS
      NAF_EPSM = EPSILON(1.D0)
      NAF_PI = ATAN2(1.D0,0.D0)*2
      NAF_UNIANG = NAF_DTOUR/(2*NAF_PI) 
      NAF_FREFON = NAF_DTOUR/(NAF_KTABS*NAF_XH)        
      NAF_TOL=1.D-4
      allocate(NAF_TFS(1:NAF_NTERM),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de NAF_TFS impossible'
        stop
      endif 
      allocate(NAF_ZAMP(1:NAF_NTERM),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de NAF_ZAMP impossible'
        stop
      endif 
      allocate(NAF_ZALP(1:NAF_NTERM,1:NAF_NTERM),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de NAF_ZALP impossible'
        stop
      endif 
      allocate(NAF_ZTABS(0:NAF_KTABS),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de NAF_ZTABS impossible'
        stop
      endif 
      allocate(TWIN(0:NAF_KTABS),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de TWIN impossible'
        stop
      endif 
      allocate(NAF_VT(0:NAF_KTABS),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de NAF_VT impossible'
        stop
      endif
      allocate(NAF_ZWORK1(0:NAF_KTABS),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de NAF_ZWORK1 impossible'
        stop
      endif 
      allocate(NAF_ZWORK2(0:NAF_KTABS),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de NAF_ZWORK2 impossible'
        stop
      endif 
!------------------ initialisation de NAF_VT
      do it=0, NAF_KTABS
       NAF_VT(it) = NAF_T0+it*NAF_XH
      enddo
      call iniwin       
      end SUBROUTINE INITNAF


      subroutine CLEANNAF
!-----------------------------------------------------------------------
!     desalloue les tableaux dynamiques
!-----------------------------------------------------------------------
 
      IF (allocated(NAF_TFS)) THEN
         deallocate(NAF_TFS)
      ENDIF
      IF (allocated(NAF_ZAMP)) THEN
         deallocate(NAF_ZAMP)
      ENDIF
      IF (allocated(NAF_ZALP)) THEN
         deallocate(NAF_ZALP)
      ENDIF
      IF (allocated(NAF_ZTABS)) THEN
         deallocate(NAF_ZTABS)
      ENDIF
      IF (allocated(TWIN)) THEN
         deallocate(TWIN)
      ENDIF
      IF (allocated(NAF_VT)) THEN
         deallocate(NAF_VT)
      ENDIF
      IF (allocated(NAF_ZWORK1)) THEN
         deallocate(NAF_ZWORK1)
      ENDIF
      IF (allocated(NAF_ZWORK2)) THEN
         deallocate(NAF_ZWORK2)
      ENDIF
      end  subroutine CLEANNAF  
      
         

      SUBROUTINE MFTNAF(NBTERM,EPS)
!-----------------------------------------------------------------------
!     MFTNAF
!                CALCULE UNE APPROXIMATION QUASI PERIODIQUE
!                DE LA FONCTION TABULEE DANS NAF_ZTABS(0:NAF_KTABS)
!
!     NBTERM               : NOMBRE DE TERMES RECHERCHES (<= NAF_NTERM)
!
!
!     EPS                 : PRECISION AVEC LAQUELLE ON RECHERCHE 
!                           LES FREQUENCES
!
!-----------------------------------------------------------------------
! TOL est MODIFIE EN NAF_TOL et SON INITIALISATION EST DEPLACE.
      IMPLICIT NONE
! (EPS)
      integer :: NBTERM
      REAL (8) :: EPS     
!
      integer :: I,IFLAG,NUMFR
      REAL (8) :: STAREP,FR,A,B,RM     
!-----------------------------------------------------------------------
      IF (NBTERM .GT.NAF_NTERM) THEN
          WRITE(*,*) 'Nbre de termes cherches trop grand'
          STOP
      ENDIF
      STAREP=ABS(NAF_FREFON)/3
      CALL INIFRE
      DO  I=1,NBTERM
         CALL FFTMAX(FR)
         CALL FREFIN(FR,A,B,RM,STAREP,EPS)
         CALL FRETES(FR,IFLAG,NAF_TOL,NUMFR)
         IF (IFLAG.EQ.0) GOTO 999
         IF (IFLAG.EQ.1) THEN
            CALL GRAMSC(FR,A,B)
            NAF_NFS=NAF_NFS+1
         ENDIF
         IF (IFLAG.EQ.-1) THEN
            CALL MODFRE(NUMFR,A,B)
         ENDIF
      END DO
999   CONTINUE
      END SUBROUTINE MFTNAF 

    
      SUBROUTINE PRTABS(KTABS,ZTABS,IPAS)
!-----------------------------------------------------------------------
!     IMPRESSION DE ZTABS
!-----------------------------------------------------------------------
  
      IMPLICIT NONE
! (KTABS,ZTABS,IPAS)
      integer :: KTABS,IPAS
      complex (8) :: ZTABS(0:KTABS)
!
      integer :: I         
!
      IF (NAF_IPRT.EQ.1) THEN
         DO  I=0,KTABS,IPAS
            WRITE(NAF_NFPRT,1000) I,DREAL(ZTABS(I)),DIMAG(ZTABS(I))
         END DO
      ENDIF
1000  FORMAT (1X,I6,2X,2D20.6)
      END SUBROUTINE PRTABS 


      SUBROUTINE FFTMAX(FR)
!-----------------------------------------------------------------------
!     FFTMAX          MODIF NICE 90 (J. LASKAR)
!                     MODIF 13/12/90
!                     VERSION STANDARD. (FOUR1)
!                     ANALYSE DE FOURIER RAPIDE
!                     RECHERCHE D'UN MAXIMUM DE SPECTRE POUR UNE FREQ. F
!***********************************************************************
!     ON UTILISERA CELA ULTERIEUREMENT (13/12/90)
!     NFR11, NFR22 NUMERO MIN ET MAX DES FREQUENCES RECHERCHEES
!                   NFR1 ET NFR2 INCLUSES
!                -NRTAB<NFR1<NFR2<NRTAB
!
!     RTAB(NRTAB,2) EN 1 AMPLITUDE DE SFREQUENCES POSITIVES
!                      2 AMPLITUDE DES FREQUENCES NEGATIVES
!     DISTANCE ENTRE DEUX LIGNES  FREFON
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
! (FR)
      REAL (8), intent(out) :: FR      
!
      integer :: KTABS2,ISG,IPAS,I,IDV,INDX,IFR
      REAL (8) :: FREFO2
      REAL (8), dimension(:),allocatable :: TAB
      REAL (8), dimension(:),allocatable :: RTAB
!      
      CALL PUISS2(NAF_KTABS+1,KTABS2)
!      
      allocate(TAB(2*KTABS2),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de FFTMAX:TAB impossible'
        stop
      endif
      allocate(RTAB(0:KTABS2-1),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de FFTMAX:RTAB impossible'
        stop
      endif
!      
      FREFO2=(NAF_FREFON*NAF_KTABS)/KTABS2
!******************   
      IF (NAF_IPRT.EQ.1) WRITE(NAF_NFPRT,*) 
     $      'KTABS2= ',KTABS2,'   FREFO2= ', FREFO2
!****************! CALCUL DES FREQUENCES 
      ISG=-1
      IDV=KTABS2
      IPAS=1
      DO 10 I=0,KTABS2-1
         TAB(2*I+1)=DREAL(NAF_ZTABS(I))*TWIN(I)
         TAB(2*I+2)=DIMAG(NAF_ZTABS(I))*TWIN(I)
10    CONTINUE
      CALL FOUR1(TAB,KTABS2,ISG)
      DO 20 I=0,KTABS2-1
         RTAB(I)=SQRT(TAB(2*I+1)**2+TAB(2*I+2)**2)/IDV
20    CONTINUE
      deallocate(TAB)
!**********************        CALL MODTAB(KTABS2,RTAB)
      INDX = MAXLOC(RTAB(0:KTABS2-1),1)-1 
      IF (INDX+1.LE.KTABS2/2) THEN
          IFR = INDX
      ELSE
          IFR = INDX-KTABS2
      ENDIF
      FR = IFR*FREFO2*IPAS
      IF (NAF_IPRT.EQ.1)  WRITE (NAF_NFPRT,*) 
     $    'IFR = ',IFR,' FR = ', FR, 'RTAB = ',RTAB(INDX)
      deallocate(RTAB)
      END SUBROUTINE FFTMAX


      SUBROUTINE FOUR1(DATA,NN,ISIGN)
!-----------------------------------------------------------------------
!************** FOURIER RAPIDE DE NUMERICAL RECIPES
!-----------------------------------------------------------------------
      implicit none
! (DATA,NN,ISIGN)
       integer :: NN,ISIGN
       REAL (8) :: DATA(2*NN)
! 
      integer :: N,I,J,M,MMAX,ISTEP
      REAL (8) :: THETA,WPR,WPI,WR,WI,TEMPR,TEMPI,WTEMP
      N=2*NN
      J=1
      DO 11 I=1,N,2
! en fait N est pair donc I<= N-1      
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*SIN(0.5D0*THETA)**2
        WPI=SIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
!**            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
!**            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END SUBROUTINE FOUR1


      SUBROUTINE PUISS2(NT,N2)
!-----------------------------------------------------------------------
!    CALCULE LA PLUS GRANDE PUISSANCE DE 2
!    CONTENUE DANS NT
!-----------------------------------------------------------------------
      implicit none
!  (NT,N2)    
      INTEGER, intent(out) ::  N2
      integer, intent(in) :: NT
!    
      integer :: n
!            
      N=NT
      IF (N.EQ.0) THEN
         N2=0
         RETURN
      ENDIF
      N2=1
100   CONTINUE
      IF (N.GE.2) THEN
         N2=N2*2
         N=N/2
         GOTO 100
      ENDIF
      END SUBROUTINE PUISS2


      
      SUBROUTINE SECANTES(X,PASS,EPS,XM,IPRT,NFPRT)
!-----------------------------------------------------------------------
!      CALCUL PAR  LA METHODE DES SECANTES D'UN ZERO
!      D'UNE FONCTION FUNC FOURNIE PAR L'UTILISATEUR
!      
!     X, PASS : LE ZERO EST SUPPOSE ETRE DANS X-PASS,X+PASS
!     EPS     : PRECISION AVEC LAQUELLE ON VA CALCULER
!               LE MAXIMUM
!     XM      : ABSCISSE DU  ZERO
!    IPRT     : -1 PAS D'IMPRESSION
!               0: IMPRESSION D'UNE EVENTUELLE ERREUR FATALE
!               1: IMPRESSION SUR LE FICHIER NFPRT
!               2: IMPRESSION DES RESULTATS INTERMEDIAIRES
!     
!    FONC(X) :FONCTION DONNEE PAR L'UTILISATEUR
!    PARAMETRES
!       NMAX  : NOMBRE D'ESSAIS
!       IENC  : 1 : test d'encadrement 0 sinon
!   VARIABLE NERROR
!       0: A PRIORI TOUT S'EST BIEN PASSE
!       1: LE ZERO N'EST PAS TROUVE A EPSI PRES
!          -->   SEULEMENT SI IENC=1
!       2: PENTRE TROP FAIBLE (DIVISION PAR PRESQUE ZERO)
!       3: ECHEC TOTAL 
!  F. Joutel 26/5/98
!  Modif 26/8/98 Enleve l'appel de FONC car dans un module
! il ne semble pas etre possible de passer une fonction
! interne en parametre, remplace FONC par FUNCP
!-----------------------------------------------------------------------
       
       IMPLICIT NONE
! (X,PASS,EPS,XM,FUNCP,IPRT,NFPRT)
       real (8) :: X, PASS,EPS,XM
       integer :: IPRT,NFPRT
!              
       integer, parameter :: NMAX=30, IENC=1
       integer :: I
       REAL (8) :: EPSI,A,B,FA,FB,DELTA,COR
!
       NAF_NERROR=0
       EPSI=MAX(NAF_EPSM,EPS)
       IF (DABS(PASS).GT.EPSI) THEN
         IF (IPRT.GE.1) THEN
           WRITE(*,*) ' AMELIORATION PAR LES SECANTES'
          ENDIF
       ENDIF
       I=0
       A=X-PASS
       B=X
       FB=FUNCP(A)
!      FB=FONC(A)
10    CONTINUE
       IF (DABS(B-A).GT.EPSI) THEN 
          FA=FB
          FB=FUNCP(B)
!           FB=FONC(B)
          DELTA= FB-FA
          IF (IPRT.GE.2) THEN
             WRITE(NFPRT,*) 'SEC: A, B, abs(B-A) ', A, B,ABS(B-A)
            WRITE(NFPRT,*) 'SEC: F(A), F(B) ', FA, FB
          ENDIF
          IF (DABS(DELTA).LE.NAF_EPSM) THEN
            NAF_NERROR=2
            IF (IPRT.GE.1) THEN
              WRITE(NFPRT,*) 'ECHEC DE LA METHODE DES SECANTES'
              WRITE(NFPRT,*) 'DIVISION PAR PRESQUEZERO'
              WRITE(NFPRT,*) 'ON CONTINUE AVEC LA VALEUR TROUVEE'
            ENDIF
            XM=B
            RETURN
          ENDIF
          COR = FB*(A-B)/DELTA
          A = B
          B = B+COR
          I=I+1 
          IF (I.GT.NMAX) THEN
             NAF_NERROR=3
             IF (IPRT.GE.0) THEN
               WRITE(NFPRT,*) ' ECHEC DE LA METHODE DES SECANTES'
               WRITE(NFPRT,*) ' BEAUCOUP TROP D ITERATIONS'
               WRITE(NFPRT,*) ' ON CONTINUE AVEC LA VALEUR INITIALE'
             ENDIF
             XM=X
             RETURN
           ENDIF
           GOTO 10
         ENDIF
!----- !fin de la boucle
       XM=B
       IF (IENC.EQ.1) THEN
            IF (FUNCP(B-EPSI)*FUNCP(B+EPSI).GT.0.D0) THEN
!           IF (FONC(B-EPSI)*FONC(B+EPSI).GT.0.D0) THEN
             NAF_NERROR=1
           ENDIF
       ENDIF
       IF (IPRT.EQ.1) THEN
           WRITE(NFPRT,*) 'POSITION SECANTES',XM,'TROUVEE A',EPSI
          IF (IENC.EQ.1) THEN
             IF (NAF_NERROR.EQ.1) THEN
                WRITE(NFPRT,*) 'SANS GARANTIE'
             ENDIF
             IF (NAF_NERROR.EQ.0) THEN
                WRITE(NFPRT,*) 'AVEC GARANTIE'
             ENDIF
           ENDIF
           WRITE (NFPRT,*)
           WRITE (NFPRT,1000) PASS,XM
       ENDIF
1000   FORMAT (1X, 1D19.9, 1F12.8)
       END SUBROUTINE SECANTES

      SUBROUTINE MAXIQUA(X,PASS,EPS,XM,YM,IPRT,NFPRT)
!-----------------------------------------------------------------------
!      CALCUL PAR INTERPOLATION QUADRATIQUE DU MAX
!      D'UNE FONCTION FONC FOURNIE PAR L'UTILISATEUR
!      
!     X, PASS : LE MAX EST SUPPOSE ETRE SI POSSIBLE 
!               DANS X-PASS,X+PASS
!     EPS     : PRECISION AVEC LAQUELLE ON VA CALCULER
!               LE MAXIMUM
!     XM      : ABSCISSE DU MAX
!     YM      : YM=FUNC(XM)
!    IPRT     : -1 PAS D'IMPRESSION 
!               0: IMPRESSION D'UNE EVENTUELLE ERREUR FATALE
!               1: IMPRESSION SUR LE FICHIER NFPRT
!               2: IMPRESSION DES RESULTATS INTERMEDIAIRES
! 
!    FONC(X)  :FONCTION DONNEE PAR L'UTILISATEUR
! 
!***********************************************************
!  PARAMETRES
!   NMAX     : NOMBRE MAXIMAL D'ITERATIONS
!  VARIABLE (COMMON ASD)
!   NERROR:
!        0 : OK
!        1 : COURBURE BIEN FAIBLE !
!        2 : LA FONCTION EST BIEN PLATE !
!        3 : OSCILLATION  DE L'ENCADREMENT
!        4 : PAS D'ENCADREMENT TROUVE
!        5 : ERREUR FATALE
!***********************************************************
!  ASD VERSION 0.9  J.LASKAR & F.JOUTEL 16/5/97
!  CHANGEMENT POUR UN EPSILON MACHINE CALCULE PAR CALCEPSM
!  JXL 22/4/98
!  Modif 26/8/98 Enleve l'appel de FONC car dans un module
! il ne semble pas etre possible de passer une fonction
! interne en parametre, remplace FONC par FUNCP
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
!  (X,PASS,EPS,XM,YM,FONC,IPRT,NFPRT)
      integer :: IPRT,NFPRT
      REAL (8) :: X,PASS,EPS,XM,YM

      integer, parameter :: NMAX=200, NFATAL=100 
      REAL (8), parameter :: RABS=1.D0
      integer :: NITER,NCALCUL
      REAL (8)  :: PAS,EPSLOC,X1,X2,X3,Y1,Y2,Y3,A,ERR,DX,TEMP
      REAL (8) :: D1,D2   
      
! 
!* ENTRE 0 ET RABS LES ERREURS SONT ENVISAGEES DE FACON ABSOLUE
!* APRES RABS, LES ERREURS SONT ENVISAGEES DE FACON RELATIVE.
!
      EPSLOC=MAX(EPS,SQRT(NAF_EPSM))
!      
!*  AU SOMMET D'UNE PARABOLE Y=Y0-A*(X-X0)**2, UN ECART DE X AUTOUR DE X0
!*  INFERIEUR A EPSLOC DONNE UN ECART DE Y INFERIEUR A A*EPSLOC**2
      IF (IPRT.GE.1) THEN
          WRITE(NFPRT,*) ' ROUTINE DE RECHERCHE DU MAXIMUM :'
      ENDIF
      PAS = PASS
      NITER =0
      NCALCUL=0
      NAF_NERROR=0
      X2 = X
      Y2 = FUNC (X2)
      X1 = X2 - PAS
      X3 = X2 + PAS
      Y1 = FUNC (X1)
      Y3 = FUNC (X3)
!* DECALAGE POUR ENCADRER LE MAXIMUM
10    IF ((NITER.GT.NMAX).OR.(NCALCUL.GT.NFATAL)) THEN
         IF (NCALCUL.GT.NFATAL) THEN
            NAF_NERROR =5
            IF (IPRT.GE.0) THEN
               WRITE(NFPRT,*) ' ERREUR FATALE'
            ENDIF
         ELSE 
           IF (PAS.GE.PASS) THEN
              NAF_NERROR = 4
              IF (IPRT.GE.0) THEN
                WRITE(NFPRT,*) ' PAS D''ENCADREMENT TROUVE'
              ENDIF
           ELSE
              NAF_NERROR =3
              IF (IPRT.GE.0) THEN
                WRITE(NFPRT,*) ' OSCILLATION DE L''ENCADREMENT ?'
              ENDIF
           ENDIF
         ENDIF
         ERR=PAS
         GOTO 100
      ENDIF
      IF ((Y1.GT.Y2).OR.(Y3.GT.Y2)) THEN
            NITER = NITER +1
            IF (Y1.GT.Y3) THEN
                  X2 = X1
                  Y2 = Y1
            ELSE
                  X2 = X3
                  Y2 = Y3
            ENDIF
            X1 = X2 - PAS
            X3 = X2 + PAS
            Y1 = FUNC (X1)
            Y3 = FUNC (X3)
            GOTO 10
      ENDIF
      D1 = Y2-Y1
      D2 = Y3-Y2
      A = ABS(0.5D0*(D2-D1)/PAS**2)
      IF (IPRT.GE.2) THEN
        WRITE(NFPRT,1000) PAS,X2,Y1,Y2,Y3,EPSLOC/SQRT(A)
      ENDIF
!* TEST POUR UN CALCUL SIGNIFICATIF DES VALEURS DE LA FONCTION
      IF ((ABS (D1-D2)).LT.2.D0*NAF_EPSM*MAX(RABS,ABS(Y2))) THEN
        NAF_NERROR= 2
        IF (IPRT.GE.1) THEN
              WRITE (NFPRT,*) ' PLATEAU DE LA FONCTION'
        ENDIF
        ERR=PAS
        GOTO 100
      ENDIF

!* TEST SUR LA FORME DE LA COURBE 
      IF (A.LT.NAF_EPSM*MAX(RABS,ABS(Y2))) THEN
         NAF_NERROR= 1
         IF (IPRT.GE.1) THEN
            WRITE(NFPRT,*) 'COURBE TROP PLATE'
         ENDIF
         ERR=PAS
         GOTO 100
      ENDIF
!* TEST POUR UN CALCUL SIGNIFICATIF AU SOMMET D'UNE PARABOLE
!* ON PEUT IGNORER LE TEST MAIS ON N'A PLUS DE GARANTIE SUR
!* LES VALEURS DE LA FONCTION PLUS TARD 
      IF (PAS.LE.EPSLOC*SQRT(MAX(RABS,ABS(Y2))/A)) THEN
          IF (IPRT.GE.1) THEN
             WRITE(NFPRT,*) 'SORTIE AU SOMMET DE LA PARABOLE' 
          ENDIF
          ERR=EPSLOC*SQRT(MAX(RABS,Y2)/A)
          GOTO 100
        ENDIF 
      DX= 0.5D0*PAS*(Y1-Y3)/(D2-D1)
!* TEST POUR UN NOUVEAU PAS SIGNIFICATIF
      IF (ABS(DX).LT.NAF_EPSM*MAX(RABS,ABS(X2))) THEN
          IF (IPRT.GE.1) THEN
             WRITE(NFPRT,*) 'CORRECTION MINUSCULE'
          ENDIF 
          ERR=PAS
          GOTO 100
      ENDIF
      PAS=ABS(DX) 
      IF (DX.GT.0.D0) THEN
         X1 = X2
         Y1 = Y2
         X2 = X2+PAS
         Y2 = FUNC(X2)
         X3 = X2 + PAS
         Y3 = FUNC(X3)
      ELSE
         X3 = X2
         Y3 = Y2
         X2 = X2-PAS
         Y2 = FUNC(X2)
         X1 = X2 - PAS
         Y1 = FUNC(X1)
        ENDIF
!* A LA SORTIE DE CE CALCUL LE PAS DOIT ETRE DIVISE AU MOINS PAR 2 
!! MAIS ON N'EST PAS ASSURE QUE LES VALEURS FINALES Y1,Y2,Y3 
!* VERIFIENT Y1<=Y2>=Y3, AUSSI, ON REPART AU DEBUT
        NCALCUL=NCALCUL+1
        GOTO 10
100   CONTINUE
      XM = X2
      YM = Y2
      IF (IPRT.GE.1) THEN
         WRITE(NFPRT,1000) PAS,XM,YM
         WRITE(NFPRT,*)' POSITION TROUVEE A',ERR,'PRES' 
         TEMP =  ABS (Y3-Y2)+ ABS(Y2-Y1)
         IF (TEMP.LT.2.D0*NAF_EPSM*MAX(RABS,ABS(Y2))) THEN
            WRITE(NFPRT,*) ' TROUVEE SOUS UN PLATEAU DE LA FONCTION'
         ENDIF 
      ENDIF
1000  FORMAT (1X, 2F12.8, 4D25.16)
      END SUBROUTINE MAXIQUA



      SUBROUTINE INIFRE
!-----------------------------------------------------------------------
!     REMET A ZERO NAF_TFS, NAF_ZAMP,NAF_ZALP et NAF_NFS
!     UTILE QUAND ON BOUCLE SUR PLUSIEURS CAS
!
!-----------------------------------------------------------------------
       
      IMPLICIT NONE
      integer :: I,J
      complex (8) :: ZERO
      ZERO=DCMPLX(0.D0,0.D0)
      DO  I = 1, NAF_NTERM
         NAF_TFS(I) = 0.D0
         NAF_ZAMP(I) = ZERO
         DO  J = 1, NAF_NTERM
            NAF_ZALP(I,J) = ZERO
         END DO
      END DO
      NAF_NFS = 0
      END SUBROUTINE INIFRE


      SUBROUTINE INIWIN
!-----------------------------------------------------------------------
!     INITIALISE LE TABLEAU TWIN POUR UTILISATION D'UNE FENETRE
!     DANS L'ANALYSE DE FOURIER.
!     FENETRE DE HANNING:
!      PHI(T) = (1+COS(PI*T))
!     MODIF LE 13/12/90 (J. LASKAR)
!     MODIF LE 27/1/96 FOR VARIOUS WINDOWS  (J. LASKAR)
!
!     WORKING AREA TWIN(*) SHOULD BE GIVEN IN 
!
!     IW IS THE WINDOW FLAG
!     IW = 0 : NO WINDOW
!     IW = N > 0   PHI(T) = CN*(1+COS(PI*T))**N
!                  WITH CN = 2^N*(N!)^2/(2N)!
!     IW = -1  EXPONENTIAL WINDOW PHI(T) = 1/CE*EXP(-1/(1-T^2))
!
!     MODIF LE 22/4/98 POUR CALCUL EN PLUS DE L'EPSILON MACHINE 
!     EPSM 
!      26/5/98 CORRECTION DE L'ORIGINE DES TABLEAUX (*m/4)
!-----------------------------------------------------------------------
!
       
      IMPLICIT NONE
!
      integer :: I,IT           
      REAL (8) :: CE,T1,T2,TM,T,CN,PIST
!
!---------------------------------------------------------------- 
      CE= 0.22199690808403971891D0
!
!      PI=ATAN2(1.D0,0.D0)*2.D0
      T1=NAF_T0
      T2=NAF_T0+NAF_KTABS*NAF_XH
      TM=(T2-T1)/2
      PIST=NAF_PI/TM
      IF (NAF_IW.EQ.0) THEN
         DO   IT=0,NAF_KTABS
            TWIN(IT)=1.D0 
         ENDDO
      ELSE IF(NAF_IW.GE.0) THEN
         CN = 1.D0
         DO I = 1,NAF_IW
            CN = CN*2.D0*I*1.D0/(NAF_IW+I)
         ENDDO       
         DO IT=0,NAF_KTABS
            T=IT*NAF_XH-TM
            TWIN(IT)=CN*(1.D0+COS(T*PIST))**NAF_IW
         ENDDO
      ELSE IF(NAF_IW.EQ.-1) THEN
         TWIN(0) =0.D0
         TWIN(NAF_KTABS) =0.D0
         DO IT=1,NAF_KTABS-1 
            T=(IT*NAF_XH-TM)/TM
            TWIN(IT)= EXP(-1.D0/(1.D0-T**2))/CE
         ENDDO
      ENDIF
      IF (NAF_IPRT.EQ.1) THEN
!----------------   IMPRESSION TEMOIN
         DO 20 IT=0,NAF_KTABS,NAF_KTABS/20
            WRITE (*,1000) TWIN(IT)*1.D6
20       CONTINUE
      ENDIF
      
      NAF_PHISEC0 = PHISEC0(NAF_IW)
      
1000  FORMAT (1X,F20.3)
      END SUBROUTINE INIWIN



      SUBROUTINE ZARDYD(ZT,N,H,ZOM)
!-----------------------------------------------------------------------
!     CALCULE L'INTEGRALE D'UNE FONCTION TATULEE PAR LA METHODE DE HARDY
!      T(0:N) TABLEAU DOUBLE PRECISION DES VALEURS DE LA FONCTION
!      N = 6*K  ENTIER
!      H PAS ENTRE DEUX VALEURS
!      SOM VALEUR DE L'INTEGRALE SUR L'INTERVALLE [X1,XN]
!      LE PROGRAMME EST EN DOUBLE PRECISION
!               REVU LE 26/9/87 J. LASKAR
!      
!-----------------------------------------------------------------------
      IMPLICIT NONE
!  (ZT,N,H,ZOM)
      integer, intent(in) :: N
      complex (8),intent(in) :: ZT(0:N)
      complex (8),intent(out) :: ZOM
      REAL (8), intent(in) :: H
!
      integer :: ITEST,K,I,INC 
!               
      ITEST=MOD(N,6)
      IF (ITEST.NE.0) THEN
         WRITE(*,*) 'N N''EST PAS UN MULTIPLE DE 6'
         STOP
      ENDIF
      K=N/6
      ZOM=41*ZT(0)+216*ZT(1)+27*ZT(2)+272*ZT(3)
     $          +27*ZT(4)+216*ZT(5)+41*ZT(N)
      INC=0
      DO 10 I=1,K-1
          INC=INC+6
         ZOM=ZOM+82*ZT(INC)+216*ZT(INC+1)+27*ZT(INC+2)+272*ZT(INC+3)
     $           +27*ZT(INC+4)+216*ZT(INC+5)
10    CONTINUE
      ZOM =ZOM*H*6.D0/840.D0
      END SUBROUTINE ZARDYD



      SUBROUTINE PROFFNU(ZFF,NU,VT,ZP)
!-----------------------------------------------------------------------
!     PROFFNU   CALCULE LE PRODUIT SCALAIRE 
!           < ZFF,  EXP(I*NU*T)>
!               SUR L'INTERVALLE [0:KTABS]  T=T0+XH*IT
!              CALCUL NUMERIQUE
!     VT(0:KTABS)   valeur du temps (VT(i) = T0+i*XH)
!     ZFF[0:KTABS]  fonction numerique complexe (in)
!     ZP        PRODUIT SCALAIRE COMPLEXE (out)
!     NU        FREQUENCES EN "/AN (in)
!               24/2/2003 jxl
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL (8), INTENT(IN) :: NU, VT(0:NAF_KTABS)
      COMPLEX (8), INTENT(IN)   :: ZFF(0:NAF_KTABS)
      COMPLEX (8), INTENT(OUT)  :: ZP
!
      real (8) :: OM,H
      complex (8) :: ZI
      
      ZI=DCMPLX(0.D0,1.D0)
!----------! FREQUENCES EN UNITE D'ANGLE PAR UNITE DE TEMPS
      OM = NU/NAF_UNIANG
      NAF_ZWORK1=ZFF*TWIN*EXP(-ZI*OM*VT)
!------------------ TAILLE DU PAS
      H=1.D0/NAF_KTABS
      CALL ZARDYD(NAF_ZWORK1,NAF_KTABS,H,ZP)
      END SUBROUTINE PROFFNU



      SUBROUTINE PRONUNU(NU1,NU2,VT,ZP)
!-----------------------------------------------------------------------
!     PRONUNU   CALCULE LE PRODUIT SCALAIRE 
!           < EXP(I*NU1*T),  EXP(I*NU2*T)>
!               SUR L'INTERVALLE [0:KTABS]  T=T0+XH*IT
!              CALCUL NUMERIQUE
!     VT(0:KTABS)   valeur du temps (VT(i) = T0+i*XH)
!     ZP        PRODUIT SCALAIRE COMPLEXE (out)
!     NU1,NU2        FREQUENCES EN "/AN (in)
!               24/2/2003 jxl
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL (8), INTENT(IN) :: NU1,NU2, VT(0:NAF_KTABS)
      COMPLEX (8), INTENT(OUT)  :: ZP
!
      real(8) :: OM,H
      complex (8) :: ZI
      
      ZI=DCMPLX(0.D0,1.D0)
!----------! FREQUENCES EN UNITE D'ANGLE PAR UNITE DE TEMPS
      OM = (NU1-NU2)/NAF_UNIANG
      NAF_ZWORK1=TWIN*EXP(ZI*OM*VT)
!------------------ TAILLE DU PAS
      H=1.D0/NAF_KTABS
      CALL ZARDYD(NAF_ZWORK1,NAF_KTABS,H,ZP)
      END SUBROUTINE PRONUNU


      SUBROUTINE PRODERFFNU(ZFF,NU,VT,ZP)
!-----------------------------------------------------------------------
!     PRODERFFNU   CALCULE LE PRODUIT SCALAIRE 
!           < ZFF,  I*T*EXP(I*NU*T)> = deriv( <ZFF,EXP(I*NU*T)>,NU)
!               SUR L'INTERVALLE [0:KTABS]  T=T0+XH*IT
!              CALCUL NUMERIQUE
!     VT(0:KTABS)   valeur du temps (VT(i) = T0+i*XH)
!     ZFF[0:KTABS]  fonction numerique complexe (in)
!     ZP        PRODUIT SCALAIRE COMPLEXE (out)
!     NU        FREQUENCES EN "/AN (in)
!               24/2/2003 jxl
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL (8), INTENT(IN) :: NU, VT(0:NAF_KTABS)
      COMPLEX (8), INTENT(IN)   :: ZFF(0:NAF_KTABS)
      COMPLEX (8), INTENT(OUT)  :: ZP
!
      real(8) :: OM,H
      complex (8) :: ZI
      
      ZI=DCMPLX(0.D0,1.D0)
!----------! FREQUENCES EN UNITE D'ANGLE PAR UNITE DE TEMPS
      OM = NU/NAF_UNIANG
      NAF_ZWORK1=-ZI*VT*ZFF*TWIN*EXP(-ZI*OM*VT)
!------------------ TAILLE DU PAS
      H=1.D0/NAF_KTABS
      CALL ZARDYD(NAF_ZWORK1,NAF_KTABS,H,ZP)
      END SUBROUTINE PRODERFFNU



      SUBROUTINE DERMODFFNU(ZFF,NU,VT,DMOD)
!-----------------------------------------------------------------------
!     DERMODFFNU   CALCULE LA DERIVEE DU MODULE DU PRODUIT SCALAIRE 
!         ZPF(NU) = < ZFF,  EXP(I*NU*T)>  
!         SUR L'INTERVALLE [0:KTABS]  T=T0+XH*IT
!           
!     on aura ZP = deriv(ZPF(nu)*conj(ZPF(nu)),nu)
!
!     VT(0:KTABS)   valeur du temps (VT(i) = T0+i*XH)
!     ZFF[0:KTABS]  fonction numerique complexe (in)
!     DMOD        derivee du module (reel)  (out)
!     NU        FREQUENCES EN "/AN (in)
!               24/2/2003 jxl
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL (8), INTENT(IN)  :: NU, VT(0:NAF_KTABS)
      complex (8), INTENT(IN)  :: ZFF(0:NAF_KTABS)
      REAL (8), INTENT(OUT)  :: DMOD
      COMPLEX(8) :: ZPF,ZPFD
!
      CALL  PROFFNU(ZFF,NU,VT,ZPF)
      CALL  PRODERFFNU(ZFF,NU,VT,ZPFD)
      
      DMOD = 2*DREAL(ZPF*CONJG(ZPFD))     
      END SUBROUTINE DERMODFFNU


      FUNCTION FUNC(X)
!--------------------------------------------------------------------------  
!     fonction utile dans MAXIQUA
!     renvoie le module du produit scalaire < f, exp i*X*t>
!--------------------------------------------------------------------------  
      IMPLICIT NONE
      REAL (8) :: X,FUNC
      COMPLEX (8) :: ZP
!
      CALL PROFFNU(NAF_ZTABS,X,NAF_VT,ZP)            
      FUNC=ABS(ZP)
      END FUNCTION FUNC
     
      FUNCTION FUNCP(X)
!--------------------------------------------------------------------------  
!     fonction utile dans SECANTE
!     renvoie la derivee du carre du module du produit scalaire < f, exp i*X*t>
!--------------------------------------------------------------------------  
      IMPLICIT NONE
      REAL (8) :: X,FUNCP
!          
      REAL (8) :: DMOD
!     
      CALL DERMODFFNU(NAF_ZTABS,X,NAF_VT,DMOD) 
      FUNCP=DMOD
      END FUNCTION FUNCP



      SUBROUTINE FREFIN (FR,A,B,RM, RPAS0,RPREC)
!-----------------------------------------------------------------------
!     FREFIN            RECHERCHE FINE DE LA FREQUENCE
!                       FR     FREQUENCE DE BASE EN "/AN  ENTREE ET SORT
!                       A      PARTIE REELLE DE L'AMPLITUDE     (SORTIE)
!                       B      PARTIE IMAGINAIRE DE L'AMPLITUDE (SORTIE)
!                       RM     MODULE DE L'AMPLITUDE            (SORTIE)
!                       RPAS0  PAS INITIAL EN "/AN
!                       RPREC  PRECISION DEMANDEE EN "/AN
!               REVU LE 26/9/87 J. LASKAR
!               MODIF LE 15/11/90 MAX PAR INTERP. QUAD.
! Modif 26/8/98 Enleve les appels aux fonctions func et funcp internes
! au module
!-----------------------------------------------------------------------
      
      IMPLICIT NONE
! (FR,A,B,RM, RPAS0,RPREC)
      REAL (8) FR,A,B,RM,RPAS0,RPREC
!
      REAL (8) X,PASS,EPS,XM,YM
      COMPLEX (8) ZP
! 
      X    = FR
      PASS = RPAS0
      EPS  = RPREC
      CALL MAXIQUA(X,PASS,EPS,XM,YM,NAF_IPRT,NAF_NFPRT)
!************! AMELIORATION DE LA RECHERCHE PAR LE METHODE DES SECANTES
!************!  APPLIQUEE A LA DERIVEE DU MODULE.
      IF (NAF_ISEC.EQ.1) THEN
         X=XM
         CALL SECANTES(X,PASS,EPS,XM,NAF_IPRT,NAF_NFPRT)
         YM=FUNC(XM)
      ENDIF
!------------------ on recalcule l'amplitude avec la frequence trouvee
      CALL PROFFNU(NAF_ZTABS,XM,NAF_VT,ZP)            
         
      FR=XM
      RM=ABS(ZP)
      A=DREAL(ZP)
      B=DIMAG(ZP) 
      IF (NAF_IPRT.EQ.1) THEN
      WRITE (NAF_NFPRT,*)
      WRITE (NAF_NFPRT,1000) FR,RM,A,B
      ENDIF
1000  FORMAT (1X,F10.6,3D19.9)
      END SUBROUTINE FREFIN
      

      SUBROUTINE FRETES (FR,IFLAG,TOL,NUMFR)
!-----------------------------------------------------------------------
!     TEST DE LA NOUVELLE FREQUENCE TROUVEE PAR RAPPORT AUX ANCIENNES.
!     LA DISTANCE ENTRE DEUX FREQ  DOIT ETRE DE NAF_FREFON
!
!     RENVOIE IFLAG =  1 SI LE TEST REUSSIT (ON PEUT CONTINUER)
!             IFLAG =  0 SI LE TEST ECHOUE (IL VAUT MIEUX S'ARRETER)
!             IFLAG = -1 SI TEST < ECART, MAIS TEST/ECART < TOL
!                        (ON RETROUVE PRATIQUEMENT LA MEME FREQUENCE 
!                         D'INDICE NFR)
!     TOL (ENTREE) TOLERANCE ( 1.D-7) EST UN BON EXEMPLE
!     NUMFR (SORTIE) INDICE DE LA FREQUENCE RETROUVEE
!
!-----------------------------------------------------------------------
  
      IMPLICIT NONE
! (FR,IFLAG,TOL,NUMFR)
      integer, intent(out) :: IFLAG,NUMFR
      REAL(8), intent(in) :: FR,TOL
!
      integer :: I
      REAL (8) :: ECART,TEST
!                 
      IFLAG = 1
      ECART = ABS(NAF_FREFON)  
      DO I = 1, NAF_NFS
         TEST = ABS(NAF_TFS(I) - FR)
         IF (TEST.LT.ECART) THEN
            IF (TEST/ECART.LT.TOL) THEN
               IFLAG = -1
               NUMFR   = I
               IF (NAF_IPRT.GE.1) THEN
               WRITE(NAF_NFPRT,*) 'TEST/ECART = ', TEST/ECART ,
     $          '  ON CONTINUE'
               ENDIF
               GOTO 999
            ELSE
               IFLAG = 0
               IF (NAF_IPRT.GE.0) THEN
               WRITE(NAF_NFPRT,*)
               WRITE(NAF_NFPRT,*) 'TEST = ', TEST, 'ECART = ',ECART
               WRITE(NAF_NFPRT,*) 'FREQUENCE   FR = ', FR,
     $           ' TROP PROCHE  DE  ', NAF_TFS(I)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
999   CONTINUE
      END SUBROUTINE FRETES


      SUBROUTINE GRAMSC(FS,A,B)
!-----------------------------------------------------------------------
!     GRAMSC   ORTHONORMALISATION  PAR GRAM-SCHIMDT DE LA BASE DE FONCTI
!               CONSTITUEE DES TERMES PERIODIQUES TROUVES PAR ANALYSE
!               SPECTRALE .
!     FS        FREQUENCE EN "/AN
!     A,B       PARTIES REELLES ET IMAGINAIRES DE L'AMPLITUDE
!
!     AUTRES PARAMETRES TRANSMIS PAR LE COMMON/CGRAM/ZAMP,ZALP,TFS,NFS
!
!     MODIFIE LE 26 9 87 POUR LES FONCTIONS REELLES   J. LASKAR
!     modif 25/2/2003 jxl nouvelles routines  PRONUNU
!-----------------------------------------------------------------------
      
!----------! NAF_TFS    EST LE TABLEAU FINAL DES FREQUENCES
!----------! NAF_ZAMP   TABLEAU DES AMPLITUDES COMPLEXES
!----------! NAF_ZALP   MATRICE DE PASSAGE DE LA BASE ORTHONORMALISEE
!----------! ZTEE       TABLEAU DE TRAVAIL
!----------! NAF_NFS    NOMBRE DE FREQ.DEJA DETERMINEES
      IMPLICIT NONE
! (KTABS,FS,A,B)
      REAL (8), intent(in) :: FS,A,B     
!
      integer :: I, J, K, NF, IT
      REAL (8) DIV
      complex (8), dimension(:),allocatable :: ZTEE
      complex (8) :: ZDIV,ZMUL,ZI,ZEX,ZINC,ZA,ZOM,ZAUX

      allocate(ZTEE(1:NAF_NTERM),stat = NAF_NERROR)
      if (NAF_NERROR.ne.0) THEN
        Write(*,*) 'Allocation de GRAMSC:ZTEE impossible'
        stop
      endif

!----------! CALCUL DE ZTEE(I)=<EN,EI>
      DO 10 I =1,NAF_NFS
        CALL PRONUNU(FS,NAF_TFS(I),NAF_VT,ZTEE(I))
10    CONTINUE
!----------! NF EST LE NUMERO DE LA NOUVELLE FREQUENCE
      NF=NAF_NFS+1
      ZTEE(NF)=DCMPLX(1.D0,0.D0)
!----------! CALCUL DE FN = EN - SOM(<EN,FI>FI) QUI EST ORTHOGONAL AUX F
      NAF_TFS(NF)=FS
      DO 20 K=1,NAF_NFS
         ZAUX=DCMPLX(0.D0,0.D0)
         DO 30 I=K,NAF_NFS
            DO 40 J=1,I
               ZAUX=ZAUX-DCONJG(NAF_ZALP(I,J))*NAF_ZALP(I,K)*ZTEE(J)
40          CONTINUE
30       CONTINUE
         NAF_ZALP(NF,K)=ZAUX
20    CONTINUE
      NAF_ZALP(NF,NF)=DCMPLX(1.D0,0.D0)
!----------! ON REND LA NORME DE FN = 1
      ZDIV=DCMPLX(0.D0,0.D0)
      DO 50 I=1,NF
         ZDIV=ZDIV+DCONJG(NAF_ZALP(NF,I))*ZTEE(I)
50    CONTINUE
      DIV=SQRT(ABS(ZDIV))
      IF (NAF_IPRT.EQ.1) WRITE(NAF_NFPRT,*) 'ZDIV= ',ZDIV,' DIV= ',DIV
      DO 60 I=1,NF
         NAF_ZALP(NF,I)=NAF_ZALP(NF,I)/DIV
60    CONTINUE
!-----------! F1,F2....., FN EST UNE BASE ORTHONORMEE
!-----------! ON RETIRE MAINTENANT A F  <F,FN>FN  (<F,FN>=<F,EN>)
      ZMUL=DCMPLX(A,B)/DIV
      ZI=DCMPLX(0.D0,1.D0)
      DO 70 I=1,NF
         ZOM=NAF_TFS(I)/NAF_UNIANG*ZI
         ZA=NAF_ZALP(NF,I)*ZMUL
!-----------! LES AMPLITUDES DES TERMES SONT CORRIGEES
!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(ZAMP(I)*EXP(ZI*TFS(I)*T) )
         NAF_ZAMP(I)=NAF_ZAMP(I)+ZA
         IF (NAF_IPRT.EQ.1)  WRITE (NAF_NFPRT,1000) NAF_TFS(I),
     $       ABS(NAF_ZAMP(I)),DREAL(NAF_ZAMP(I)),
     $    DIMAG(NAF_ZAMP(I)),
     $    ATAN2(DIMAG(NAF_ZAMP(I)),DREAL(NAF_ZAMP(I)))
1000     FORMAT (1X,F14.9,4D16.8)
!-----------! ON RETIRE LA CONTRIBUTION TOTALE DU TERME TFS(I) DANS TABS
       NAF_ZWORK1 = ZA*EXP(ZOM*NAF_VT)
       IF (NAF_ICPLX.EQ.1) THEN
           NAF_ZTABS=NAF_ZTABS- NAF_ZWORK1 
       ELSE
           NAF_ZTABS=NAF_ZTABS- DREAL(NAF_ZWORK1) 
       ENDIF
70    CONTINUE
      deallocate(ZTEE)
      END SUBROUTINE GRAMSC


      SUBROUTINE MODFRE(NUMFR,A,B)
!-----------------------------------------------------------------------
!     PERMET DE MODIFIER UNE AMPLITUDE DEJA CALCULEE QUAND
!     ON RETROUVE LE MEME TERME, A TOL PRES DANS LE DEVELOPPEMENT
! 
!     NUMFR     INDICE DE LA FREQUENCE EXISTANT DEJA
!     A,B       PARTIES REELLES ET IMAGINAIRES DE L'AMPLITUDE DE LA MODIF 
!               A APPORTER A L'AMPLITUDE DE LA FREQUENCE NUMFR
!
!     AUTRES PARAMETRES TRANSMIS PAR LE COMMON/CGRAM/ZAMP,ZALP,TFS,NFS
!
!     30/4/91
!-----------------------------------------------------------------------
      
!----------! NAF_TFS EST LE TABLEAU FINAL DES FREQUENCES
!----------! NAF_ZAMP   TABLEAU DES AMPLITUDES COMPLEXES
!----------! NAF_ZALP   MATRICE DE PASSAGE DE LA BASE ORTHONORMALISEE
!----------! NAF_NFS    NOMBRE DE FREQ.DEJA DETERMINEES
      IMPLICIT NONE
! (NUMFR,A,B)
      integer, intent(in) :: NUMFR
      REAL (8), intent(in) :: A,B
!            
      integer :: IT
      complex (8) :: ZI,ZOM,ZA,ZEX,ZINC
!      
      ZI = DCMPLX(0.D0,1.D0)
      ZOM=NAF_TFS(NUMFR)/NAF_UNIANG*ZI
      ZA=DCMPLX(A,B)
      IF (NAF_IPRT.EQ.1)  then
            WRITE(*,*) 'CORRECTION DE  IFR = ',NUMFR,
     $               'AMPLITUDE  = ', ABS(ZA)
      ENDIF
!-----------! L' AMPLITUDES DU TERMES EST CORRIGEES
!-----------! ATTENTION ICI (CAS REEL) ON AURA AUSSI LE TERME CONJUGUE
!-----------! QU'ON NE CALCULE PAS. LE TERME TOTAL EST
!-----------!       2*RE(ZAMP(I)*EXP(ZI*TFS(I)*T) )
      NAF_ZAMP(NUMFR)=NAF_ZAMP(NUMFR)+ZA
      IF (NAF_IPRT.EQ.1)  WRITE (NAF_NFPRT,1000) NAF_TFS(NUMFR),
     $   ABS(NAF_ZAMP(NUMFR)),DREAL(NAF_ZAMP(NUMFR)),
     $   DIMAG(NAF_ZAMP(NUMFR)),
     $   ATAN2(DIMAG(NAF_ZAMP(NUMFR)),DREAL(NAF_ZAMP(NUMFR)))
1000     FORMAT (1X,F14.9,4D16.8)
!---------! ON RETIRE LA CONTRIBUTION DU TERME NAF_TFS(NUMFR) DANS TABS
      IF (NAF_ICPLX.EQ.1) THEN
         ZEX = ZA*EXP(ZOM*(NAF_T0-NAF_XH))
         ZINC=EXP(ZOM*NAF_XH)
         CALL  ZTPOW (NAF_KTABS,64,NAF_ZWORK1,ZINC,ZEX)
         DO  IT=0,NAF_KTABS
            NAF_ZTABS(IT)=NAF_ZTABS(IT)- NAF_ZWORK1(IT) 
         END DO
      ELSE
         ZEX = ZA*EXP(ZOM*(NAF_T0-NAF_XH))
         ZINC=EXP(ZOM*NAF_XH)
         CALL  ZTPOW (NAF_KTABS,64,NAF_ZWORK1,ZINC,ZEX)
         DO  IT=0,NAF_KTABS
            NAF_ZTABS(IT)=NAF_ZTABS(IT)- DREAL(NAF_ZWORK1(IT)) 
         END DO
      ENDIF
      END SUBROUTINE MODFRE


      SUBROUTINE CORRECTION(NU0,DNU)
!-----------------------------------------------------------------------
!     NU0 est  donnee comme le zero de LA DERIVEE DU 
!         carre du MODULE DU PRODUIT SCALAIRE 
!         ZPF(NU) = < ZFF,  EXP(I*NU*T)>  
!         SUR L'INTERVALLE [0:KTABS]   
!     
!     On calcule alors une meilleure estimation NU de la frequence principale
!
!     le temps est dans naf_VT(0:KTABS)
!     la fonction est dans NAF_ZTABS(0:KTABS)           
!
!     NU0        FREQUENCES EN "/AN (in)
!     DNU        ecart. La nouvelle frequence est  NU=NU0+DNU
!               24/2/2003 jxl
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL (8), intent(in) :: NU0
      REAL (8), intent(out) :: DNU
!
      COMPLEX (8) :: ZA,ZI,ZD
      REAL(8)  :: OM
!----------------- calcul de l'amplitude complexe
      CALL PROFFNU(NAF_ZTABS,NU0,NAF_VT,ZA)            
!----------------- calcul de la nouvelle fonction dans NAF_ZWORK2(0:KTABS)
      ZI=DCMPLX(0.D0,1.D0)
      OM = NU0/NAF_UNIANG
      NAF_ZWORK1=EXP(ZI*OM*NAF_VT)      
      NAF_ZWORK2=NAF_ZTABS-ZA*NAF_ZWORK1
!----------------- nouvelle frequence
      CALL  PRODERFFNU(NAF_ZWORK2,NU0,NAF_VT,ZD)
      DNU = DREAL(DCONJG(ZA)*ZD)/(ZA*DCONJG(ZA)*NAF_PHISEC0)
      END SUBROUTINE CORRECTION



      real(8) function PHISEC0(IW)
!------------------------------------------------------------------
!   calcul de phi''(0) fenetre d'ordre IW  
!
!  jxl 23/2/2003
!------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: IW
      integer IT
      IF (IW.EQ.-1) then
         PHISEC0=-0.035100738376487704994D0
      elseif (IW.GE.0) then
         PHISEC0=0.D0
         DO IT=IW,1,-1
            PHISEC0=PHISEC0+2.D0/IT**2 
         END DO
         PHISEC0=PHISEC0/NAF_PI**2-1.D0/3.D0
      else
         write (*,*) ' probleme dans PHISEC0  IW=',IW
         stop
      endif
      end function PHISEC0
      
      
      SUBROUTINE ZTPOW (N,N1,ZT,ZA,ZAST)
!-----------------------------------------------------------------------
!     ZTPOW   CALCULE  ZT(I) = ZAST*ZA**I EN VECTORIEL
!             ZT(0:N)
!-----------------------------------------------------------------------
      IMPLICIT NONE
! (N,N1,ZT,ZA,ZAST)
      integer :: N,N1
      complex (8) :: ZT(0:N),ZA,ZAST
!           
      integer :: I,INC,NX,IT,NT
      complex (8) ZT1,ZINC
!      
      IF (N.LT.N1-1) THEN
         WRITE (*,*) 'DANS ZTPOW, N = ', N
         RETURN
      END IF
!----------! 
      ZT(0) = ZAST*ZA
      DO  I = 1,N1-1
         ZT(I) = ZT(I-1)*ZA
      END DO
      ZT1 = ZT(N1-1)/ZAST
      ZINC= 1
      INC =0
      NT = (N+1)/N1  
      DO  IT = 2, NT
         ZINC = ZINC*ZT1
         INC  = INC + N1
         DO  I = 0, N1-1
            ZT(INC +I) = ZT(I)*ZINC
         END DO
      END DO
      ZINC = ZINC*ZT1
      INC  = INC + N1
      NX = N+1-NT*N1
      DO 40 I = 0, NX-1
         ZT(INC +I) = ZT(I)*ZINC
40    CONTINUE
      END SUBROUTINE ZTPOW

      end module naff

