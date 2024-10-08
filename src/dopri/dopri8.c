/* dopri8 -- numerical integration of ODE using DOPRI method

    Copyright (C) 1995-2014, IMCCE,CNRS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
    Authors : J. Laskar, F. Joutel, M. Gastineau
    e-mail : laskar@imcce.fr

    History:
        DOPRI8.C (version de dopri 1.5) Version 0.96 
        F. JOUTEL Bureau des Longitudes, Paris (version originale fortran)
        02/09/98 : M. GASTINEAU, rewritten in C 
        17/09/98 : M. GASTINEAU, update of s_integnum_evalarray, 
          integnum_array and dopri8_mem
        24/11/98 : M. GASTINEAU, large update of all functions  
*/
/*------------------------------------------------------------------------*/
  
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "dopri8.h"


#define nd 2
#define Myyerror printf
#define MAX(x,y) (((x) >= (y)) ? (x) : (y))
#define MIN(x,y) (((x) <= (y)) ? (x) : (y))
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* Common Block Declarations                                              */

/*informations relatives a l'execution de dopri */
struct dop8st_ 
{
    int  nfcn;   /* nombre d'evaluation */
    int  nstep;  /* nombre de pas calcule */
    int  naccpt; /* nombre de pas accepte */
    int  nrejct; /* nombre de pas rejete */
};

static struct dop8st_ dop8st_1;


struct coefrk_ {
    double c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, a21, a31, 
	    a32, a41, a43, a51, a53, a54, a61, a64, a65, a71, a74, a75, a76, 
	    a81, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98, a101, a104,
	     a105, a106, a107, a108, a109, a111, a114, a115, a116, a117, a118,
	     a119, a1110, a121, a124, a125, a126, a127, a128, a129, a1210, 
	    a1211, a131, a134, a135, a136, a137, a138, a139, a1310, a1311, b1,
	     b6, b7, b8, b9, b10, b11, b12, b13, bh1, bh6, bh7, bh8, bh9, 
	    bh10, bh11, bh12;
} ;

/*equivaut a coefst dans dopri.f */
static struct coefrk_ coefrk_1= {
        /*C2=*/ 1.E0/18.E0 ,
        /*C3=*/ 1.E0/12.E0 ,
        /*C4=*/ 1.E0/8.E0 ,
        /*C5=*/ 5.E0/16.E0 ,
        /*C6=*/ 3.E0/8.E0 ,
        /*C7=*/ 59.E0/400.E0 ,
        /*C8=*/ 93.E0/200.E0 ,
        /*C9=*/ 5490023248.E0/9719169821.E0 ,
        /*C10=*/ 13.E0/20.E0 ,
        /*C11=*/ 1201146811.E0/1299019798.E0 ,
        /*C12=*/ 1.E0 ,
        /*C13=*/ 1.E0 ,
        /*A21= C2 */ 1.E0/18.E0 ,
        /*A31=*/ 1.E0/48.E0 ,
        /*A32=*/ 1.E0/16.E0 ,
        /*A41=*/ 1.E0/32.E0 ,
        /*A43=*/ 3.E0/32.E0 ,
        /*A51=*/ 5.E0/16.E0 ,
        /*A53=*/ -75.E0/64.E0 ,
        /*A54= -A53 */ 75.E0/64.E0,
        /*A61=*/ 3.E0/80.E0 ,
        /*A64=*/ 3.E0/16.E0 ,
        /*A65=*/ 3.E0/20.E0 ,
        /*A71=*/ 29443841.E0/614563906.E0 ,
        /*A74=*/ 77736538.E0/692538347.E0 ,
        /*A75=*/ -28693883.E0/1125.E6 ,
        /*A76=*/ 23124283.E0/18.E8 ,
        /*A81=*/ 16016141.E0/946692911.E0 ,
        /*A84=*/ 61564180.E0/158732637.E0 ,
        /*A85=*/ 22789713.E0/633445777.E0 ,
        /*A86=*/ 545815736.E0/2771057229.E0 ,
        /*A87=*/ -180193667.E0/1043307555.E0 ,
        /*A91=*/ 39632708.E0/573591083.E0 ,
        /*A94=*/ -433636366.E0/683701615.E0 ,
        /*A95=*/ -421739975.E0/2616292301.E0 ,
        /*A96=*/ 100302831.E0/723423059.E0 ,
        /*A97=*/ 790204164.E0/839813087.E0 ,
        /*A98=*/ 800635310.E0/3783071287.E0 ,
        /*A101=*/ 246121993.E0/1340847787.E0 ,
        /*A104=*/ -37695042795.E0/15268766246.E0 ,
        /*A105=*/ -309121744.E0/1061227803.E0 ,
        /*A106=*/ -12992083.E0/490766935.E0 ,
        /*A107=*/ 6005943493.E0/2108947869.E0 ,
        /*A108=*/ 393006217.E0/1396673457.E0 ,
        /*A109=*/ 123872331.E0/1001029789.E0 ,
        /*A111=*/ -1028468189.E0/846180014.E0 ,
        /*A114=*/ 8478235783.E0/508512852.E0 ,
        /*A115=*/ 1311729495.E0/1432422823.E0 ,
        /*A116=*/ -10304129995.E0/1701304382.E0 ,
        /*A117=*/ -48777925059.E0/3047939560.E0 ,
        /*A118=*/ 15336726248.E0/1032824649.E0 ,
        /*A119=*/ -45442868181.E0/3398467696.E0 ,
        /*A1110=*/ 3065993473.E0/597172653.E0 ,
        /*A121=*/ 185892177.E0/718116043.E0 ,
        /*A124=*/ -3185094517.E0/667107341.E0 ,
        /*A125=*/ -477755414.E0/1098053517.E0 ,
        /*A126=*/ -703635378.E0/230739211.E0 ,
        /*A127=*/ 5731566787.E0/1027545527.E0 ,
        /*A128=*/ 5232866602.E0/850066563.E0 ,
        /*A129=*/ -4093664535.E0/808688257.E0 ,
        /*A1210=*/ 3962137247.E0/1805957418.E0 ,
        /*A1211=*/ 65686358.E0/487910083.E0 ,
        /*A131=*/ 403863854.E0/491063109.E0 ,
        /*A134=*/ -5068492393.E0/434740067.E0 ,
        /*A135=*/ -411421997.E0/543043805.E0 ,
        /*A136=*/ 652783627.E0/914296604.E0 ,
        /*A137=*/ 11173962825.E0/925320556.E0 ,
        /*A138=*/ -13158990841.E0/6184727034.E0 ,
        /*A139=*/ 3936647629.E0/1978049680.E0 ,
        /*A1310=*/ -160528059.E0/685178525.E0 ,
        /*A1311=*/ 248638103.E0/1413531060.E0 ,
        /*B1=*/ 14005451.E0/335480064.E0 ,
        /*B6=*/ -59238493.E0/1068277825.E0 ,
        /*B7=*/ 181606767.E0/758867731.E0 ,
        /*B8=*/ 561292985.E0/797845732.E0 ,
        /*B9=*/ -1041891430.E0/1371343529.E0 ,
        /*B10=*/ 760417239.E0/1151165299.E0 ,
        /*B11=*/ 118820643.E0/751138087.E0 ,
        /*B12=*/ -528747749.E0/2220607170.E0 ,
        /*B13=*/ 1.E0/4.E0 ,
        /*BH1=*/ 13451932.E0/455176623.E0 ,
        /*BH6=*/ -808719846.E0/976000145.E0 ,
        /*BH7=*/ 1757004468.E0/5645159321.E0 ,
        /*BH8=*/ 656045339.E0/265891186.E0 ,
        /*BH9=*/ -3867574721.E0/1518517206.E0 ,
        /*BH10=*/ 465885868.E0/322736535.E0 ,
        /*BH11=*/ 53011238.E0/667516719.E0 ,
        /*BH12=*/ 2.E0/45.E0 
};
/*------------------------------------------------------------------------*/
/* defines locales propres a DOPRI                                        */
#define NMAX  (2000000)
#define FACSEC (.03E0)

/*macro identique a DSIGN du fortran */
/*v0.96 M. GASTINEAU 26/11/98 : optimisation*/
/*#define DSIGN(x,y) (((y)>=0)?fabs(x):-fabs(x))*/ /*remplacee par:*/
#define DSIGN(x,y) (copysign(x,y))



/*------------------------------------------------------------*/
/* dopri8 : routine d'integration numerique.                  */
/* x, hmax, h  sont modifies.                                 */
/* retourne 0  en cas de succes et 1 ou 2 ou 3en cas d'erreur.*/
/*------------------------------------------------------------*/
int gdopri8( int *n,void (*fcn)(int n, double x, double *y, double *f),double *y0,
                   double *x, double *xend, double *eps, double *hmax, double *h)
{
/*#define UROUND DBL_EPSILON*/
#define UROUND 2.23E-16

    double hmaxloc; /* =abs(hmax) */
    bool   reject;  /* indique si le pas est rejete */
    bool   last;    /* indique si on a atteint xend */
    int    nerrtemp;
    double gloc, hnew;
   /* bool   bResultat=true; *//*v0.96 M. GASTINEAU 18/11/98 : mise en commentaire */
    int    nerror=0; /* nombre d'erreur */
    
    double denom;
    double *k1, *k2, *k3, *k4, *k5, *k6, 
	    *k7, *y1, *y;

    double posneg, epsloc, facteur, err, y11s, y12s, xph ;

    int i;

    /* System generated locals */
    double d__1;  
    int naccpt = 0;
    int nrejct = 0;
    int  nfcn = 0;
    int nstep = 0;
    
/* --------------------------------------------------------- */
/*       NUMERICAL SOLUTION OF A SYSTEM OF FIRST ORDER */
/*       ORDINARY DIFFERRENTIAL EQUATIONS Y'=F(X,Y). */
/*       THIS IS AN EMBEDDED RUNGE-KUTTA METHOD OF ORDER (7)8 */
/*       DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL). */
/*       C.F. SECTION II.6 */

/*       PAGE 435 HAIRER, NORSETT & WANNER */

/*       INPUT PARAMETERS */
/*       ---------------- */
/*       N            DIMENSION OF THE SYSTEM ( N.LE.NDMAX) */
/*       FCN          NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
/*                    FIRST DERIVATIVE F(X,Y): */
/*                      SUBROUTINE FCN(N,X,Y,F) */
/*                      REAL*8 X,Y(N),F(N) */
/*                      F(1)=....  ETC. */
/*       X            INITIAL X-VALUE */
/*       XEND         FINAL X-VALUE (XEND-X POSITIVE OR NEGATIVE) */
/*       Y(N)         INITIAL VALUES FOR Y */
/*       EPS          LOCAL TOLERANCE */
/*       HMAX         MAXIMAL STEPSIZE */
/*       H            INITIAL STEPSIZE GUESS */
/*       OUTPUT PARAMETERS */
/*       ----------------- */
/*       Y(N) SOLUTION AT XEND */

/*       EXTERNAL SUBROUTINE (TO BE SUPPLIED BY THE USER) */
/*       ------------------- */
/*       SOLDOPRI       THIS SUBROUTINE IS CALLED AFTER EVERY */
/*                    STEP */
/*                       SUBROUTINE SOLDOPRI(NR,X,Y,N) */
/*                       REAL*8 X,Y(N) */
/*                    FURNISHES THE SOLUTION Y AT THE NR-TH */
/*                    GRID-POINT X (THE INITIAL VALUE IS CON- */
/*                    SIDERED AS THE FIRST GRID-POINT). */
/*                    SUPPLIED A DUMMY SUBROUTINE, IF THE SOLUTION */
/*                    IS NOT DESIRED AT THE INTERMEDIATE POINTS. */
/* -------------------------------------------------------------------- */
/*   Routine reentrante (Utilisation de la variable FIRST) */
/*   Modifie le 9/5/95 pour une recuperation de l'echec (NERROR) */
/*   et la prise en compte du pas XEND-X (variable LAST) */
/*   Modif le 24/5/95 pour un bug sur l'emploi de LAST */
/*   Modif le 31/5/95 pour plus de modularite et suppression des */
/*   boucles numerotees */
/*  Modif le 26/6/95 pour prendre en compte l'u.d.p. local dans la recherc
he*/
/*   du dernier pas */
/*   (*m/4) */
/*   Modif le 10/7/97 pour remplacer 51 par ndmax */
/*   ASD version 1.05 (10/7/97) */
/* -------------------------------------------------------------------- */
/* -------COMMON DOP8ST CAN BE USED FOR STATISTICS */
/* -------   NFCN     NUMBER OF FUNCTION EVALUATIONS */
/* -------   NSTEP    NUMBER OF COMPUTEF STEPS */
/* -------   NACCPT   NUMBER OF ACCEPTED STEPS */
/* -------   NREJCT   NUMBER OF REJECTED STEPS */
/* *******   NERROR   ERROR NUMBER */
/*                      + 0 : OK */
/*                      + 1 : MORE THAN NMAX STEPS */
/*                      + 2 : STEP IS TOO SMALL */

    /* Parameter adjustments */
    y = y0+1;
    y1 = (double*)malloc(sizeof(double)*(*n));
    k1 = (double*)malloc(sizeof(double)*(*n));
    k2 = (double*)malloc(sizeof(double)*(*n));
    k3 = (double*)malloc(sizeof(double)*(*n));
    k4 = (double*)malloc(sizeof(double)*(*n));
    k5 = (double*)malloc(sizeof(double)*(*n));
    k6 = (double*)malloc(sizeof(double)*(*n));
    k7 = (double*)malloc(sizeof(double)*(*n));

    /* Function Body */
/* -------  NMAX    MAXIMAL NUMBER OF STEPS */
/* -------  UROUND  SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.0 */
/* -------          (TO BE ADAPTED BY THE USER) */
/**--- Modif *m/4 (31/5/95) pour adapter le programme a des pas plus petit
s*/
/* -------  FACSEC  Facteur de securite pour ne pas prendre un pas */
/* -------          plus petit que (u.d.p. local)/FACSEC (FACSEC<=1) */

    posneg = DSIGN(1E0, *xend - *x);
/* ------- INITIAL PREPARATIONS */
    hmaxloc = fabs(*hmax);
/* Computing MIN */
/* Computing MAX */
    *h = MIN( MAX(1E-10,fabs(*h)) , hmaxloc);
    *h = DSIGN(*h, posneg);
/* Computing MAX */
    epsloc = MAX((*eps),(UROUND * 13.E0));
    reject = false;
    last = false;
    naccpt = 0;
    nrejct = 0;
    nfcn = 0;
    nstep = 0;
    nerrtemp = 0;
/* *** retire le premier appel */
/*      CALL SOLDOPRI(NACCPT,X,Y,N,h) */
/* -------BASIC INTEGRATION STEP */
  do /*equivaut a label 1 (L1) */
  {
    if (nstep > NMAX || *x + FACSEC * *h == *x) 
    {
	 goto dopri8_fail;
    }
/* *------modif *m/4 (21/6/95) pour respecter l'u.d.p. local */
/* Computing MAX */
/* Computing MIN */
    gloc = MAX( MIN( fabs(*xend) , fabs(*x) ) , 1.E0);
    
    if ((fabs(*xend - *x)) <= UROUND * gloc) 
    {
	 goto dopri8_sortie;
    }
/* *      IF((X-XEND)*POSNEG+UROUND.GT.0D0) GO TO 100 */
/* *------end modif */
/* *------modif *m/4 (9/5/95) pour prendre en compte le dernier pas */
    if ((*x + *h - *xend) * posneg >= 0.) 
    {
	 *h = *xend - *x;
	 last = true;
    }
/* *- -----end modif */
    fcn(*n,*x, y-1, k1-1);
/* *------Modif *m/4 (9/5/95) Le test est rejete au cas de rejet du pas */
/* *       pour eviter de passer dessus lors du dernier pas */
/* * 2       if (nstep.gt.nmax.or.x+0.03D0*h.eq.x) goto 79 */
/* *------End Modif */
L2:
    ++nstep;
/* -----THE FIRST 9 STAGES */
    for (i = 0; i < *n; ++i)
    {
	 y1[i] = y[i] + *h * coefrk_1.a21 * k1[i];
    }
    d__1 = *x + coefrk_1.c2 * *h;
    fcn(*n, d__1, y1-1, k2-1);
    for (i = 0; i < *n; ++i) {
	y1[i ] = y[i] + *h * (coefrk_1.a31 * k1[i] + coefrk_1.a32 * k2[ i ]);
    }
    
    d__1 = *x + coefrk_1.c3 * *h;
    fcn(*n,d__1, y1-1, k3-1);
    for (i = 0; i < *n; ++i) 
    {
	y1[i ] = y[i] + *h * (coefrk_1.a41 * k1[i] + coefrk_1.a43 * k3[i]);
    }
    
    d__1 = *x + coefrk_1.c4 * *h;
    fcn(*n, d__1, y1-1, k4-1);
    for (i = 0; i <  *n; ++i) {
	y1[i] = y[i] + *h * (coefrk_1.a51 * k1[i] + coefrk_1.a53 * k3[ i ] + coefrk_1.a54 * k4[i ]);
    }
    
    d__1 = *x + coefrk_1.c5 * *h;
    fcn(*n, d__1, y1-1, k5-1);
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (coefrk_1.a61 * k1[i] + coefrk_1.a64 * k4[ i ] + coefrk_1.a65 * k5[i]);
    }
    
    d__1 = *x + coefrk_1.c6 * *h;
    fcn(*n, d__1, y1-1, k6-1);
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (coefrk_1.a71 * k1[i ] + coefrk_1.a74 * k4[i] + coefrk_1.a75 * k5[i] 
	        + coefrk_1.a76 * k6[i ]);
    }
    
    d__1 = *x + coefrk_1.c7 * *h;
    fcn(*n, d__1, y1-1, k7-1);
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (coefrk_1.a81 * k1[i ] + coefrk_1.a84 * k4[
		i ] + coefrk_1.a85 * k5[i ] + coefrk_1.a86 * k6[i] + coefrk_1.a87 * k7[i ]);
    }
    
    d__1 = *x + coefrk_1.c8 * *h;
    fcn(*n, d__1, y1-1, k2-1);
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (coefrk_1.a91 * k1[i] + coefrk_1.a94 * k4[
		i] + coefrk_1.a95 * k5[i] + coefrk_1.a96 * k6[i] + coefrk_1.a97 * k7[i] + coefrk_1.a98 * k2[i]);
    }
    
    d__1 = *x + coefrk_1.c9 * *h;
    fcn(*n, d__1, y1-1, k3-1);
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (coefrk_1.a101 * k1[i] + coefrk_1.a104 * k4[i] 
		+ coefrk_1.a105 * k5[i] + coefrk_1.a106 * k6[i ]
		+ coefrk_1.a107 * k7[i] + coefrk_1.a108 * k2[i] 
		+ coefrk_1.a109 * k3[i]);
    }
/* ------- COMPUTE INTERMEDIATE SUMS TO SAVE MEMORY */
    for (i = 0; i < *n; ++i) {
	y11s = coefrk_1.a111 * k1[i] + coefrk_1.a114 * k4[i] + 
		coefrk_1.a115 * k5[i] + coefrk_1.a116 * k6[i] + 
		coefrk_1.a117 * k7[i] + coefrk_1.a118 * k2[i] + 
		coefrk_1.a119 * k3[i];
	y12s = coefrk_1.a121 * k1[i] + coefrk_1.a124 * k4[i] + 
		coefrk_1.a125 * k5[i] + coefrk_1.a126 * k6[i] + 
		coefrk_1.a127 * k7[i ] + coefrk_1.a128 * k2[i] + 
		coefrk_1.a129 * k3[i];
	k4[i] = coefrk_1.a131 * k1[i] + coefrk_1.a134 * k4[i] + 
		coefrk_1.a135 * k5[i] + coefrk_1.a136 * k6[i] + 
		coefrk_1.a137 * k7[i] + coefrk_1.a138 * k2[i] + 
		coefrk_1.a139 * k3[i ];
	k5[i] = coefrk_1.b1 * k1[i] + coefrk_1.b6 * k6[i] + 
		coefrk_1.b7 * k7[i] + coefrk_1.b8 * k2[i] + 
		coefrk_1.b9 * k3[i];
	k6[i] = coefrk_1.bh1 * k1[i] + coefrk_1.bh6 * k6[i] + 
		coefrk_1.bh7 * k7[i ] + coefrk_1.bh8 * k2[i] + 
		coefrk_1.bh9 * k3[i];
	k2[i] = y11s;
	k3[i] = y12s;
    }
/* -----THE LAST 4 STAGES */
    d__1 = *x + coefrk_1.c10 * *h;
    fcn(*n, d__1, y1-1, k7-1);
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (k2[i ] + coefrk_1.a1110 * k7[i ]);
    }
    d__1 = *x + coefrk_1.c11 * *h;
    fcn(*n, d__1, y1-1, k2-1);
    xph = *x + *h;
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (k3[i] + coefrk_1.a1210 * k7[i ] + 
		coefrk_1.a1211 * k2[i ]);
    }
    fcn(*n, xph, y1-1, k3-1);
    for (i = 0; i < *n; ++i) {
	y1[i] = y[i] + *h * (k4[i] + coefrk_1.a1310 * k7[i] + 
		coefrk_1.a1311 * k2[i]);
    }
    fcn(*n, xph, y1-1, k4-1);
    nfcn += 13;
    for (i = 0; i < *n; ++i) {
	k5[i] = y[i] + *h * (k5[i ] + coefrk_1.b10 * k7[i] + 
		coefrk_1.b11 * k2[i ] + coefrk_1.b12 * k3[i] + 
		coefrk_1.b13 * k4[i ]);
	k6[i ] = y[i] + *h * (k6[i ] + coefrk_1.bh10 * k7[i ] + 
		coefrk_1.bh11 * k2[i ] + coefrk_1.bh12 * k3[i]);
    }
/* -----ERROR ESTIMATION */
    err = 0.E0;
    for (i = 0; i < *n; ++i) 
    {
/* Computing MAX */
	 denom = MAX( MAX(1.E-6, fabs( k5[i])) , MAX(fabs(y[i]), UROUND * 2.E0 / epsloc) );
/* Computing 2nd power */
	 d__1 = (k5[i] - k6[i]) / denom;
	 err += d__1 * d__1;
    }
    err = sqrt(err / (double) (*n));
/* -----COMPUTATION OF HNEW */
/* -----WE REQUIRE .333 <=HNEW/W<=6. */
/* Computing MAX */
/* Computing MIN */
    facteur = MAX((1.E0/6.E0),MIN(3.E0, pow( err / epsloc, 1.E0/8.E0) / .9E0 ) );
    hnew = *h / facteur;
    if (err > epsloc)
    {
	 goto L51;
    }
/* -----STEP IS ACCEPTED */
    ++naccpt;
    for (i = 0; i < *n; ++i) 
    {
	 y[i] = k5[i];
    }
    *x = xph;
  /*  i__1 = dop8st_1.naccpt + 1;
      soldopri_(&i__1, x, &y[1], n, h);*/
#ifdef DBUG
     if (niveau>=3)
     {
      for (i=0; i<*n; ++i)
      {
       printf("x=%g,h=%g, y[%d]=%g \n", *x, *h,i, y[i]); 
      } 
     }      
#endif /*DBUG*/ 
    if (fabs(hnew) > hmaxloc)
    {
	 hnew = posneg * hmaxloc;
    }
    if (reject==true)
    {
/* Computing MIN */
	 hnew = posneg * MIN(fabs(hnew),fabs(*h));
    }
    reject = false;
    *h = hnew;
   } while (last==false);
   /*le dernier poin a ete traite => sortir */
   goto dopri8_sortie;

/* -----STEP IS REJECT */
L51:
    reject = true;
/* *----Modif *m/4 (24/5/95) pour corriger l'erreur sur le dernier pas */
/* *----  Meme si c'etait le dernier, cela ne l'est plus */
    last = false;
/* *----endmodif */
    *h = hnew;
    if (naccpt >= 1)
    {
	 ++nrejct;
    }
/* *---- Modif *m/4 (9/5/95) pour prendre en compte le pas XEND-X */
    if (nstep > NMAX || *x + FACSEC * *h == *x)
    {
	 goto dopri8_fail;
    }
/* *---- End Modif */
    --nfcn;
    goto L2;
/* -----FAIL EXIT */
/* ******** a partir d'ici sortie avec erreur. */
/*           On peut eventuellement */
/*          recuperer nerror par le common pour traiter */
/*          l'erreur dans le programme appelant. */
dopri8_fail:
/* *--- Modif *m/4 (9/5/95) pour recuperer l'erreur */
    if (nstep > NMAX)
    {
	 nerrtemp = 1;
    }
    if (*x + FACSEC * *h == *x)
    {
	 nerrtemp = 2;
    }
/* sortie de dopri8 */
dopri8_sortie:
   nerror += nerrtemp;
   if (nerrtemp > 0)
   {
	switch (nerrtemp)
	{
	 case 1:
	  Myyerror("Integration numerique impossible. probleme dans dopri8.\nSortie pour X=%g . TROP D'ETAPES\n",*x);
	  break;
	 case 2: 
	  Myyerror("Integration numerique impossible. probleme dans dopri8.\nSortie pour X=%g\n LE PAS EST DEVENU TROP PETIT, PAS : %g\n",*x,*h);
	  break;
	 default: 
	  Myyerror("Integration numerique impossible. probleme dans dopri8.\nSortie pour X=%g\n",*x);
	  break;
	}
   }
#ifdef DBUG
	if (niveau >= 2)
	{
	   fprintf(stdout, "Routine DOPRI8 - sortie\n");
	   if (niveau >= 3)
	   {
		fprintf(stdout, "RESULTATS\n");
		fprintf(stdout, "X=%g\n",*x);
		fprintf(stdout, "VALEURS DE Y:\n");
		for (i = 0; i < *n; ++i) 
		{
		  fprintf(stdout, "%g\n", y[i]);
		}
	   }
	}
#endif /*DBUG*/
/* *--- End Modif */
    free(y1);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(k7);
   return nerror;
} /* dopri8_ */


