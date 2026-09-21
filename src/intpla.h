/*
               AAA                PPPPPPPPPPPPPPPPP    TTTTTTTTTTTTTTTTTTTTTTT IIIIIIIIII DDDDDDDDDDDDD                   AAA                LLLLLLLLLLL             
              A:::A               P::::::::::::::::P   T:::::::::::::::::::::T I::::::::I D::::::::::::DDD               A:::A               L:::::::::L             
             A:::::A              P::::::PPPPPP:::::P  T:::::::::::::::::::::T I::::::::I D:::::::::::::::DD            A:::::A              L:::::::::L             
            A:::::::A             PP:::::P     P:::::P T:::::TT:::::::TT:::::T II::::::II DDD:::::DDDDD:::::D          A:::::::A             LL:::::::LL             
           A:::::::::A              P::::P     P:::::P TTTTTT  T:::::T  TTTTTT   I::::I     D:::::D    D:::::D        A:::::::::A              L:::::L               
          A:::::A:::::A             P::::P     P:::::P         T:::::T           I::::I     D:::::D     D:::::D      A:::::A:::::A             L:::::L               
         A:::::A A:::::A            P::::PPPPPP:::::P          T:::::T           I::::I     D:::::D     D:::::D     A:::::A A:::::A            L:::::L               
        A:::::A   A:::::A           P:::::::::::::PP           T:::::T           I::::I     D:::::D     D:::::D    A:::::A   A:::::A           L:::::L               
       A:::::A     A:::::A          P::::PPPPPPPPP             T:::::T           I::::I     D:::::D     D:::::D   A:::::A     A:::::A          L:::::L               
      A:::::AAAAAAAAA:::::A         P::::P                     T:::::T           I::::I     D:::::D     D:::::D  A:::::AAAAAAAAA:::::A         L:::::L               
     A:::::::::::::::::::::A        P::::P                     T:::::T           I::::I     D:::::D     D:::::D A:::::::::::::::::::::A        L:::::L               
    A:::::AAAAAAAAAAAAA:::::A       P::::P                     T:::::T           I::::I     D:::::D    D:::::D A:::::AAAAAAAAAAAAA:::::A       L:::::L         LLLLLL
   A:::::A             A:::::A    PP::::::PP                 TT:::::::TT       II::::::II DDD:::::DDDDD:::::D A:::::A             A:::::A    LL:::::::LLLLLLLLL:::::L
  A:::::A               A:::::A   P::::::::P                 T:::::::::T       I::::::::I D:::::::::::::::DD A:::::A               A:::::A   L::::::::::::::::::::::L
 A:::::A                 A:::::A  P::::::::P                 T:::::::::T       I::::::::I D::::::::::::DDD  A:::::A                 A:::::A  L::::::::::::::::::::::L
AAAAAAA                   AAAAAAA PPPPPPPPPP                 TTTTTTTTTTT       IIIIIIIIII DDDDDDDDDDDDD    AAAAAAA                   AAAAAAA LLLLLLLLLLLLLLLLLLLLLLLL
*/

/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    intpla.h                                                    ********/
/******** @brief   Header file to intpla.c                                     ********/
/******** @author  Jérémy COUTURIER <jeremycouturier.com>                      ********/
/********                                                                      ********/
/******** @section LICENSE                                                     ********/
/******** Copyright (c) 2026 Jérémy COUTURIER                                  ********/
/********                                                                      ********/
/******** Aptidal is free software. You can redistribute it and/or modify      ********/
/******** it under the terms of the GNU General Public License as published by ********/
/******** the Free Software Foundation, either version 3 of the License, or    ********/
/******** (at your option) any later version.                                  ********/
/********                                                                      ********/
/******** Aptidal is distributed in the hope that it will be useful,           ********/
/******** but WITHOUT ANY WARRANTY; without even the implied warranty of       ********/
/******** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ********/
/******** GNU General Public License for more details.                         ********/
/********                                                                      ********/
/******** You should have received a copy of the GNU General Public License    ********/
/******** along with Aptidal. If not, see <http://www.gnu.org/licenses/>.      ********/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

#ifndef _INTPLA_H_
#define _INTPLA_H_

#include "parameters.h"
#include "structure.h"

extern typ nu_fast; //The fast frequency
extern typ nu_reso; //The frequency of the periapses
extern typ avgs[1 + 4*how_many_planet]; //Index 4*i - 3 (resp. -2, -1, -0) contains the average of phi_i (resp. v_i, Phi_i, u_i)


void ell2cart(typ a, typ e, typ i, typ E, typ varpi, typ Omega, typ mu, typ * cart);


void cart2ell(typ * cart, typ * alkhqp, typ mu);


typ mean2eccentric(typ l, typ k, typ h);


void prods(typ a, typ c, const typ * X, typ b, typ d, const typ * Y, typ * R, typ * S);


void newt(typ DM, typ A, typ B, typ * const p_X, typ * const p_C, typ * const p_S, typ * const EXP, typ * const CM1, typ * const SMX, int bexitonerror);


void kepsaut(typ * cart, typ mu, typ dt);


void exp_tau_LB(typ tau, typ * X_cart);


#if tides_bool
void exp_tau_LHt(typ * X_cart, typ tau, int planet);
#endif


void SABAn(typ tau, typ T, int output_step, typ * X_old, int n);


void SABAH1064(typ tau, typ T, int output_step, typ * X_old);


typ Hamiltonian(typ * X_cart);


void FundamentalFrequency(typ tau, typ T, typ * X_old, int n, int n_freq, typ freq[][how_many_planet], int Hanning_order);


void get_frequencies(typ tau, typ T, typ * X_old, int n);


typ UnaveragedSABAn_NAFF(typ tau, typ T, int Hanning_order, typ * X_uv, typ * X_old, int n, int how_many_harmonics);


int UnaveragedSABAn_amplitude(typ tau, typ T, typ * X_new_min, typ * X_new_max, typ * X_old, int n);


void Renormalization(typ * X_old);


void get_averaged_n(typ * n);


void get_n(typ * X_cart, typ * n);


void ConstantParameter(typ * X_new, typ * X_uv);


void LibrationCenterFind(typ * X_old, int precision);


void LibrationCenterNAFF(typ * X_old, typ tau, typ T, int Hf, int N, int Hr);


void LibrationCenterFollow(typ * X_old, typ dG, int Npoints, int precision);


void PeriodicOrbitFind(typ * X_old);


void LibrationCenterGradientDescent(typ * X_old, typ tau, typ T, typ LearningRate, int N, int sigma);
#endif
