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
/******** @file    calculus.h                                                  ********/
/******** @brief   Header file to calculus.c                                   ********/
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

#ifndef _CALCULUS_H_
#define _CALCULUS_H_

#include "parameters.h"
#include "structure.h"

extern typ X_old_t0[4*how_many_planet + 1]; 
extern typ X_new_t0[4*how_many_planet + 1];
extern typ X_uv_t0 [4*how_many_planet + 1];


typ fast_pow(typ x, int power);


void dHdold(typ * dH, typ * X_old, int KP);


void dHdnew(typ * dH_polar, typ * dH_rect, typ * dH_old, typ * X_new, typ * X_uv);


void old2new(typ * X_old, typ * X_new, typ * X_uv);


void new2old(typ * X_old, typ * X_new, typ * X_uv);


void canonical2nonCanonical(typ * X_cart);


void nonCanonical2canonical(typ * X_cart);


#if (toInvar_bool && _3D_bool)
void toInvar(typ * X_cart);
#endif


void X_init(typ * X_old);


void nonDofReinit(typ * X_new, typ * X_uv);


void exp_tau_LB_Ralston(typ tau, typ * X_old);


void AveragedSABAn(typ tau, typ T, int output_step, typ * X_old, int n);


void SABAn_average(typ tau, typ T, int Hanning_order, typ * X_uv_mean, typ * X_old, int n);


void RK2(typ tau, typ T, int output_step);


typ AveragedHamiltonian(typ * X_old);


void PointPrint(typ * X_old, int iter);


int EquilibriumFind(typ * X_old, int precision);


void EquilibriumFollow(typ * X_old, typ dG, int Npoints, int precision);


void EquilibriumFindUntil(typ * X_old, int precision);
#endif
