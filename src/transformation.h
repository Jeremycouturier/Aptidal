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
/******** @file    transformation.h                                            ********/
/******** @brief   Header file to transformation.c                             ********/
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

#ifndef _TRANSFORMATION_H_
#define _TRANSFORMATION_H_

#include "parameters.h"
#include "structure.h"

extern struct rational transformation[2*how_many_planet + 1][2*how_many_planet + 1]; //matrix M such that (phi, sigma) = M(lambda, -varpi)
extern struct rational transpose_inv [2*how_many_planet + 1][2*how_many_planet + 1]; //tM^-1, that is, (Phi,D) = (tM^-1)(Lambda, D)
extern typ             Transformation[2*how_many_planet + 1][2*how_many_planet + 1]; //matrix M such that (phi, sigma) = M(lambda, -varpi)
extern typ             Transpose_inv [2*how_many_planet + 1][2*how_many_planet + 1]; //tM^-1, that is, (Phi,D) = (tM^-1)(Lambda, D)
extern struct rational rat_l_ij      [  how_many_planet - 1][  how_many_planet + 1]; //Coefficients l_ij
extern struct rational NoverD        [  how_many_planet - 1][  how_many_planet - 1][how_many_planet + 1]; //Coefficients N_r,s^(i) and D_r,s^(i)
extern struct rational rat_c_i       [  how_many_planet - 1];                        //Coefficients c_i
extern int dof                       [  how_many_planet + 1];                        //The indexes j of the degrees of freedom (phi_j; Phi_j)
extern int how_many_dof;                                                             //The number of degrees of freedom of the form (phi_j; Phi_j)
extern int nondof                    [  how_many_planet + 1];                        //The indexes j of the non-degrees of freedom (phi_j; Phi_j)
extern int how_many_nondof;                                                          //The number of non-degrees of freedom of the form (phi_j; Phi_j)


void l_ij_init();


void NoverD_init();


void rat_c_i_init();


void transformation_init();


void transpose_inv_init();


void verification();


void matrix_fill();


void transformation_display();


void Hamiltonian_display();


int GCD(int * As, int k);


int LCM(int * As, int k);
#endif
