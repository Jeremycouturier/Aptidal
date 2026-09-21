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
/******** @file    coefficients.h                                              ********/
/******** @brief   Header file to coefficients.c                               ********/
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

#ifndef _COEFFICIENTS_H_
#define _COEFFICIENTS_H_

/******** Coefficients of the Hamiltonian ********/
extern typ C_ppp1_1, C_ppp1_2, C_ppp1_3, C_ppp1_4, C_ppp1_5, C_ppp1_6, C_ppp1_7, C_ppp1_8, C_ppp1_9, C_ppp1_10, C_ppp1_11, C_ppp1_12, C_ppp1_13, C_ppp1_14, C_ppp1_15; //Resonance p:p+1
extern typ C_ppp2_1, C_ppp2_2, C_ppp2_3; //Resonance p:p+2
extern typ C_ppp3_1, C_ppp3_2, C_ppp3_3, C_ppp3_4; //Resonance p:p+3
extern typ C_00_1, C_00_2, C_00_3, C_00_4, C_00_5, C_00_6, C_00_7, C_00_8, C_00_9, C_00_10; //Resonance 0:0

/******** Array of pointers towards the functions resonance_pq ********/
extern void (*resonances[10][10])(typ alp, typ mi);
extern typ  Cppq[how_many_planet + 1][how_many_planet + 1][32];

#if second_mass_bool
extern int all2pla[16][2];
#endif

void resonance_init();


void Cppq_init();


void resonance_00(typ alp, typ mi);


void resonance_12(typ alp, typ mi);


void resonance_23(typ alp, typ mi);


void resonance_34(typ alp, typ mi);


void resonance_45(typ alp, typ mi);


void resonance_56(typ alp, typ mi);


void resonance_67(typ alp, typ mi);


void resonance_78(typ alp, typ mi);


void resonance_89(typ alp, typ mi);


void resonance_13(typ alp, typ mi);


void resonance_35(typ alp, typ mi);


void resonance_57(typ alp, typ mi);


void resonance_79(typ alp, typ mi);


void resonance_14(typ alp, typ mi);


void resonance_25(typ alp, typ mi);


void resonance_47(typ alp, typ mi);


void resonance_58(typ alp, typ mi);
#endif
