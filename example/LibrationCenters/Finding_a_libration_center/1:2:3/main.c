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
/******** @file    main.c                                                      ********/
/******** @brief   The main file of Aptidal                                    ********/
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"
#include "calculus.h"
#include "intpla.h"


int main(){

      init();
      transformation_display();
      //Hamiltonian_display();
            
      int i, j;
      typ X[Nd*how_many_planet + 1];
      
      X_init(X);
      
      /******** Finding a fixed point of the averaged model ********/
      EquilibriumFind(X, 1);
      
      /******** Finding a libration center of the complete Hamiltonian ********/
      LibrationCenterFind(X, .03125, 40000., 5, 2, 130, 2.);   //Integrating from t=-40000 to t=40000 with a SABA2, a Hanning filter of order 5, with 130 harmonics, until amplitude < 2°                                                          
      LibrationCenterFind(X, .03125, 40000., 5, 2, 200, .001); //Refining to amplitude < 0.001° with 200 harmonics
      
      printf("\nAfter normalization of the semi-major axes, the libration center is at (this is one particular point on the periodic orbit):\n");
      Renormalization(X); //Normalizing to a_1 = 1
      PointPrint(X, 0);

      SABAH1064(.04, 16000., 1, X); //Integrating the libration center.
      deallocation();
      return 0;
}
