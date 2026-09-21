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
//#include <omp.h>    //To be removed
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"
#include "calculus.h"
#include "intpla.h"


int main(){

      init();
      //transformation_display();
      //Hamiltonian_display();
            
      int i, j;
      typ X[Nd*how_many_planet + 1];

      
      X_init(X);
      //EquilibriumFind(X, 1);
      //AveragedSABAn(2., 50000., 1, X_old, 6);
      //X_init(X_old);
      //SABAn(.04, 20000., 4, X, 10);
      //SABAH1064(3.518764, 1.2*365256000., 10380, X);
      SABAH1064(-.08, -2329., 4, X);
      //SABAH1064(3.518764, 1000., 4, X);
      //X_init(X_old);
      //LibrationCenterFind(X, 1);
      /*LibrationCenterNAFF(X, .0625, 64000., 3, 2, 60);
      SABAn(.2, 64000., 1, X, 9);
      LibrationCenterNAFF(X, .0625, 64000., 3, 2, 60);
      SABAn(.17, 64000., 1, X, 8);
      LibrationCenterNAFF(X, .0625, 64000., 3, 2, 60);
      SABAn(.14, 64000., 1, X, 7);
      LibrationCenterNAFF(X, .0625, 64000., 3, 2, 60);
      SABAn(.12, 64000., 1, X, 6);*/
      //LibrationCenterNAFF(X, .0625, 64000., 3, 2, 60);
      //Renormalization(X);
      //PointPrint(X, -1);
      //SABAn(.22, 64000., 1, X, 9);
      //SABAn(.26, 11814714968., 262144, X_old, 10);
      //EquilibriumFollow(X_old, -epsilon/40000., 20000, 1);
      //SABAn(.25390625, 4.e10, 262144, X_old, 10);
      
      //SABAn(.25390625, 50000., 8, X_old, 10);
      //LibrationCenterFind(X_old, 2);
      //SABAn(.0625, 20000., 1, X_old, 4);
      
      #if 0
      /******** Trying to make a stability map of the 1:2:3 resonance chain. Sigma - delta on Y-axis ********/
      typ Lbd1, Lbd2, Lbd3, lbd1, lbd2, lbd3, g1, g2, g3, D1, D2, D3, Phi, Gamma, Upsilon, Phi_lc, delta_Phi, nu3, nu, delta, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2, B, P;
      typ freq[2][how_many_planet];
      int n_vertical = 503;
      typ * diffusion_rate;
      typ * n1n2;
      typ * n3n2;
      diffusion_rate = (typ *)malloc((2*n_vertical + 1)*sizeof(typ));
      n1n2 = (typ *)malloc((2*n_vertical + 1)*sizeof(typ));
      n3n2 = (typ *)malloc((2*n_vertical + 1)*sizeof(typ));
      if (diffusion_rate == NULL || n1n2 == NULL || n3n2 == NULL){
            fprintf(stderr, "\nError : Cannot allocate memory for array.\n");
            abort();
      }
      typ p = 1.;
      typ q = 3.;
      char LibCenPath[200];
      char out_diffusion[200];
      char out_n1n2[200];
      char out_n3n2[200];
      strcpy(out_diffusion, "/home/atipique/Documents/git/K2138/Stability_map/123/");
      strcat(out_diffusion, "diffusion.txt");
      strcpy(out_n1n2, "/home/atipique/Documents/git/K2138/Stability_map/123/");
      strcat(out_n1n2, "n1n2.txt");
      strcpy(out_n3n2, "/home/atipique/Documents/git/K2138/Stability_map/123/");
      strcat(out_n3n2, "n3n2.txt");
      strcpy(LibCenPath, "/home/atipique/Documents/git/K2138/Stability_map/");
      strcat(LibCenPath, "LibCen_123_aptidal_complete.txt");
      FILE * file_diffusion = fopen(out_diffusion, "w");
      fprintf(file_diffusion, "The columns are delta, B, delta_Phi, dPhi, Phi_lc, Gamma, Upsilon, ..., diff rate at Phi_lc-dPhi, diff rate at Phi_lc, diff rate at Phi_lc+dPhi, ...\n\n");
      FILE * file_n1n2 = fopen(out_n1n2, "w");
      FILE * file_n3n2 = fopen(out_n3n2, "w");
      typ * data = NULL;
      int n_data;
      data = readFromFile(LibCenPath, &n_data);
      if (n_data % 24 != 0){
            fprintf(stderr, "\nError : The file must have exactly 24 columns per line.\n");
            abort();
      }
      int n_lines = n_data/24;
      typ delta_before = 0.;
      int pixel = 0;
      int i_before = -10;
      for (i = 0; i < n_lines; i ++){
            delta = masses[1]*masses[2]*data[24*i + 4]/(m0*m0);
            if (delta - delta_before >= .0003 && i - i_before >= 2){
                  delta_before = delta;
                  i_before = i;
                  pixel ++;
                  printf("\ni = %d, pixel = %d\n\n", i, pixel);
                  // Getting the old coordinates at the libration centers
                  lbd1 = fmod(data[24*i + 12], 2.*M_PI); g1 = fmod(data[24*i + 13], 2.*M_PI); Lbd1 = data[24*i + 14]; D1 = data[24*i + 15];
                  lbd2 = fmod(data[24*i + 16], 2.*M_PI); g2 = fmod(data[24*i + 17], 2.*M_PI); Lbd2 = data[24*i + 18]; D2 = data[24*i + 19];
                  lbd3 = fmod(data[24*i + 20], 2.*M_PI); g3 = fmod(data[24*i + 21], 2.*M_PI); Lbd3 = data[24*i + 22]; D3 = data[24*i + 23];
                  delta_Phi = data[24*i + 5];
                  
                  nu3 = data[24*i + 10];
                  nu  = data[24*i + 11];
                  P = 2.5*max(fabs(2.*M_PI/nu3), fabs(2.*M_PI/nu)); //Integration length. Will never be enough at the separatrix
                  //nu = data[24*i + 10];
                  //P = 2.5*fabs(2.*M_PI/nu); //Integration length. Will never be enough at the separatrix
                  
                  B = data[24*i + 3];
                  Phi_lc = Lbd1/p;  Gamma = (p+q)*Lbd1/p + Lbd2;  Upsilon = Lbd1 + Lbd2 + Lbd3;
                  fprintf(file_diffusion, "%.16g %.16g %.16g %.16g %.16g %.16g %.16g", delta, B, delta_Phi, 3.5*Gamma/((typ) n_vertical)*delta_Phi, Phi_lc, Gamma, Upsilon);
                  // Getting the coordinates of the 1 dof model at the libration centers
                  printf("    j = ");
                  #pragma omp parallel for num_threads(36) private(Phi, Lbd1, Lbd2, Lbd3, X_old, freq, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2) shared(diffusion_rate, n1n2, n3n2)
                  for (j = -n_vertical; j <= n_vertical; j ++){
                        printf("%d,", j);
                        Phi  = Phi_lc + 3.5*Gamma*((typ) j)/((typ) n_vertical)*delta_Phi;
                        Lbd1 = p*Phi;
                        Lbd2 = Gamma - (p+q)*Phi;
                        Lbd3 = q*Phi - Gamma + Upsilon;
                        X_old[1] = lbd1; X_old[2]  = g1; X_old[3]  = Lbd1; X_old[4]  = D1;
                        X_old[5] = lbd2; X_old[6]  = g2; X_old[7]  = Lbd2; X_old[8]  = D2; 
                        X_old[9] = lbd3; X_old[10] = g3; X_old[11] = Lbd3; X_old[12] = D3;                       
                        FundamentalFrequency(0.046875, P, X_old, 2, 2, freq, 2);
                        n1_1 = freq[0][0];
                        n1_2 = freq[1][0];
                        n2_1 = freq[0][1];
                        n2_2 = freq[1][1];
                        n3_1 = freq[0][2];
                        n3_2 = freq[1][2];
                        diffusion_rate[j + n_vertical] = (fabs((n1_1-n1_2)/n1_1) + fabs((n2_1-n2_2)/n2_1) + fabs((n3_1-n3_2)/n3_1))/3.;
                        n1n2[j + n_vertical] = (n1_1 + n1_2)/(n2_1 + n2_2);
                        n3n2[j + n_vertical] = (n3_1 + n3_2)/(n2_1 + n2_2);
                  }
                  printf("\n");
                  for (j = 0; j <= 2*n_vertical; j ++){
                        fprintf(file_diffusion, " %.16g", diffusion_rate[j]);
                        fprintf(file_n1n2, " %.16g", n1n2[j]);
                        fprintf(file_n3n2, " %.16g", n3n2[j]);
                  }
                  fprintf(file_diffusion, "\n");  fprintf(file_n1n2, "\n");  fprintf(file_n3n2, "\n");
            }
      }
      printf("Horizontal pixels = %d\n", pixel);
      free(diffusion_rate); diffusion_rate = NULL;
      free(n1n2); n1n2 = NULL;
      free(n3n2); n3n2 = NULL;
      free(data); data = NULL;
      fclose(file_diffusion); fclose(file_n1n2); fclose(file_n3n2);
      #endif
      
      #if 0
      /******** Trying to make a stability map of the 4:6:9 resonance chain. Sigma - delta on Y-axis ********/
      typ Lbd1, Lbd2, Lbd3, lbd1, lbd2, lbd3, g1, g2, g3, D1, D2, D3, Phi, Gamma, Upsilon, Phi_lc, delta_Phi, nu3, nu, delta, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2, B, P;
      typ freq[2][how_many_planet];
      int n_vertical = 503;
      typ * diffusion_rate;
      typ * n1n2;
      typ * n3n2;
      diffusion_rate = (typ *)malloc((2*n_vertical + 1)*sizeof(typ));
      n1n2 = (typ *)malloc((2*n_vertical + 1)*sizeof(typ));
      n3n2 = (typ *)malloc((2*n_vertical + 1)*sizeof(typ));
      if (diffusion_rate == NULL || n1n2 == NULL || n3n2 == NULL){
            fprintf(stderr, "\nError : Cannot allocate memory for array.\n");
            abort();
      }
      typ p = 2.;
      typ q = 3.;
      char LibCenPath[200];
      char out_diffusion[200];
      char out_n1n2[200];
      char out_n3n2[200];
      strcpy(out_diffusion, "/home/atipique/Documents/git/K2138/Stability_map/469/");
      strcat(out_diffusion, "diffusion.txt");
      strcpy(out_n1n2, "/home/atipique/Documents/git/K2138/Stability_map/469/");
      strcat(out_n1n2, "n1n2.txt");
      strcpy(out_n3n2, "/home/atipique/Documents/git/K2138/Stability_map/469/");
      strcat(out_n3n2, "n3n2.txt");
      strcpy(LibCenPath, "/home/atipique/Documents/git/K2138/Stability_map/");
      strcat(LibCenPath, "LibCen_469_aptidal.txt");
      FILE * file_diffusion = fopen(out_diffusion, "w");
      fprintf(file_diffusion, "The columns are delta, B, delta_Phi, dPhi, Phi_lc, Gamma, Upsilon, ..., diff rate at Phi_lc-dPhi, diff rate at Phi_lc, diff rate at Phi_lc+dPhi, ...\n\n");
      FILE * file_n1n2 = fopen(out_n1n2, "w");
      FILE * file_n3n2 = fopen(out_n3n2, "w");
      typ * data = NULL;
      int n_data;
      data = readFromFile(LibCenPath, &n_data);
      if (n_data % 24 != 0){
            fprintf(stderr, "\nError : The file must have exactly 24 columns per line.\n");
            abort();
      }
      int n_lines = n_data/24;
      typ delta_before = 0.;
      int pixel = 0;
      int i_before = -10;
      for (i = 0; i < n_lines; i ++){
            delta = masses[1]*masses[2]*data[24*i + 4]/(m0*m0);
            if (delta - delta_before >= .00007967 && i - i_before >= 2){
                  delta_before = delta;
                  i_before = i;
                  pixel ++;
                  printf("\ni = %d, pixel = %d\n\n", i, pixel);
                  // Getting the old coordinates at the libration centers
                  lbd1 = fmod(data[24*i + 12], 2.*M_PI); g1 = fmod(data[24*i + 13], 2.*M_PI); Lbd1 = data[24*i + 14]; D1 = data[24*i + 15];
                  lbd2 = fmod(data[24*i + 16], 2.*M_PI); g2 = fmod(data[24*i + 17], 2.*M_PI); Lbd2 = data[24*i + 18]; D2 = data[24*i + 19];
                  lbd3 = fmod(data[24*i + 20], 2.*M_PI); g3 = fmod(data[24*i + 21], 2.*M_PI); Lbd3 = data[24*i + 22]; D3 = data[24*i + 23];
                  delta_Phi = data[24*i + 5];
                  
                  //nu3 = data[24*i + 10];
                  //nu  = data[24*i + 11];
                  //P = 2.5*max(fabs(2.*M_PI/nu3), fabs(2.*M_PI/nu)); //Integration length. Will never be enough at the separatrix
                  nu = data[24*i + 10];
                  P = 2.5*fabs(2.*M_PI/nu); //Integration length. Will never be enough at the separatrix
                  
                  B = data[24*i + 3];
                  Phi_lc = Lbd1/p;  Gamma = (p+q)*Lbd1/p + Lbd2;  Upsilon = Lbd1 + Lbd2 + Lbd3;
                  fprintf(file_diffusion, "%.16g %.16g %.16g %.16g %.16g %.16g %.16g", delta, B, delta_Phi, 3.5*Gamma/((typ) n_vertical)*delta_Phi, Phi_lc, Gamma, Upsilon);
                  // Getting the coordinates of the 1 dof model at the libration centers
                  printf("    j = ");
                  #pragma omp parallel for num_threads(36) private(Phi, Lbd1, Lbd2, Lbd3, X_old, freq, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2) shared(diffusion_rate, n1n2, n3n2)
                  for (j = -n_vertical; j <= n_vertical; j ++){
                        printf("%d,", j);
                        Phi  = Phi_lc + 3.5*Gamma*((typ) j)/((typ) n_vertical)*delta_Phi;
                        Lbd1 = p*Phi;
                        Lbd2 = Gamma - (p+q)*Phi;
                        Lbd3 = q*Phi - Gamma + Upsilon;
                        X_old[1] = lbd1; X_old[2]  = g1; X_old[3]  = Lbd1; X_old[4]  = D1;
                        X_old[5] = lbd2; X_old[6]  = g2; X_old[7]  = Lbd2; X_old[8]  = D2; 
                        X_old[9] = lbd3; X_old[10] = g3; X_old[11] = Lbd3; X_old[12] = D3;                       
                        FundamentalFrequency(0.046875, P, X_old, 2, 2, freq, 2);
                        n1_1 = freq[0][0];
                        n1_2 = freq[1][0];
                        n2_1 = freq[0][1];
                        n2_2 = freq[1][1];
                        n3_1 = freq[0][2];
                        n3_2 = freq[1][2];
                        diffusion_rate[j + n_vertical] = (fabs((n1_1-n1_2)/n1_1) + fabs((n2_1-n2_2)/n2_1) + fabs((n3_1-n3_2)/n3_1))/3.;
                        n1n2[j + n_vertical] = (n1_1 + n1_2)/(n2_1 + n2_2);
                        n3n2[j + n_vertical] = (n3_1 + n3_2)/(n2_1 + n2_2);
                  }
                  printf("\n");
                  for (j = 0; j <= 2*n_vertical; j ++){
                        fprintf(file_diffusion, " %.16g", diffusion_rate[j]);
                        fprintf(file_n1n2, " %.16g", n1n2[j]);
                        fprintf(file_n3n2, " %.16g", n3n2[j]);
                  }
                  fprintf(file_diffusion, "\n");  fprintf(file_n1n2, "\n");  fprintf(file_n3n2, "\n");
            }
      }
      printf("Horizontal pixels = %d\n", pixel);
      free(diffusion_rate); diffusion_rate = NULL;
      free(n1n2); n1n2 = NULL;
      free(n3n2); n3n2 = NULL;
      free(data); data = NULL;
      fclose(file_diffusion); fclose(file_n1n2); fclose(file_n3n2);
      #endif
      
      #if 0
      /******** Trying to make a stability map of the 1:2:3 resonance chain. a3 - a3_lc on Y-axis ********/
      typ Lbd1, Lbd2, Lbd3, lbd1, lbd2, lbd3, g1, g2, g3, D1, D2, D3, Phi, Gamma, Upsilon, Phi_lc, a1, a2, a3, a3_lc, delta_Phi, nu3, nu, delta, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2, B, P, n1_n2;
      typ freq[2][how_many_planet];
      //int every = 14;
      int n_vertical = 431;
      typ * diffusion_rate;
      diffusion_rate = (typ *)malloc((2*n_vertical + 1)*sizeof(typ));
      if (diffusion_rate == NULL){
            fprintf(stderr, "\nError : Cannot allocate memory for array.\n");
            abort();
      }
      typ p = 1.;
      typ q = 3.;
      char LibCenPath[200];
      char out_diffusion[200];
      strcpy(out_diffusion, "/home/atipique/Documents/git/K2138/Stability_map/");
      strcat(out_diffusion, "diffusion.txt");
      strcpy(LibCenPath, "/home/atipique/Documents/git/K2138/Stability_map/");
      strcat(LibCenPath, "LibCen_123_aptidal_complete.txt");
      FILE * file_diffusion = fopen(out_diffusion, "w");
      fprintf(file_diffusion, "The columns are n1/n2, B, delta_Phi, da3, a3_lc, a1, a2, ..., diff rate at a3_lc-da3, diff rate at a3_lc, diff rate at a3_lc+da3, ...\n\n");
      typ * data = NULL;
      int n_data;
      data = readFromFile(LibCenPath, &n_data);
      if (n_data % 24 != 0){
            fprintf(stderr, "\nError : The file must have exactly 24 columns per line.\n");
            abort();
      }
      int n_lines = n_data/24;
      typ n1n2_before = 0.;
      int pixel = 0;
      int i_before = -10;
      typ beta1, beta2, beta3, mu1, mu2, mu3;
      beta1 = m0*masses[1]/(m0 + masses[1]);  beta2 = m0*masses[2]/(m0 + masses[2]);  beta3 = m0*masses[3]/(m0 + masses[3]);
      mu1 = G*(m0 + masses[1]);               mu2 = G*(m0 + masses[2]);               mu3 = G*(m0 + masses[3]);
      for (i = 0; i < n_lines; i ++){
            n1_n2 = data[24*i];
            if (n1_n2 - n1n2_before >= .000003 && i <= 5300 && i - i_before >= 2){
                  n1n2_before = n1_n2;
                  i_before = i;
                  pixel ++;
                  printf("\ni = %d, pixel = %d\n\n", i, pixel);
                  // Getting the old coordinates at the libration centers
                  lbd1 = fmod(data[24*i + 12], 2.*M_PI); g1 = fmod(data[24*i + 13], 2.*M_PI); Lbd1 = data[24*i + 14]; D1 = data[24*i + 15];
                  lbd2 = fmod(data[24*i + 16], 2.*M_PI); g2 = fmod(data[24*i + 17], 2.*M_PI); Lbd2 = data[24*i + 18]; D2 = data[24*i + 19];
                  lbd3 = fmod(data[24*i + 20], 2.*M_PI); g3 = fmod(data[24*i + 21], 2.*M_PI); Lbd3 = data[24*i + 22]; D3 = data[24*i + 23];
                  a1 = Lbd1*Lbd1/(beta1*beta1*mu1);
                  a2 = Lbd2*Lbd2/(beta2*beta2*mu2);
                  a3 = Lbd3*Lbd3/(beta3*beta3*mu3);
                  a3_lc = a3;
                  delta_Phi = data[24*i + 5];
                  nu3 = data[24*i + 10];
                  nu  = data[24*i + 11];
                  P = 2.5*max(fabs(2.*M_PI/nu3), fabs(2.*M_PI/nu)); //Integration length. Will never be enough at the separatrix
                  B = data[24*i + 3];
                  Phi_lc = Lbd1/p;  Gamma = (p+q)*Lbd1/p + Lbd2;  Upsilon = Lbd1 + Lbd2 + Lbd3;
                  fprintf(file_diffusion, "%.16g %.16g %.16g %.16g %.16g %.16g %.16g", n1_n2, B, delta_Phi, .002/((typ) n_vertical)*a3_lc, a3_lc, a1, a2);
                  // Getting the coordinates of the 1 dof model at the libration centers
                  printf("    j = ");
                  #pragma omp parallel for num_threads(36) private(a3, Lbd3, X_old, freq, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2) shared(diffusion_rate)
                  for (j = -n_vertical; j <= n_vertical; j ++){
                        printf("%d,", j);
                        a3 = a3_lc + .002*((typ) j)/((typ) n_vertical)*a3_lc;
                        Lbd3 = beta3*sqrt(mu3*a3);
                        X_old[1] = lbd1; X_old[2]  = g1; X_old[3]  = Lbd1; X_old[4]  = D1;
                        X_old[5] = lbd2; X_old[6]  = g2; X_old[7]  = Lbd2; X_old[8]  = D2; 
                        X_old[9] = lbd3; X_old[10] = g3; X_old[11] = Lbd3; X_old[12] = D3;
                        FundamentalFrequency(0.046875, P, X_old, 2, 2, freq, 2);
                        n1_1 = freq[0][0];
                        n1_2 = freq[1][0];
                        n2_1 = freq[0][1];
                        n2_2 = freq[1][1];
                        n3_1 = freq[0][2];
                        n3_2 = freq[1][2];
                        diffusion_rate[j + n_vertical] = (fabs((n1_1-n1_2)/n1_1)+fabs((n2_1-n2_2)/n2_1)+fabs((n3_1-n3_2)/n3_1))/3.;
                  }
                  printf("\n");
                  for (j = 0; j <= 2*n_vertical; j ++){
                        fprintf(file_diffusion, " %.16g", diffusion_rate[j]);
                  }
                  fprintf(file_diffusion, "\n");
            }
      }
      printf("Horizontal pixels = %d\n", pixel);
      free(diffusion_rate); diffusion_rate = NULL;
      free(data); data = NULL;
      fclose(file_diffusion);
      #endif

      /*typ frequencies[2];
      typ nu1, nu2;
      FundamentalFrequency(0.0625, 40000., X_old, 2, 1, 2, frequencies, 2);
      nu1 = *frequencies;
      nu2 = *(frequencies + 1);
      printf("nu1_1 = %.16g, nu1_2 = %.16g, diffusion index = %.16g\n", nu1, nu2, log10(fabs((nu1-nu2)/nu1)));*/
      
      
      //SABAn(.25390625, 20000000000., 262144, X_old, 10);
      //SABAH1064(.0625, 50000., 40, X_old);
      
      //EquilibriumFind(X_old, 1);
      //SABAn(0.2,  50000., 1, X_old, 9);
      //LibrationCenterNAFF(X_old, 0.0078125, 100000., 2, 2, 40);
      //SABAn(0.2, 50000., 1, X_old, 8);
      //LibrationCenterFollow(X_old, epsilon/400000., 10000, 2);
      //EquilibriumFollow(X_old, epsilon/400000., 200, 2);
      
      /*for (int __ = 1; __ <= Nd*how_many_planet; __ ++){
            X_buff[__] = X_old[__];
      }
      SABAn(0.125, 16000., 1, X_buff, 4);
      for (int _ = 5; _ < 8; _ ++){
            LibrationCenterNAFF(X_old, .0078125, 16000., 5, 2, 55);
            for (int __ = 1; __ <= Nd*how_many_planet; __ ++){
                  X_buff[__] = X_old[__];
            }
            SABAn(0.125, 16000., 1, X_buff, _);
      }*/
      
      //LibrationCenterFind(X_old, 2);
      //SABAn(0.25, 64000., 1, X_old, 10);
      //AveragedSABAn(2., 64000., 1., X_old, 4);
      
      /*
      for (int _ = 0; _ < 3; _ ++){
            LibrationCenterNAFF(X_old, 0.0625,  15000., 5, 1, 50);
      }
      for (int _ = 0; _ < 8; _ ++){
            LibrationCenterNAFF(X_old, 0.03125, 45000., 5, 1, 50);
      }
      SABAn(0.125, 15000., 1, X_old, 7);
      */
      
      //AveragedSABAn(0.125, 25000., 4, X_old, 6);
      
      //SABAn(0.25, 450158004., 2540*10, X_old, 10);
      
      //SABAH1064(.25, 45000., 2, X_old);
      //SABAH1064(.125, 1600000., 64, X_old);
      //SABAH84(.0625, 1600000., 128, X_old);
      /*X_init(X_old);
      for (int _ = 0; _ < 20; _ ++){
            LibrationCenterGradientDescent(X_old, 0.125, 40000., 0.00000018, 4, 1);
      }
      UnaveragedSABAn(0.125, 40000., 1, X_old, 6);*/

  
      deallocation();
      return 0;
}
