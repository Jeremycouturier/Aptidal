#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <gsl/gsl_linalg.h>
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"


int main(){


      /*typ alpha = 0.43;
      struct pairOfReal b_12_01 = b_12(alpha, 1.0e-15);
      struct pairOfReal b_32_01 = b_k2(alpha, b_12_01, 3);
      struct pairOfReal b_52_01 = b_k2(alpha, b_32_01, 5);
      struct pairOfReal b_72_01 = b_k2(alpha, b_52_01, 7);
      
      typ C1 = -0.5*b_12_01.fst;                                               //Constant term in the inequality (0,0)
      typ C2 = -3.0/8.0*alpha*b_32_01.fst+(0.25+0.25*alpha*alpha)*b_32_01.snd; //Term propto XXb'+XbX' in the inequality (0,0)
      typ C3 = -0.125*alpha*b_32_01.snd;                                       //Term propto XXb+X'Xb' in the inequality (0,0)
      typ C13 = -9.0/128.0*alpha*alpha*b_52_01.fst+3.0/64.0*alpha*(1.0-alpha*alpha)*b_52_01.snd; //Term propto X^2Xb^2
      typ Cbeaucoup = 5.0/384.0*alpha*(1.0-alpha*alpha*alpha*alpha)*b_72_01.snd+55.0/512.0*alpha*alpha*alpha*b_72_01.snd
      -35.0/768.0*alpha*alpha*b_72_01.fst-85.0/768.0*alpha*alpha*alpha*alpha*b_72_01.fst; //Term propto X^3Xb^3
      
      printf("C1 = %.15lf\nC2 = %.15lf\nC3 = %.15lf\nC13 = %.15lf\nCbeaucoup = %.15lf\n", C1, C2, C3, C13, Cbeaucoup);*/
      
      
      initialization();
      chain_validity();
      array_init();
      array2rational();
      resonance_init();
      if (how_many_resonant > 2){
            l_ij_init();
            NoverD_init();
            rat_c_i_init();
      }
      transformation_init();
      transpose_inv_init();
      verification();
      transformation_display();
      Cppq_init();
      Hamiltonian_display();
      typ alpha = 0.894;
      resonances[1][2](alpha,masses[1]);
      /*resonance_00(alpha);
      printf("C_00_1 = %.14lf\nC_00_2 = %.14lf\nC_00_3 = %.14lf\nC_00_4 = %.14lf\nC_00_5 = %.14lf\nC_00_6 = %.14lf\nC_00_7 = %.14lf\nC_00_8 = %.14lf\nC_00_9 = %.14lf\nC_00_10 = %.14lf\n"
      ,C_00_1, C_00_2, C_00_3, C_00_4, C_00_5, C_00_6, C_00_7, C_00_8, C_00_9, C_00_10);*/
      /*resonance_57(alpha,masses[1]);
      printf("C_ppp2_1 = %.14lf\nC_ppp2_2 = %.14lf\nC_ppp2_3 = %.14lf\n" ,C_ppp2_1, C_ppp2_2, C_ppp2_3);*/
      /*resonance_58(alpha,masses[1]);
      printf("C_ppp3_1 = %.14lf\nC_ppp3_2 = %.14lf\nC_ppp3_3 = %.14lf\nC_ppp3_4 = %.14lf\n" ,C_ppp3_1, C_ppp3_2, C_ppp3_3, C_ppp3_4);*/
      
      int i,j;
      for (i = 1; i < how_many_resonant + 1; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  if (i != j){
                        printf("k_%d%d = %d\n", i, j, k_ij[i][j]);
                  }
            }
      }
      for (i = 1; i < how_many_resonant; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  if (j > i){
                        printf("q_%d%d = %d\n", i, j, q_ij[i][j]);
                  }
            }
      }
      for (i = 1; i < how_many_resonant; i++){
            printf("k_%d = %d\n", i, k_i[i]);
            printf("k_star_%d = %d\n", i, k_star_i[i]);
            printf("q_%d = %d\n", i, q_i[i]);
      }
      
      for (i = 1; i < how_many_resonant + 1; i++){
            printf("subchain[%d] = %d\n", i, subchain[i]);
      }
      for (i = 1; i < how_many_missed + 1; i++){
            printf("chain_missed[%d] = %d\n", i, chain_missed[i]);
      }
      for (i = 1; i < how_many_resonant + 1; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  if (i != j){
                        printf("k_%d%d = ",i,j);
                        ratprint(rat_k_ij[i][j]);
                        printf("       ");
                  }
            }
            printf("\n");
      }
      for (i = 1; i < how_many_resonant - 1; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  printf("l_%d%d = ",i,j);
                  ratprint(rat_l_ij[i][j]);
                  printf("       ");
            }
            printf("\n");
      }
      int r,s;
      for (i = 1; i < how_many_resonant - 1; i++){
            for (r = 1; r <= i; r++){
                  for (s = r + 1; s < how_many_resonant + 1; s++){
                        printf("NoverD_%d%d^(%d) = ", r, s, i);
                        ratprint(NoverD[i][r][s]);
                        printf("     ");
                  }
            }
            printf("\n");
      }
      for (i = 1; i < how_many_resonant - 1; i++){
            printf("c_%d = ", i);
            ratprint(rat_c_i[i]);
            printf("     ");
      }
      printf("\n");
      for (i = 1; i < 2*how_many_planet + 1; i++){
            for (j = 1; j < 2*how_many_planet + 1; j++){
                  printf("M_%d,%d=", i, j);
                  ratprint(transformation[i][j]);
                  printf("    ");
            }
            printf("\n");
      }
      printf("\n");
      for (i = 1; i < 2*how_many_planet + 1; i++){
            for (j = 1; j < 2*how_many_planet + 1; j++){
                  printf("Q_%d,%d=", i, j);
                  ratprint(transpose_inv[i][j]);
                  printf("    ");
            }
            printf("\n");
      }
      
      
      /*printf("\nPert = \n");
      printf("       %.14lf * ", C_ppp1_1); printf(string1);printf("\n");
      printf("       %.14lf * ", C_ppp1_2); printf(string2);printf("\n");
      printf("       %.14lf * ", C_ppp1_3); printf(string3);printf("\n");
      printf("       %.14lf * ", C_ppp1_4); printf(string4);printf("\n");
      printf("       %.14lf * ", C_ppp1_5); printf(string5);printf("\n");
      printf("       %.14lf * ", C_ppp1_6); printf(string6);printf("\n");
      printf("       %.14lf * ", C_ppp1_7); printf(string7);printf("\n");
      printf("       %.14lf * ", C_ppp1_8); printf(string8);printf("\n");
      printf("       %.14lf * ", C_ppp1_9); printf(string9);printf("\n");
      printf("       %.14lf * ", C_ppp1_10); printf(string10);printf("\n");
      printf("       %.14lf * ", C_ppp1_11); printf(string11);printf("\n");
      printf("       %.14lf * ", C_ppp1_12); printf(string12);printf("\n");
      printf("       %.14lf * ", C_ppp1_13); printf(string13);printf("\n");
      printf("       %.14lf * ", C_ppp1_14); printf(string14);printf("\n");
      printf("       %.14lf * ", C_ppp1_15); printf(string15);printf("\n");*/
      
      
      /******** GSL test ********/
      double a_data[] = { 1.0, 0.6, 0.0,
                         0.0, 1.5, 1.0,
                         0.0, 1.0, 1.0 };
     /*
      * Inverse is
      *    1  -1.2   1.2
      *    0   2.0  -2.0
      *    0  -2.0   3.0
      */
      double inva[9];
     
      gsl_matrix_view m
          = gsl_matrix_view_array(a_data, 3, 3);
      gsl_matrix_view inv
          = gsl_matrix_view_array(inva,3,3);
      gsl_permutation * p = gsl_permutation_alloc (3);
 
      printf("The matrix is\n");
      for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                  printf(j==2?"%6.3f\n":"%6.3f ", gsl_matrix_get(&m.matrix,i,j));
 
      gsl_linalg_LU_decomp (&m.matrix, p, &s);    
      gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
 
      printf("The inverse is\n");
      for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                  printf(j==2?"%6.3f\n":"%6.3f ",gsl_matrix_get(&inv.matrix,i,j));
 
      gsl_permutation_free (p);

      deallocation();
      return 0;

}
