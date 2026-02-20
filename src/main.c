#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h> //To be removed
#include <omp.h>    //To be removed
//#include <gsl/gsl_linalg.h>
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"
#include "calculus.h"
#include "intpla.h"


int main(){

      init();
            
      int i, j;
      typ n       [   how_many_planet + 1];
      typ X_old   [Nd*how_many_planet + 1];
      typ X_buff  [Nd*how_many_planet + 1];
      typ X_new   [Nd*how_many_planet + 1];
      typ X_uv    [Nd*how_many_planet + 1];
      typ xvXu    [Nd*how_many_planet + 1];
      typ X_cart  [Nd*how_many_planet + 1];
      typ dH_old  [Nd*how_many_planet + 1];
      typ dH_polar[Nd*how_many_planet + 1];
      typ dH_rect [Nd*how_many_planet + 1];
      typ epsilon = 0.;
      
      for (i = 1; i <= how_many_planet; i ++){
            epsilon += masses[i]/m0;
      }
      
      X_old_init(X_old);
      
      /******** Trying to make a stability map of the 1:2:3 resonance chain ********/
      typ Lbd1, Lbd2, Lbd3, lbd1, lbd2, lbd3, g1, g2, g3, D1, D2, D3, Phi, Gamma, Upsilon, Phi_lc, delta_Phi, nu3, nu, delta, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2, B, P;
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
      fprintf(file_diffusion, "The columns are delta, B, delta_Phi, dPhi, Phi_lc, Gamma, Upsilon, ..., diff rate at Phi_lc-dPhi, diff rate at Phi_lc, diff rate at Phi_lc+dPhi, ...\n\n");
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
                  B = data[24*i + 3];
                  Phi_lc = Lbd1/p;  Gamma = (p+q)*Lbd1/p + Lbd2;  Upsilon = Lbd1 + Lbd2 + Lbd3;
                  fprintf(file_diffusion, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf", delta, B, delta_Phi, 2.*Gamma/((typ) n_vertical)*delta_Phi, Phi_lc, Gamma, Upsilon);
                  // Getting the coordinates of the 1 dof model at the libration centers
                  printf("    j = ");
                  #pragma omp parallel for num_threads(36) private(Phi, Lbd1, Lbd2, Lbd3, X_old, freq, n1_1, n2_1, n3_1, n1_2, n2_2, n3_2) shared(diffusion_rate)
                  for (j = -n_vertical; j <= n_vertical; j ++){
                        printf("%d,", j);
                        Phi  = Phi_lc + 2.*Gamma*((typ) j)/((typ) n_vertical)*delta_Phi;
                        Lbd1 = p*Phi;
                        Lbd2 = Gamma - (p+q)*Phi;
                        Lbd3 = q*Phi - Gamma + Upsilon;
                        X_old[1] = lbd1; X_old[2]  = g1; X_old[3]  = Lbd1; X_old[4]  = D1;
                        X_old[5] = lbd2; X_old[6]  = g2; X_old[7]  = Lbd2; X_old[8]  = D2; 
                        X_old[9] = lbd3; X_old[10] = g3; X_old[11] = Lbd3; X_old[12] = D3;                       
                        FundamentalFrequency(0.046875, P, X_old, 2, 2, freq, 2);
                        n1_1 = freq[0][0];
                        n1_2 = freq[1][0];
                        //n1_3 = freq[2][0];
                        n2_1 = freq[0][1];
                        n2_2 = freq[1][1];
                        //n2_3 = freq[2][1];
                        n3_1 = freq[0][2];
                        n3_2 = freq[1][2];
                        //n3_3 = freq[2][2];
                        //diffusion_rate[j + n_vertical] = (fabs((n1_1-n1_2)/n1_1)+fabs((n1_2-n1_3)/n1_2)+fabs((n2_1-n2_2)/n2_1)+fabs((n2_2-n2_3)/n2_2)+fabs((n3_1-n3_2)/n3_1)+fabs((n3_2-n3_3)/n3_2))/6.;
                        diffusion_rate[j + n_vertical] = (fabs((n1_1-n1_2)/n1_1)+fabs((n2_1-n2_2)/n2_1)+fabs((n3_1-n3_2)/n3_1))/3.;
                  }
                  printf("\n");
                  for (j = 0; j <= 2*n_vertical; j ++){
                        fprintf(file_diffusion, " %.16lf", diffusion_rate[j]);
                  }
                  fprintf(file_diffusion, "\n");
            }
      }
      printf("Horizontal pixels = %d\n", pixel);
      free(diffusion_rate); diffusion_rate = NULL;
      free(data); data = NULL;
      fclose(file_diffusion);

      /*typ frequencies[2];
      typ nu1, nu2;
      FundamentalFrequency(0.0625, 40000., X_old, 2, 1, 2, frequencies, 2);
      nu1 = *frequencies;
      nu2 = *(frequencies + 1);
      printf("nu1_1 = %.16lf, nu1_2 = %.16lf, diffusion index = %.6lf\n", nu1, nu2, log10(fabs((nu1-nu2)/nu1)));*/
      
      
      //SABAn(0.125, 10000., 1, X_old, 6);
      //SABAH1064(.0625, 50000., 40, X_old);
      
      //EquilibriumFind(X_old, 1);
      //SABAn(0.2,  50000., 1, X_old, 9);
      //LibrationCenterNAFF(X_old, 0.0078125, 100000., 2, 2, 40);
      //SABAn(0.2, 50000., 1, X_old, 8);
      //LibrationCenterFollow(X_old, -epsilon/400000., 1000, 2);
      //EquilibriumFollow(X_old, -epsilon/120000., 200, 0);
      
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
      /*X_old_init(X_old);
      for (int _ = 0; _ < 20; _ ++){
            LibrationCenterGradientDescent(X_old, 0.125, 40000., 0.00000018, 4, 1);
      }
      UnaveragedSABAn(0.125, 40000., 1, X_old, 6);*/

  
      deallocation();
      return 0;
}
