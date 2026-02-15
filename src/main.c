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
      typ Lbd1, Lbd2, Lbd3, lbd1, lbd2, lbd3, g1, g2, g3, D1, D2, D3, Phi, Gamma, Upsilon, Phi_lc, delta_Phi, nu3, delta, nu1, nu2, B;
      typ frequencies[2];
      int every = 10;
      int n_vertical = 342;
      typ diffusion_rate[n_vertical];
      typ p = 1.;
      typ q = 3.;
      char LibCenPath[800];
      char out_path[800];
      strcpy(out_path, "/home/atipique/Documents/git/K2138/Stability_map/");
      strcat(out_path, "fundamental_frequencies.txt");
      strcpy(LibCenPath, "/home/atipique/Documents/git/K2138/Stability_map/");
      strcat(LibCenPath, "LibCen_123_aptidal.txt");
      FILE * file = fopen(out_path, "w");
      fprintf(file, "The columns are delta, B, delta_Phi, dPhi, Phi_lc, Gamma, Upsilon, diffusion rate at Phi_lc, diffusion rate at Phi_lc + dPhi, diffusion rate at Phi_lc + 2*dPhi...\n\n");
      typ * data = NULL;
      int n_data;
      data = readFromFile(LibCenPath, &n_data);
      if (n_data % 23 != 0){
            fprintf(stderr, "Error : The file must have exactly 21 columns per line.\n");
            abort();
      }
      int n_lines = n_data/23;
      for (i = 0; i < n_lines; i ++){
            if (i % every == 0){
                  printf("i = %d\n", i);
                  /******** Getting the old coordinates at the libration centers ********/
                  lbd1 = fmod(data[23*i + 11], 2.*M_PI); g1 = fmod(data[23*i + 12], 2.*M_PI); Lbd1 = data[23*i + 13]; D1 = data[23*i + 14];
                  lbd2 = fmod(data[23*i + 15], 2.*M_PI); g2 = fmod(data[23*i + 16], 2.*M_PI); Lbd2 = data[23*i + 17]; D2 = data[23*i + 18];
                  lbd3 = fmod(data[23*i + 19], 2.*M_PI); g3 = fmod(data[23*i + 20], 2.*M_PI); Lbd3 = data[23*i + 21]; D3 = data[23*i + 22];
                  delta_Phi = data[23*i + 5];
                  nu3 = data[23*i + 10];
                  delta = data[23*i + 4];
                  B = data[23*i + 3];
                  Phi_lc = Lbd1/p;  Gamma = (p+q)*Lbd1/p + Lbd2;  Upsilon = Lbd1 + Lbd2 + Lbd3;
                  fprintf(file, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf", delta, B, delta_Phi, 1.32*Gamma/((typ) n_vertical)*delta_Phi, Phi_lc, Gamma, Upsilon);
                  /******** Getting the coordinates of the 1 dof model at the libration centers ********/
                  #pragma omp parallel for num_threads(19) private(Phi, Lbd1, Lbd2, Lbd3, X_old, frequencies, nu1, nu2) shared(diffusion_rate)
                  for (j = 0; j < n_vertical; j ++){
                        printf("    j = %d\n", j);
                        Phi  = Phi_lc + 1.32*Gamma*((typ) j)/((typ) n_vertical)*delta_Phi;
                        Lbd1 = p*Phi;
                        Lbd2 = Gamma - (p+q)*Phi;
                        Lbd3 = q*Phi - Gamma + Upsilon;
                        X_old[1] = lbd1; X_old[2]  = g1; X_old[3]  = Lbd1; X_old[4]  = D1;
                        X_old[5] = lbd2; X_old[6]  = g2; X_old[7]  = Lbd2; X_old[8]  = D2; 
                        X_old[9] = lbd3; X_old[10] = g3; X_old[11] = Lbd3; X_old[12] = D3;
                        FundamentalFrequency(0.0625, fabs(2.*M_PI/nu3*4.), X_old, 2, 1, 2, frequencies, 2);
                        nu1 = *frequencies;
                        nu2 = *(frequencies + 1);
                        diffusion_rate[j] = fabs((nu1 - nu2)/nu1);
                  }
                  for (j = 0; j < n_vertical; j ++){
                        fprintf(file, " %.16lf", diffusion_rate[j]);
                  }
                  fprintf(file, "\n");
            }
      }
      
      free(data);  data = NULL;
      fclose(file);

      /*typ frequencies[2];
      typ nu1, nu2;
      FundamentalFrequency(0.0625, 40000., X_old, 2, 1, 2, frequencies, 2);
      nu1 = *frequencies;
      nu2 = *(frequencies + 1);
      printf("nu1_1 = %.16lf, nu1_2 = %.16lf, diffusion index = %.6lf\n", nu1, nu2, log10(fabs((nu1-nu2)/nu1)));*/
      
      
      //SABAn(0.25, 10000000., 1000, X_old, 10);
      //SABAH1064(.0625, 50000., 40, X_old);
      
      //EquilibriumFind(X_old, 1);
      //SABAn(0.2,  50000., 1, X_old, 9);
      //LibrationCenterNAFF(X_old, 0.0078125, 100000., 2, 2, 40);
      //SABAn(0.2, 50000., 1, X_old, 8);
      //LibrationCenterFollow(X_old, epsilon/400000., 20000, 2);
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
      
      //LibrationCenterFind(X_old, 1);
      //SABAn(0.2, 50000., 1, X_old, 9);
      
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
