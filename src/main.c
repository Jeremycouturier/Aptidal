#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
//#include <gsl/gsl_linalg.h>
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"
#include "calculus.h"
#include "intpla.h"


int main(){

      init();
            
      int i;
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
      //SABAn(0.25, 50000., 1, X_old, 10);
      LibrationCenterFollow(X_old, epsilon/400000., 20000, 2);
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
