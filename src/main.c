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
      UnaveragedSABAn(0.0625, 7500., 1, X_old, 4);
      /*X_old_init(X_old);
      for (int _ = 0; _ < 20; _ ++){
            LibrationCenterGradientDescent(X_old, 0.125, 40000., 0.00000018, 4, 1);
      }
      UnaveragedSABAn(0.125, 40000., 1, X_old, 6);*/

  
      deallocation();
      return 0;
}
