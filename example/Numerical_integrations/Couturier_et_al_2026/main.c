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
      transformation_display(); //Displaying the change of variable from the natural Poincaré's coordinates to the coordinates relevant from the resonance chain
      
      typ X[Nd*how_many_planet + 1];
      X_init(X);
      
      SABAn(.26, 4.e10, 524288, X, 10); //Integrating with a SABA10 and a timestep 0.26 for 4.e10 periods and one output every 524288 timestep. Should take a few days to terminate.
  
      deallocation();
      return 0;
}
