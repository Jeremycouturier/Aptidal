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
      
      typ X[Nd*how_many_planet + 1];
      X_init(X);
      
      SABAH1064(3.518764, 438307200., 10380, X); //Integrating with a SABAH1064 and a timestep 3.518764 days for 1200000 years (438307200 days) and one output every ~100 years. Should take ~15 minutes
      //SABAH1064(-3.518764, -438307200., 10380, X); //To go in the past
  
      deallocation();
      return 0;
}
