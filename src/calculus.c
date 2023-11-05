#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"
#include "calculus.h"



/******** Defining some global variables ********/
typ grad_old[4 * how_many_planet + 4];
typ * grad_new;
typ * grad_uv;
typ X_old   [4 * how_many_planet + 4];
typ * X_new;
typ * X_uv;


void dH_old(int k){


}


void gradH_old(){

      /******** Computes the gradient of the Hamiltonian in the old variables (lbd, -vrp; Lbd, D) ********/
      /******** and stores it in the array grad_old, in such a way that the indexes 4*i to 4*i+3  ********/
      /******** contain dH/dlbd_i, dH/d(-vrp_i), dH/dLbd_i and dH/dD_i, respectively.             ********/

}
