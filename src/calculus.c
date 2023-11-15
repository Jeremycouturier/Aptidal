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
typ grad_old[4 * how_many_planet + 1];
typ * grad_new;
typ * grad_uv;
typ X_old   [4 * how_many_planet + 1];
typ * X_new;
typ * X_uv;



typ fast_pow(typ x, int power){

      /******** Raises x to the integer power power where 0 <= power <= 4 ********/
      
      if (power == 0){
            return 1.0;
      }
      else if (power == 1){
            return x;
      }
      else if (power == 2){
            return x*x;
      }
      else if (power == 3){
            return x*x*x;
      }
      else if (power == 4){
            typ x2 = x*x;
            return x2*x2;
      }
      else{
            fprintf(stderr, "\nAptidal error : The power must be an integer between 0 and 4 in function fast_pow.\n");
            abort();
      }
}

typ dH_old(int K){

      /******** Returns the partial derivative dH/dX_old[K] where 1 <= K <= 4 * how_many_planet ********/
      
      if (K < 1 || K > 4 * how_many_planet){
            fprintf(stderr, "\nAptidal error : K must be between 1 and 4 * how_many_planet in function dH_old.\n");
            abort();
      }
      
      int i, j, k, pi, pj, the_gcd, I, J;
      typ p, q, l, r;
      int n, m;
      typ * Cppqk;
      typ C, angle, sq2DIsLI, sq2DJsLJ;
      typ to_be_returned = 0.0;
      typ lbd_I, lbd_J, g_I, g_J, Lbd_I, Lbd_J, D_I, D_J;
      typ mi, beta_i, mu_i, Lbd_i, aI, aJ, mI, mJ, Di;
      
      /******** Getting the value of i ********/
      if (K % 4 == 0){
            i = K/4;
      }
      else{
            i = K/4 + 1;
      }
      
      /******** Retrieving mi, beta_i and mu_i ********/
      mi     = masses[i];
      mu_i   = G*(m0 + mi);
      beta_i = m0*mi/(m0 + mi);
      
      /******** If differenciating with respect to Lbd_i, then no need to loop over j ********/
      if (K % 4 == 3){
            Lbd_i = X_old[4*i - 1];
            return beta_i*beta_i*beta_i*mu_i*mu_i/(Lbd_i*Lbd_i*Lbd_i);
      }
      
      /******** If differenciating with respect to D_i and D_i is zero, we return 0.0 ********/
      Di = X_old[4*i];
      if (K % 4 == 0 && absolute(Di) < 1.0e-14){
            return 0.0;   
      }
      
      /******** Looping over j ********/
      for (j = 1; j <= how_many_planet; j++){
            if (i != j){
                  I     = min(i,j);
                  J     = max(i,j);
                  Cppqk = Cppq[I][J];
                  
                  /******** Getting the old coordinates ********/
                  aI    = sma[I];
                  aJ    = sma[J];
                  mI    = masses[I];
                  mJ    = masses[J];
                  lbd_I = X_old[4*I - 3];
                  lbd_J = X_old[4*J - 3];
                  g_I   = X_old[4*I - 2];   // -vrp_I
                  g_J   = X_old[4*J - 2];   // -vrp_J
                  Lbd_I = mI*sqrt(G*m0*aI); // Using value at t=0 since the perturbation is expanded at order 0 in semi-major axes
                  Lbd_J = mJ*sqrt(G*m0*aJ); // Using value at t=0 since the perturbation is expanded at order 0 in semi-major axes
                  D_I   = X_old[4*I];
                  D_J   = X_old[4*J];
                  
                  /******** Getting sqrt(2*D_I/Lbd_I) and sqrt(2*D_J/Lbd_J) ********/
                  sq2DIsLI = sqrt(2.0*D_I/Lbd_I);
                  sq2DJsLJ = sqrt(2.0*D_J/Lbd_J);
                  
                  /******** Retrieving the value of p for a resonance p:p+q ********/
                  pi      = p_i[I];
                  pj      = p_i[J];
                  if (pi == 0 && pj == 0){
                        p = 0.0;
                  }
                  else{
                        the_gcd = gcd(pi, pj);
                        p       = (typ) (pi / the_gcd);
                  }
                  
                  /******** Looping over secular and non-co-orbital MMR resonances ********/
                  for (k = 1; k < 32; k++){
                        C = Cppqk[k];  //Coefficient of the corresponding term in the Hamiltonian
                        if (C != 0.0){ //If the term actually exists
                              q         = qnmlr[k][0];  n = (int) qnmlr[k][1];  m = (int) qnmlr[k][2];  l = qnmlr[k][3];  r = qnmlr[k][4];
                              angle     = l*p*lbd_I - l*(p + q)*lbd_J - r*g_I - (l*q - r)*g_J; //The angle inside the cosine of that term.
                              
                              if (K % 4 == 1){      //Derivation with respect to lbd_i
                                    if (i == I){
                                          to_be_returned -= C*l*p*fast_pow(sq2DIsLI,m)*fast_pow(sq2DJsLJ,n - m)*sin(angle);
                                    }
                                    else if(i == J){
                                          to_be_returned += C*l*(p + q)*fast_pow(sq2DIsLI,m)*fast_pow(sq2DJsLJ,n - m)*sin(angle);
                                    }
                              }
                              else if (K % 4 == 2){ //Derivation with respect to -vrp_i
                                    if (i == I){
                                          to_be_returned += C*r*fast_pow(sq2DIsLI,m)*fast_pow(sq2DJsLJ,n - m)*sin(angle);
                                    }
                                    else if(i == J){
                                          to_be_returned += C*(l*q - r)*fast_pow(sq2DIsLI,m)*fast_pow(sq2DJsLJ,n - m)*sin(angle);
                                    }
                              }
                              else if (K % 4 == 0){ //Derivation with respect to D_i
                                    if (i == I){
                                          to_be_returned += ((typ) m)*C*fast_pow(sq2DIsLI,m)*fast_pow(sq2DJsLJ,n - m)*cos(angle)/(2.0*D_I);
                                    }
                                    else if(i == J){
                                          to_be_returned += ((typ) (n - m))*C*fast_pow(sq2DIsLI,m)*fast_pow(sq2DJsLJ,n - m)*cos(angle)/(2.0*D_J);
                                    }
                              }    
                        }
                  } 
            }
      }
      return to_be_returned;
}


void gradH_old(){

      /******** Computes the gradient of the Hamiltonian in the old variables (lbd, -vrp; Lbd, D) ********/
      /******** and stores it in the array grad_old, in such a way that the indexes 4*i to 4*i+3  ********/
      /******** contain dH/dlbd_i, dH/d(-vrp_i), dH/dLbd_i and dH/dD_i, respectively.             ********/

}
