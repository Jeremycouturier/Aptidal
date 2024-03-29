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

void dH_old(int K, typ * dHk){

      /******** Returns into dHk the four partial derivatives dH/dX_old[j] for 4*K - 3 <= j <= 4*K ********/
      /******** If j = 4*K, then returns sqrt(2*D_j*Lbd_j)*dH/dD_j instead of dH/dD_j              ********/
      /******** The array dHk is assumed to be given initialized to 0                              ********/
      
      if (K < 1 || K > how_many_planet){
            fprintf(stderr, "\nAptidal error : K must be between 1 and how_many_planet in function dH_old.\n");
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
      i = K;
      
      /******** Retrieving mi, beta_i and mu_i ********/
      mi     = masses[i];
      mu_i   = G*(m0 + mi);
      beta_i = m0*mi/(m0 + mi);
      
      /******** When differenciating with respect to Lbd, then no need to loop over j (H_p does not depend on Lbd) ********/
      Lbd_I = X_old[4*K - 1];
      dHk[2] = beta_i*beta_i*beta_i*mu_i*mu_i/(Lbd_i*Lbd_i*Lbd_i);
      
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
                  Lbd_I = mI*sqrt(G*m0*aI); // Using value at t = 0 since the perturbation is expanded at order 0 in semi-major axes
                  Lbd_J = mJ*sqrt(G*m0*aJ); // Using value at t = 0 since the perturbation is expanded at order 0 in semi-major axes
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

                              /******** Derivation with respect to lbd_i ********/
                              if (i == I){
                                    *dHk -= C*l*p*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle);
                              }
                              else if(i == J){
                                    *dHk += C*l*(p + q)*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle);
                              }

                              /******** Derivation with respect to -vrp_i ********/
                              if (i == I){
                                    dHk[1] += C*r*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle);
                              }
                              else if(i == J){
                                    dHk[1] += C*(l*q - r)*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle);
                              }

                              /******** Derivation with respect to D_i ********/
                              if (i == I && m != 0){
                                    dHk[3] += ((typ) m)*C*fast_pow(sq2DIsLI, m - 1)*fast_pow(sq2DJsLJ, n - m)*cos(angle);
                              }
                              else if(i == J && n - m != 0){
                                    dHk[3] += ((typ) (n - m))*C*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m - 1)*cos(angle);
                              }    
                        }
                  }

                  /******** Looping over co-orbital MMR resonances ********/
            }
      }
}


void gradH_old(){

      /******** Computes the gradient of the Hamiltonian in the old variables (lbd, -vrp; Lbd, D) ********/
      /******** and stores it in the array grad_old, in such a way that the indexes 4*i to 4*i+3  ********/
      /******** contain dH/dlbd_i, dH/d(-vrp_i), dH/dLbd_i and dH/dD_i, respectively.             ********/

}


typ average_test(typ A, typ dt, typ T){

      /******** Test function to numerically compute the average A of the quasi-periodic function ********/
      /******** f(x) = cos(x) + 3*sin(sqrt(2)*x) + A. To be removed later                         ********/

      int j;
      int N   = (int) (T/dt);
      typ * f = (typ *)malloc(N*sizeof(typ));
      typ t, ft, to_be_returned;
      
      for (j = 0; j < N; j++){
            t  = dt * (typ) j;
            ft = cos(t) + 3.0*sin(sqrt(2.0)*t) + A;
            *(f + j) = ft;
      }
      
      
      
      free(f);
      f = NULL;
      return to_be_returned;
}





