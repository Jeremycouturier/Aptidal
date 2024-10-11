#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"
#include "calculus.h"


typ X_old_t0[4*how_many_planet + 1];
typ X_new_t0[4*how_many_planet + 1];
typ X_uv_t0 [4*how_many_planet + 1];


typ fast_pow(typ x, int power){

      /******** Raises x to the integer power power where 0 <= power <= 5 ********/
      
      if (power == 0){
            return 1.;
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
      else if (power == 5){
            typ x2 = x*x;
            return x2*x2*x;
      }
      else{
            fprintf(stderr, "\nAptidal error : The power must be an integer between 0 and 5 in function fast_pow.\n");
            abort();
      }
}


void dHdold(typ * dH, typ * X_old, int KP){

      /******** Fills the array dH with the gradient of the Keplerian part of the Hamiltonian if KP is 0 ********/
      /******** and with the perturbative part if KP is 1. Indexes 4*i - 3 to 4*i contain dH/dlbd_i,     ********/
      /******** dH/d(-vrp_i), dH/dLbd_i and dH/dsq2DsLi=sqrt(2*D_i*Lbd_i0)*dH/dD_i, respectively.        ********/
      /******** The total gradient (of H_K+H_P) is computed if KP is 2                                   ********/
      
      int i, j, k, pi, pj, the_gcd, is_coorbital;
      typ p, q, l, r;
      int n, m;
      typ * Cppqk;
      typ C, angle, sq2DIsLI, sq2DJsLJ;
      typ lbd_I, lbd_J, g_I, g_J, Lbd_I, Lbd_J, D_I, D_J;
      typ mi, beta_i, mu_i, Lbd_i, mI, mJ, aJ, factor;
      typ Delta, Delta5, xi, dHdxi, dHdDelta, Dvp, dHdDvp, dHdsq2DIsLI, dHdsq2DJsLJ;
      
      /******** Initializing dH to 0 ********/
      for (i = 1; i <= 4*how_many_planet; i ++){
            *(dH + i) = 0.;
      }
      
      for (i = 1; i <= how_many_planet; i ++){
      
            if (KP == 0 || KP == 2){
                  /******** Retrieving mi, beta_i and mu_i ********/
                  mi     = masses[i];
                  mu_i   = G*(m0 + mi);
                  beta_i = m0*mi/(m0 + mi);
      
                  /******** When differenciating with respect to Lbd, then no need to loop over j (H_p does not depend on Lbd) ********/
                  Lbd_I       = X_old[4*i - 1];
                  dH[4*i - 1] = beta_i*beta_i*beta_i*mu_i*mu_i/(Lbd_I*Lbd_I*Lbd_I);
            }
            if (KP == 1 || KP == 2){
                  /******** Looping over j ********/
                  for (j = i + 1; j <= how_many_planet; j ++){
                        Cppqk = Cppq[i][j];
                  
                        /******** Getting the old coordinates ********/
                        aJ     = sma[j];
                        mI     = masses[i];
                        mJ     = masses[j];
                        factor = G*mI*mJ/aJ;
                        lbd_I  = X_old[4*i - 3];
                        lbd_J  = X_old[4*j - 3];
                        g_I    = X_old[4*i - 2]; // -vrp_I
                        g_J    = X_old[4*j - 2]; // -vrp_J
                        Lbd_I  = Lbd_0[i];       // Using nominal value since the perturbation is expanded at order 0 in semi-major axes
                        Lbd_J  = Lbd_0[j];       // Using nominal value since the perturbation is expanded at order 0 in semi-major axes
                        D_I    = X_old[4*i];
                        D_J    = X_old[4*j];
                  
                        /******** Getting sqrt(2*D_I/Lbd_I) and sqrt(2*D_J/Lbd_J) ********/
                        sq2DIsLI = sqrt(2.0*D_I/Lbd_I);
                        sq2DJsLJ = sqrt(2.0*D_J/Lbd_J);
                  
                        /******** Retrieving the value of p for a resonance p:p+q ********/
                        pi           = p_i[i];
                        pj           = p_i[j];
                        is_coorbital = (pi == pj && pi != 0) ? 1 : 0; 
                        if (pi == 0 || pj == 0){
                              p = 0.;
                        }
                        else{
                              the_gcd = gcd(pi, pj);
                              p       = (typ) (pi/the_gcd);
                        }
                  
                        /******** Looping over secular and non-co-orbital MMR resonances ********/
                        if (!is_coorbital){
                              for (k = 1; k < 32; k ++){
                                    C = Cppqk[k]; //Coefficient of the corresponding term in the Hamiltonian
                                    if (C != 0.){ //If the term actually exists
                                          q         = qnmlr[k][0];  n = (int) qnmlr[k][1];  m = (int) qnmlr[k][2];  l = qnmlr[k][3];  r = qnmlr[k][4];
                                          angle     = l*p*lbd_I - l*(p + q)*lbd_J - r*g_I - (l*q - r)*g_J; //The angle inside the cosine of that term.

                                          /******** Derivation with respect to lbd_i ********/
                                          *(dH + 4*i - 3) -= C*l*p      *fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle)*factor;
                                          /******** Derivation with respect to lbd_j ********/
                                          *(dH + 4*j - 3) += C*l*(p + q)*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle)*factor;

                                          /******** Derivation with respect to -vrp_i ********/
                                          *(dH + 4*i - 2) += C*r        *fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle)*factor;
                                          /******** Derivation with respect to -vrp_j ********/
                                          *(dH + 4*j - 2) += C*(l*q - r)*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*sin(angle)*factor;

                                          /******** Derivation with respect to sqrt(2*D_i/Lbd_i0) ********/
                                          if (m != 0){
                                                *(dH + 4*i) += ((typ) m)      *C*fast_pow(sq2DIsLI, m - 1)*fast_pow(sq2DJsLJ, n - m)*cos(angle)*factor;
                                          }
                                          /******** Derivation with respect to sqrt(2*D_j/Lbd_j0) ********/
                                          if(n - m != 0){
                                                *(dH + 4*j) += ((typ) (n - m))*C*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m - 1)*cos(angle)*factor;
                                          }    
                                    }
                              }
                        }

                        /******** The pair (i,j) is co-orbital ********/
                        else{
                              xi          = lbd_I - lbd_J;
                              Delta       = sqrt(2. - 2.*cos(xi));
                              Dvp         = g_J - g_I;
                              Delta5      = Delta*Delta*Delta*Delta*Delta;
                              dHdxi       = -sin(xi);
                              dHdxi      += max_deg >= 2 ? -((sq2DIsLI*sq2DJsLJ*(-(3.*sin(3.*xi-Dvp))-32.*sin(2.*xi-Dvp)-9.*sin(xi+Dvp)+26.*sin(xi-Dvp)))/(8.*Delta5))
                                                          -2.*sq2DIsLI*sq2DJsLJ*sin(2.*xi-Dvp)+((sq2DJsLJ*sq2DJsLJ/2.+sq2DIsLI*sq2DIsLI/2.)*(-(10.*sin(2.*xi))-8.*sin(xi)))/(4.*Delta5)
                                                          -(-(sq2DJsLJ*sq2DJsLJ/2.)-sq2DIsLI*sq2DIsLI/2.)*sin(xi) : 0.;
                              dHdDelta    = 1./(Delta*Delta);
                              dHdDelta   += max_deg >= 2 ? (5.*sq2DIsLI*sq2DJsLJ*(cos(3.*xi-Dvp)+16.*cos(2.*xi-Dvp)+9.*cos(xi+Dvp)-26.*cos(xi-Dvp)))/(8.*Delta5*Delta)-
                                                          (5.*(sq2DJsLJ*sq2DJsLJ/2.+sq2DIsLI*sq2DIsLI/2.)*(5.*cos(2.*xi)+8.*cos(xi)-13.))/(4.*Delta5*Delta) : 0.;
                              dHdxi      += sin(xi)/Delta*dHdDelta;
                              dHdDvp      = max_deg >= 2 ? sq2DIsLI*sq2DJsLJ*sin(2.*xi-Dvp)
                                                         -(sq2DIsLI*sq2DJsLJ*(sin(3.*xi-Dvp)+16.*sin(2.*xi-Dvp)-9.*sin(xi+Dvp)-26.*sin(xi-Dvp)))/(8.*Delta5) : 0.;
                              dHdsq2DIsLI = max_deg >= 2 ? -((sq2DJsLJ*(cos(3.*xi-Dvp)+16.*cos(2.*xi-Dvp)+9.*cos(xi+Dvp)-26.*cos(xi-Dvp)))/(8.*Delta5))+sq2DJsLJ*cos(2.*xi-Dvp)+
                                                             (sq2DIsLI*(5.*cos(2.*xi)+8.*cos(xi)-13.))/(4.*Delta5)-sq2DIsLI*cos(xi) : 0.;
                              dHdsq2DJsLJ = max_deg >= 2 ? -((sq2DIsLI*(cos(3.*xi-Dvp)+16.*cos(2.*xi-Dvp)+9.*cos(xi+Dvp)-26.*cos(xi-Dvp)))/(8.*Delta5))+sq2DIsLI*cos(2.*xi-Dvp)+
                                                             (sq2DJsLJ*(5.*cos(2.*xi)+8.*cos(xi)-13.))/(4.*Delta5)-sq2DJsLJ*cos(xi) : 0.;
                        
                              /******** Derivation with respect to lbd_i  ********/
                              *(dH + 4*i - 3) += dHdxi*factor;
                              /******** Derivation with respect to lbd_j  ********/
                              *(dH + 4*j - 3) -= dHdxi*factor;
                        
                              /******** Derivation with respect to -vrp_i  ********/
                              *(dH + 4*i - 2) -= dHdDvp*factor;
                              /******** Derivation with respect to -vrp_j  ********/
                              *(dH + 4*j - 2) += dHdDvp*factor;
                        
                              /******** Derivation with respect to sqrt(2*D_i/Lbd_i0)  ********/
                              *(dH + 4*i) += dHdsq2DIsLI*factor;
                              /******** Derivation with respect to sqrt(2*D_j/Lbd_j0)  ********/
                              *(dH + 4*j) += dHdsq2DJsLJ*factor;
                        }
                  }
            }
      }
}


void dHdnew(typ * dH_polar, typ * dH_rect, typ * dH_old, typ * X_new, typ * X_uv){

      /******** Fills arrays dH_polar and dH_rect with the gradient of the new polar and     ********/
      /******** rectangular variables. Indexes 4*i - 3 to 4*i of dH_polar contain dH/dphi_i, ********/
      /******** dH/dsig_i, dH/dPhi_i and dH/dsq2DsLi=sqrt(2*D_i*Lbd_i0)*dH/dD_i, while those ********/
      /******** of dH_rect contain dH/dphi_i, dH/dv_i, dH/dPhi_i and dH/du_i                 ********/
      /******** The array dH_old was initialized by function dHdold.                         ********/

      int i, k;
      int N = how_many_planet;
      typ factor, Dk, uuvv;

      /******** Initializing to 0 ********/
      for (i = 1; i <= 4*how_many_planet; i ++){
            *(dH_polar + i) = 0.;
            *(dH_rect  + i) = 0.;
      }
      
      /******** Polar case ********/
      for (k = 1; k <= N; k ++){
            Dk     = X_new[4*k];
            factor = sqrt(2.*Lbd_0[k]*Dk);
            for (i = 1; i <= N; i ++){
                  dH_polar[4*k - 3] +=        Transpose_inv [k]    [i]*dH_old[4*i - 3] + Transpose_inv [k]    [i + N]*dH_old[4*i - 2];
                  dH_polar[4*k - 2] +=        Transpose_inv [k + N][i]*dH_old[4*i - 3] + Transpose_inv [k + N][i + N]*dH_old[4*i - 2];
                  dH_polar[4*k - 1] +=        Transformation[k]    [i]*dH_old[4*i - 1] + Transformation[k]    [i + N]*dH_old[4*i]/factor;
                  dH_polar[4*k]     += factor*Transformation[k + N][i]*dH_old[4*i - 1] + Transformation[k + N][i + N]*dH_old[4*i];
            }
      }
      
      /******** Rectangular case ********/
      for (k = 1; k <= N; k ++){
            Dk               = X_new[4*k];
            uuvv             = X_uv[4*k - 2]*X_uv[4*k - 2] + X_uv[4*k]*X_uv[4*k];
            dH_rect[4*k - 3] = dH_polar [4*k - 3];
            dH_rect[4*k - 2] = sin(X_new[4*k - 2])*dH_polar[4*k] + X_uv[4*k]    /uuvv*dH_polar[4*k - 2];
            dH_rect[4*k - 1] = dH_polar [4*k - 1];
            dH_rect[4*k]     = cos(X_new[4*k - 2])*dH_polar[4*k] - X_uv[4*k - 2]/uuvv*dH_polar[4*k - 2];
      }
      
      /******** Verifying that the gradient of non-degrees of freedom is close to 0. To be removed when the code is robust ********/
      for (k = 1; k <= how_many_nondof; k ++){
            i = nondof[k];
            if (fabs(dH_rect[4*i - 3]) > 1.e-11){ //If dH/dphi_i != 0 where (phi_i; Phi_i) is not a degree of freedom
                  fprintf(stderr, "Error : dH/d_phi_%d is not 0 despite (phi_%d; Phi_%d) not being a degree of freedom. dH/d_phi_%d = %.13lf.\n", i, i, i, i, dH_rect[4*i - 3]);
                  abort();
            }
      }
}


void old2new(typ * X_old, typ * X_new, typ * X_uv){

      /******** Computes and stores inside X_new and X_uv the new coordinates           ********/
      /******** (polar and rectangular) from the knowlegde of the old coordinates X_old ********/
      /******** Indexes 4*i-3 to 4*i contain :                                          ********/
      /******** (lbd_i, -vrp_i; Lbd_i, D_i) for X_old                                   ********/
      /******** (phi_i,  sig_i; Phi_i, D_i) for X_new                                   ********/
      /******** (phi_i,    v_i; Phi_i, u_i) for X_uv                                    ********/

      int i, k;
      int N = how_many_planet;
      typ D;
      
      /******** Computing X_new ********/
      for (k = 1; k <= N; k ++){
            X_new[4*k - 3] = 0.;  X_new[4*k - 2] = 0.;  X_new[4*k - 1] = 0.;  X_new[4*k] = 0.;  D = 0.;
            for (i = 1; i <= N; i ++){
                  X_new[4*k - 3] += Transformation[k]    [i]*X_old[4*i - 3] + Transformation[k]    [i + N]*X_old[4*i - 2];
                  X_new[4*k - 2] += Transformation[k + N][i]*X_old[4*i - 3] + Transformation[k + N][i + N]*X_old[4*i - 2];
                  X_new[4*k - 1] += Transpose_inv [k]    [i]*X_old[4*i - 1] + Transpose_inv [k]    [i + N]*X_old[4*i];
                  D              += Transpose_inv [k + N][i]*X_old[4*i - 1] + Transpose_inv [k + N][i + N]*X_old[4*i];
            }
            if (fabs(D - X_old[4*k]) > 1.e-10){
                  fprintf(stderr, "Error : D_new = %.12lf significantly differs from D_old = %.12lf in function old2new.\n", D, X_old[4*k]);
                  abort();
            }
            X_new[4*k] = X_old[4*k]; //The D coordinate is unchanged
      }
      
      /******** Computing X_uv ********/
      for (k = 1; k <= N; k ++){
            X_uv[4*k - 3] = X_new[4*k - 3];
            X_uv[4*k - 2] = sqrt(2.*X_new[4*k]/Lbd_0[k])*sin(X_new[4*k - 2]);
            X_uv[4*k - 1] = X_new[4*k - 1];
            X_uv[4*k]     = sqrt(2.*X_new[4*k]/Lbd_0[k])*cos(X_new[4*k - 2]);
      }
}


void new2old(typ * X_old, typ * X_new, typ * X_uv){

      /******** Computes and stores inside X_new and X_old the new polar coordinates and the     ********/
      /******** old polar coordinates from the knowledge of the new rectangular coordinates X_uv ********/
      /******** Indexes 4*i-3 to 4*i contain :                                                   ********/
      /******** (lbd_i, -vrp_i; Lbd_i, D_i) for X_old                                            ********/
      /******** (phi_i,  sig_i; Phi_i, D_i) for X_new                                            ********/
      /******** (phi_i,    v_i; Phi_i, u_i) for X_uv                                             ********/

      int i, k;
      int N = how_many_planet;
      typ D;
      
      /******** Computing X_new ********/
      for (k = 1; k <= N; k ++){
            X_new[4*k - 3] = X_uv[4*k - 3];
            X_new[4*k - 2] = atan2(X_uv[4*k - 2], X_uv[4*k]);
            X_new[4*k - 1] = X_uv[4*k - 1];
            X_new[4*k]     = 0.5*Lbd_0[k]*(X_uv[4*k - 2]*X_uv[4*k - 2] + X_uv[4*k]*X_uv[4*k]);
      }
      
      /******** Computing X_old ********/
      for (k = 1; k <= N; k ++){
            X_old[4*k - 3] = 0.;  X_old[4*k - 2] = 0.;  X_old[4*k - 1] = 0.;  X_old[4*k] = 0.;  D = 0.;
            for (i = 1; i <= N; i ++){
                  X_old[4*k - 3] += Transpose_inv [i][k]    *X_new[4*i - 3] + Transpose_inv [i + N][k]    *X_new[4*i - 2];
                  X_old[4*k - 2] += Transpose_inv [i][k + N]*X_new[4*i - 3] + Transpose_inv [i + N][k + N]*X_new[4*i - 2];
                  X_old[4*k - 1] += Transformation[i][k]    *X_new[4*i - 1] + Transformation[i + N][k]    *X_new[4*i];
                  D              += Transformation[i][k + N]*X_new[4*i - 1] + Transformation[i + N][k + N]*X_new[4*i];
            }
            if (fabs(D - X_new[4*k]) > 1.e-10){
                  fprintf(stderr, "Error : D_old = %.12lf significantly differs from D_new = %.12lf in function new2old.\n", D, X_new[4*k]);
                  abort();
            }
            X_old[4*k] = X_new[4*k]; //The D coordinate is unchanged
      }
}


void X_old_init(typ * X_old){

      /******** Initializes the array X_old of the old variables such that indexes ********/
      /******** 4i-3 to 4i contains lbd_i, -vrp_i, Lbd_i and D_i, respectively     ********/

      int j;
      typ m, beta, mu, Lbd;

      for (j = 1; j <= how_many_planet; j ++){
            m    = masses[j];
            beta = m0*m/(m0 + m);
            mu   = G*(m0 + m);
            Lbd  = beta*sqrt(mu*sma[j]);
            X_old[4*j - 3] = lbd[j];
            X_old[4*j - 2] = -vrp[j];
            X_old[4*j - 1] = Lbd;
            X_old[4*j]     = max(Lbd*(1. - sqrt(1. - ecc[j]*ecc[j])), 1.e-17);
      }
}


void nonDofReinit(typ * X_new, typ * X_uv){

      /******** Reinitializes non-degrees-of-freedom to their initial values ********/

      int i, k;
      
      for (i = 1; i <= how_many_nondof; i ++){
            k = nondof[i];
            X_new[4*k - 3] = X_new_t0[4*k - 3];
            X_new[4*k - 1] = X_new_t0[4*k - 1];
            X_uv [4*k - 3] = X_uv_t0 [4*k - 3];
            X_uv [4*k - 1] = X_uv_t0 [4*k - 1];
      }
}


void exp_tau_LB_Ralston(typ tau, typ * X_old){

      /******** Performs a step of size tau on the perturbative part of the Hamiltonian. The operator exp(tau*L_B) is   ********/
      /******** approximated to (Id + tau*L_B + tau^2/2*L_B^2) by a Runge-Kutta 2 (Ralston) because B is not integrable ********/

      int i;
      typ two_third = 2./3.;
      typ dH_old  [4*how_many_planet + 1];
      typ dH_polar[4*how_many_planet + 1];
      typ dH_rect [4*how_many_planet + 1];
      typ K2      [4*how_many_planet + 1];
      typ X_new   [4*how_many_planet + 1];
      typ X_uv    [4*how_many_planet + 1];
      typ X_K2    [4*how_many_planet + 1];

      dHdold(dH_old, X_old, 1);                       //Computing gradient of the perturbative Hamiltonian
      old2new(X_old, X_new, X_uv);                    //(lbd, -vrp; Lbd, D) -> (phi, sig; Phi, D) -> (phi, v; Phi, u)
      dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv); //Gradient of the perturbation in the new coordinates. dH_rect is k1 of Ralston method
      for (i = 1; i <= how_many_planet; i ++){        //Computing the argument passed to the vector field to compute k2
            X_K2[4*i - 3] = X_uv[4*i - 3] + two_third*tau*dH_rect[4*i - 1];
            X_K2[4*i - 2] = X_uv[4*i - 2] + two_third*tau*dH_rect[4*i]    /Lbd_0[i];
            X_K2[4*i - 1] = X_uv[4*i - 1] - two_third*tau*dH_rect[4*i - 3];
            X_K2[4*i]     = X_uv[4*i]     - two_third*tau*dH_rect[4*i - 2]/Lbd_0[i];
      }
      new2old(X_old, X_new, X_K2);
      dHdold(dH_old, X_old, 1);
      dHdnew(dH_polar, K2, dH_old, X_new, X_K2);      //K2 is k2 of Ralston method
      for (i = 1; i <= how_many_planet; i ++){        //Completing Ralston method
            X_uv[4*i - 3] += 0.25*tau*(dH_rect[4*i - 1] + 3.*K2[4*i - 1]);
            X_uv[4*i - 2] += 0.25*tau*(dH_rect[4*i]     + 3.*K2[4*i])     /Lbd_0[i];
            X_uv[4*i - 1] -= 0.25*tau*(dH_rect[4*i - 3] + 3.*K2[4*i - 3]);
            X_uv[4*i]     -= 0.25*tau*(dH_rect[4*i - 2] + 3.*K2[4*i - 2]) /Lbd_0[i];
      }
      new2old(X_old, X_new, X_uv); //(phi, v; Phi, u) -> (phi, sig; Phi, D) -> (lbd, -vrp; Lbd, D)
}


void SABA1(typ tau, typ T, int output_step, typ * X_old){

      /******** Integrates the Hamiltonian with a SABA1 method         ********/
      /******** (Leapfrog) for a time T with a timestep tau. Outputs   ********/
      /******** every output_step timestep to file pth/SABA1.txt.      ********/
      /******** The Keplerian part A is integrated exactly, but the    ********/
      /******** perturbation B is integrated with a Ralston method     ********/
      /******** Since Ralston method is order 2, that means that in    ********/
      /******** exp(tau/2*L_A)*exp(tau*L_B)*exp(tau/2*L_A), the middle ********/
      /******** exponential is replaced by (Id+tau*L_B+tau^2/2*L_B^2)  ********/
      /******** This yields L_K=L_A+L_B-tau^2/24*L_{A,{A,B}}+tau^2/12* ********/
      /******** L_{{A,B},B}-tau^2/6*L_B^3. Therefore, the conservative ********/
      /******** error should be in O(epsilon*tau^2), whereas the non-  ********/
      /******** conservative error should be in O(epsilon^3*tau^2),    ********/
      /******** where tau is the timestep and epsilon the ratio        ********/
      /******** between the Keplerian part and the perturbative part.  ********/
      
      char file_path[800];
      FILE * file;
      int N_step, iter, i;
      typ e, a, sig, H;
      typ dH_old  [4*how_many_planet + 1];
      typ X_new   [4*how_many_planet + 1];
      typ X_uv    [4*how_many_planet + 1];
      typ old_buff[4*how_many_planet + 1];
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "SABA1.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file SABA1.txt in function SABA1.\n");
            abort();
      }
      
      /******** Integrating ********/
      N_step = (int) ceil(T/tau);
      for (iter = 0; iter < N_step; iter ++){
      
            /******** Writing to file ********/
            if (iter%output_step == 0){
                  for (i = 1; i <= how_many_planet; i ++){
                        old_buff[4*i - 3] = X_old[4*i - 3];
                        old_buff[4*i - 2] = X_old[4*i - 2];
                        old_buff[4*i - 1] = X_old[4*i - 1];
                        old_buff[4*i]     = X_old[4*i];
                        if (iter){ //Need to perform a step exp(-tau/2*L_A) on the buffer before outputting
                              old_buff[4*i - 3] -= 0.5*tau*dH_old[4*i - 1]; //lbd_i = lbd_i - tau/2*dH_K/dLbd_i
                        }
                  }
                  old2new(old_buff, X_new, X_uv);
                  /******** Outputting ********/
                  if (!iter){
                        fprintf(file, "Numerical integration of the Hamiltonian with a SABA1 integrator in the coordinates (phi, v; Phi, u) defined as: \n");
                        fprintf(file, "u_j = sqrt(2*D_j/Lbd_j0)*cos(sig_j) and v_j = sqrt(2*D_j/Lbd_j0)*sin(sig_j) where Lbd_j0 = m_j*sqrt(G*m_0*a_j0) and a_j0 is the initial\n");
                        fprintf(file, "semi-major axis given in the file parameters.h.\n");
                        fprintf(file, "\n");
                        fprintf(file, "The SABA1 integrator is a symplectic integrator defined by the operator (see Laskar & Robutel, 2001):\n");
                        fprintf(file, "exp(tau/2*L_A)*exp(tau*L_B)*exp(tau/2*L_A) where A is the Keplerian part of the Hamiltonian, B its perturbative part and tau is the timestep.\n");
                        fprintf(file, "L_A = {A, .} and L_B = {B, .} denote the Lie derivatives, that is, L_A(f) is df/dt along the trajectories of A.\n");
                        fprintf(file, "\n");
                        fprintf(file, "Since B is not integrable, it cannot be integrated exactly and the operator exp(tau*L_B) is unknown. Instead, B was integrated\n");
                        fprintf(file, "approximately with a Runge-Kutta method of order 2 (Ralston method), which is equivalenty to replacing the operator exp(tau*L_B) with the\n");
                        fprintf(file, "operator Id + tau*L_B + 1/2*tau^2*L_B^2. Writing exp(tau/2*L_A)*(Id + tau*L_B + 1/2*tau^2*L_B^2)*exp(tau/2*L_A) = exp(tau*L_K), we obtain the\n");
                        fprintf(file, "equivalent Lie derivative L_K = L_A + L_B - tau^2/24*L_{A,{A,B}} + tau^2/12*L_{{A,B},B} - tau^2/6*L_B^3 + O(tau^3).\n");
                        fprintf(file, "\n");
                        fprintf(file, "The term in L_B^3 cannot be put in Hamiltonian form and corresponds to the non-conservative error that is commited by the approximate\n");
                        fprintf(file, "integration of B. Since B is epsilon small with respect to A, this SABA1 has a non-conservative error of size O(tau^2*epsilon^3)\n");
                        if (how_many_planet == 2){
                              fprintf(file, "and a conservative error of size O(tau^2*espilon) where epsilon = (m_1 + m_2)/m_0.\n");
                        }
                        else{
                              fprintf(file, "and a conservative error of size O(tau^2*espilon) where epsilon = (m_1 + ... + m_%d)/m_0.\n", how_many_planet);
                        }
                        fprintf(file, "\n");
                        fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 2 + 7*how_many_planet, how_many_planet);
                        fprintf(file, "Time, Hamiltonian, phi_j, v_j, Phi_j, u_j, a_j, e_j, sig_j.\n");
                        fprintf(file, "\n");
                  }
                  H = Hamiltonian(old_buff);
                  fprintf(file, "%.12lf %.12lf", tau*(typ) iter, H);
                  for (i = 1; i <= how_many_planet; i ++){
                        sig = atan2(X_uv[4*i - 2], X_uv[4*i]);
                        e   = sqrt(1. - (1. - old_buff[4*i]/old_buff[4*i - 1])*(1. - old_buff[4*i]/old_buff[4*i - 1]));
                        a   = old_buff[4*i - 1]*old_buff[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                        fprintf(file, " %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", X_uv[4*i - 3], X_uv[4*i - 2], X_uv[4*i - 1], X_uv[4*i], a, e, sig);
                  }
                  fprintf(file, "\n");
            }
      
      
            /******** Step exp(tau/2*L_A). For the first iteration only ********/
            if (!iter){
                  dHdold(dH_old, X_old, 0); //Computing gradient of the Keplerian Hamiltonian
                  for (i = 1; i <= how_many_planet; i ++){
                        X_old[4*i - 3] += 0.5*tau*dH_old[4*i - 1]; //lbd_i = lbd_i + tau/2*dH_K/dLbd_i
                  }
            }
            
            
            /******** Step exp(tau*L_B). Approximated to (Id + tau*L_B + tau^2/2*L_B^2) by a Ralston method because B is not integrable ********/
            exp_tau_LB_Ralston(tau, X_old);
            /*dHdold(dH_old, X_old, 1);                       //Computing gradient of the perturbative Hamiltonian
            old2new(X_old, X_new, X_uv);                    //(lbd, -vrp; Lbd, D) -> (phi, sig; Phi, D) -> (phi, v; Phi, u)
            dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv); //Gradient of the perturbation in the new coordinates. dH_rect is k1 of Ralston method
            for (i = 1; i <= how_many_planet; i ++){        //Computing the argument passed to the vector field to compute k2
                  X_K2[4*i - 3] = X_uv[4*i - 3] + two_third*tau*dH_rect[4*i - 1];
                  X_K2[4*i - 2] = X_uv[4*i - 2] + two_third*tau*dH_rect[4*i]    /Lbd_0[i];
                  X_K2[4*i - 1] = X_uv[4*i - 1] - two_third*tau*dH_rect[4*i - 3];
                  X_K2[4*i]     = X_uv[4*i]     - two_third*tau*dH_rect[4*i - 2]/Lbd_0[i];
            }
            new2old(X_old, X_new, X_K2);
            dHdold(dH_old, X_old, 1);
            dHdnew(dH_polar, K2, dH_old, X_new, X_K2);      //K2 is k2 of Ralston method
            for (i = 1; i <= how_many_planet; i ++){        //Completing Ralston method
                  X_uv[4*i - 3] += 0.25*tau*(dH_rect[4*i - 1] + 3.*K2[4*i - 1]);
                  X_uv[4*i - 2] += 0.25*tau*(dH_rect[4*i]     + 3.*K2[4*i])     /Lbd_0[i];
                  X_uv[4*i - 1] -= 0.25*tau*(dH_rect[4*i - 3] + 3.*K2[4*i - 3]);
                  X_uv[4*i]     -= 0.25*tau*(dH_rect[4*i - 2] + 3.*K2[4*i - 2]) /Lbd_0[i];
            }
            new2old(X_old, X_new, X_uv); //(phi, v; Phi, u) -> (phi, sig; Phi, D) -> (lbd, -vrp; Lbd, D)*/
            
            /******** Step exp(tau*L_A) ********/
            dHdold(dH_old, X_old, 0);    //Computing gradient of the Keplerian Hamiltonian
            for (i = 1; i <= how_many_planet; i ++){
                  X_old[4*i - 3] += tau*dH_old[4*i - 1]; //lbd_i = lbd_i + tau*dH_K/dLbd_i
            }
      }

      /******** Closing output file ********/
      fclose(file);
}


void SABA1_average(typ tau, typ T, int Hanning_order, typ * X_uv_mean, typ * X_old){

      /******** Same as above but computes and stores the average of X_uv ********/
      /******** instead of writing to file. Hanning_order is the order of ********/
      /******** the Hanning filter used to compute the averages. The time ********/
      /******** step is tau/2 to use a Simpson 3pt method to compute      ********/
      /******** integrals. Fills X_uv_mean with the average of (phi_j,    ********/
      /******** v_j; Phi_j, u_j)                                          ********/
      
      int N_step, iter, i;
      typ dH_old  [4*how_many_planet + 1];
      typ X_new   [4*how_many_planet + 1];
      typ X_uv_0  [4*how_many_planet + 1]; //Left   of the segment
      typ X_uv_1  [4*how_many_planet + 1]; //Middle of the segment
      typ X_uv_2  [4*how_many_planet + 1]; //Right  of the segment
      typ old_buff[4*how_many_planet + 1];
      typ facto   [11] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3628800.};
      typ two_to_p[6]  = {1., 2., 4., 8., 16., 32.};
      typ Cp, Hanning0, Hanning1, Hanning2, t0, t1, t2;
      
      if (Hanning_order > 5){
            fprintf(stderr, "Error: The order of the Hanning filter cannot exceed 5.\n");
            abort();
      }
      Cp = two_to_p[Hanning_order]*facto[Hanning_order]*facto[Hanning_order]/facto[2*Hanning_order];
      
      /******** Initializing averages to 0 ********/
      for (i = 1; i <= 4*how_many_planet; i ++){
            X_uv_mean[i] = 0.;
      }
      
      /******** Integrating ********/
      tau   /= 2.;
      N_step = (int) ceil(2.*T/tau);
      N_step+= N_step%2 == 0 ? 0 : 1; //Making sure that the integration lasts for an even number of steps
      T      = 0.5*tau* (typ) N_step;
      for (iter = 1; iter <= N_step; iter ++){
      
            /******** Step exp(tau/2*L_A). For the first iteration only ********/
            if (iter == 1){
                  old2new(X_old, X_new, X_uv_0);
                  dHdold(dH_old, X_old, 0); //Computing gradient of the Keplerian Hamiltonian
                  for (i = 1; i <= how_many_planet; i ++){
                        X_old[4*i - 3] += 0.5*tau*dH_old[4*i - 1]; //lbd_i = lbd_i + tau/2*dH_K/dLbd_i
                  }
            }
            
            
            /******** Step exp(tau*L_B). Approximated to (Id + tau*L_B + tau^2/2*L_B^2) by a Ralston method because B is not integrable ********/
            exp_tau_LB_Ralston(tau, X_old);
            /*dHdold(dH_old, X_old, 1);                       //Computing gradient of the perturbative Hamiltonian
            old2new(X_old, X_new, X_uv);                    //(lbd, -vrp; Lbd, D) -> (phi, sig; Phi, D) -> (phi, v; Phi, u)
            dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv); //Gradient of the perturbation in the new coordinates. dH_rect is k1 of Ralston method
            for (i = 1; i <= how_many_planet; i ++){        //Computing the argument passed to the vector field to compute k2
                  X_K2[4*i - 3] = X_uv[4*i - 3] + two_third*tau*dH_rect[4*i - 1];
                  X_K2[4*i - 2] = X_uv[4*i - 2] + two_third*tau*dH_rect[4*i]    /Lbd_0[i];
                  X_K2[4*i - 1] = X_uv[4*i - 1] - two_third*tau*dH_rect[4*i - 3];
                  X_K2[4*i]     = X_uv[4*i]     - two_third*tau*dH_rect[4*i - 2]/Lbd_0[i];
            }
            new2old(X_old, X_new, X_K2);
            dHdold(dH_old, X_old, 1);
            dHdnew(dH_polar, K2, dH_old, X_new, X_K2); //K2 is k2 of Ralston method
            for (i = 1; i <= how_many_planet; i ++){   //Completing Ralston method
                  X_uv[4*i - 3] += 0.25*tau*(dH_rect[4*i - 1] + 3.*K2[4*i - 1]);
                  X_uv[4*i - 2] += 0.25*tau*(dH_rect[4*i]     + 3.*K2[4*i])     /Lbd_0[i];
                  X_uv[4*i - 1] -= 0.25*tau*(dH_rect[4*i - 3] + 3.*K2[4*i - 3]);
                  X_uv[4*i]     -= 0.25*tau*(dH_rect[4*i - 2] + 3.*K2[4*i - 2]) /Lbd_0[i];
            }
            new2old(X_old, X_new, X_uv); //(phi, v; Phi, u) -> (phi, sig; Phi, D) -> (lbd, -vrp; Lbd, D)*/
            
            
            /******** Step exp(tau*L_A) ********/
            dHdold(dH_old, X_old, 0);    //Computing gradient of the Keplerian Hamiltonian
            for (i = 1; i <= how_many_planet; i ++){
                  X_old[4*i - 3] += tau*dH_old[4*i - 1]; //lbd_i = lbd_i + tau*dH_K/dLbd_i
            }
            
            
            /******** Averaging ********/
            for (i = 1; i <= how_many_planet; i ++){
                        old_buff[4*i - 3] = X_old[4*i - 3] - 0.5*tau*dH_old[4*i - 1];
                        old_buff[4*i - 2] = X_old[4*i - 2];
                        old_buff[4*i - 1] = X_old[4*i - 1];
                        old_buff[4*i]     = X_old[4*i];
            }
            if (iter % 2 == 0){
                  old2new(old_buff, X_new, X_uv_2);
                  t0       = tau* (typ) (iter - 2);
                  t1       = tau* (typ) (iter - 1);
                  t2       = tau* (typ) iter;
                  Hanning0 = Cp*fast_pow(1. - cos(M_PI*t0/T), Hanning_order);
                  Hanning1 = Cp*fast_pow(1. - cos(M_PI*t1/T), Hanning_order);
                  Hanning2 = Cp*fast_pow(1. - cos(M_PI*t2/T), Hanning_order);
                  for (i = 1; i <= how_many_planet; i ++){
                        X_uv_mean[4*i - 3] += tau*(Hanning0*X_uv_0[4*i - 3] + 4.*Hanning1*X_uv_1[4*i - 3] + Hanning2*X_uv_2[4*i - 3])/(6.*T);
                        X_uv_mean[4*i - 2] += tau*(Hanning0*X_uv_0[4*i - 2] + 4.*Hanning1*X_uv_1[4*i - 2] + Hanning2*X_uv_2[4*i - 2])/(6.*T);
                        X_uv_mean[4*i - 1] += tau*(Hanning0*X_uv_0[4*i - 1] + 4.*Hanning1*X_uv_1[4*i - 1] + Hanning2*X_uv_2[4*i - 1])/(6.*T);
                        X_uv_mean[4*i]     += tau*(Hanning0*X_uv_0[4*i]     + 4.*Hanning1*X_uv_1[4*i]     + Hanning2*X_uv_2[4*i])    /(6.*T);
                  }
                  for (i = 1; i <= 4*how_many_planet; i ++){
                        X_uv_0[i] = X_uv_2[i];
                  }
            }
            else{
                  old2new(old_buff, X_new, X_uv_1);
            }
      }
}


void RK2(typ tau, typ T, int output_step){

      /******** Integrates the Hamiltonian with a Runge-Kutta 2        ********/
      /******** method (Ralston) for a time T with a timestep tau      ********/
      /******** Outputs every output_step timestep to file pth/RK2.txt ********/
      /******** The non-conservative error is in O(tau^2)              ********/
      
      char file_path[800];
      FILE * file;
      int N_step, iter, i;
      typ two_third = 2./3.;
      typ e, sig, parenthesis, H;
      typ dH_old  [4*how_many_planet + 1];
      typ dH_polar[4*how_many_planet + 1];
      typ dH_rect [4*how_many_planet + 1];
      typ K2      [4*how_many_planet + 1];
      typ X_old   [4*how_many_planet + 1];
      typ X_new   [4*how_many_planet + 1];
      typ X_uv    [4*how_many_planet + 1];
      typ X_K2    [4*how_many_planet + 1];
      X_old_init(X_old);
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "RK2.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file SABA1.txt in function RK2.\n");
            abort();
      }
      
      /******** Integrating ********/
      N_step = (int) ceil(T/tau);
      for (iter = 0; iter < N_step; iter ++){
      
            /******** Writing to file ********/
            if (iter%output_step == 0){
            
                  if (iter){
                        new2old(X_old, X_new, X_uv);
                  }
                  else{
                        old2new(X_old, X_new, X_uv);
                        nonDofReinit(X_new, X_uv);
                        fprintf(file, "Numerical integration of the Hamiltonian with a RK2 Ralston integrator in the coordinates (phi, v; Phi, u) defined as: \n");
                        fprintf(file, "u_j = sqrt(2*D_j/Lbd_j0)*cos(sig_j) and v_j = sqrt(2*D_j/Lbd_j0)*sin(sig_j) where Lbd_j0 = m_j*sqrt(G*m_0*a_j0) and a_j0 is the initial\n");
                        fprintf(file, "semi-major axis given in the file parameters.h.\n");
                        fprintf(file, "\n");
                        fprintf(file, "Ralston method is non-symplectic and has a non-conservative error of size O(tau^2) where tau is the timestep.\n");
                        fprintf(file, "\n");
                        fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 2 + 6*how_many_planet, how_many_planet);
                        fprintf(file, "Time, Hamiltonian, phi_j, v_j, Phi_j, u_j, e_j, sig_j.\n");
                        fprintf(file, "\n");
                  }
                  H = Hamiltonian(X_old);
                  fprintf(file, "%.12lf %.12lf", tau*(typ) iter, H);
                  for (i = 1; i <= how_many_planet; i ++){
                        sig         = atan2(X_uv[4*i - 2], X_uv[4*i]);
                        parenthesis = 1. - 0.5*Lbd_0[i]/X_old[4*i - 1]*(X_uv[4*i - 2]*X_uv[4*i - 2] + X_uv[4*i]*X_uv[4*i]);
                        e           = sqrt(1. - parenthesis*parenthesis);
                        fprintf(file, " %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", X_uv[4*i - 3], X_uv[4*i - 2], X_uv[4*i - 1], X_uv[4*i], e, sig);
                  }
                  fprintf(file, "\n");
            }
            
            if (iter){
                  new2old(X_old, X_new, X_uv);
            }
            dHdold(dH_old, X_old, 2);                       //Computing gradient of the perturbative Hamiltonian
            old2new(X_old, X_new, X_uv);
            nonDofReinit(X_new, X_uv);
            dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv); //Gradient of the perturbation in the new coordinates. dH_rect is k1 of Ralston method
            for (i = 1; i <= how_many_planet; i ++){        //Computing the argument passed to the vector field to compute k2
                  X_K2[4*i - 3] = X_uv[4*i - 3] + two_third*tau*dH_rect[4*i - 1];
                  X_K2[4*i - 2] = X_uv[4*i - 2] + two_third*tau*dH_rect[4*i]    /Lbd_0[i];
                  X_K2[4*i - 1] = X_uv[4*i - 1] - two_third*tau*dH_rect[4*i - 3];
                  X_K2[4*i]     = X_uv[4*i]     - two_third*tau*dH_rect[4*i - 2]/Lbd_0[i];
            }
            nonDofReinit(X_new, X_K2);                      //Reinitializing non-degrees-of-freedom
            new2old(X_old, X_new, X_K2);
            nonDofReinit(X_new, X_K2);                      //Reinitializing non-degrees-of-freedom
            dHdold(dH_old, X_old, 2);
            dHdnew(dH_polar, K2, dH_old, X_new, X_K2);      //K2 is k2 of Ralston method
            for (i = 1; i <= how_many_planet; i ++){        //Completing Ralston method
                  X_uv[4*i - 3] += 0.25*tau*(dH_rect[4*i - 1] + 3.*K2[4*i - 1]);
                  X_uv[4*i - 2] += 0.25*tau*(dH_rect[4*i]     + 3.*K2[4*i])     /Lbd_0[i];
                  X_uv[4*i - 1] -= 0.25*tau*(dH_rect[4*i - 3] + 3.*K2[4*i - 3]);
                  X_uv[4*i]     -= 0.25*tau*(dH_rect[4*i - 2] + 3.*K2[4*i - 2]) /Lbd_0[i];
            }
            nonDofReinit(X_new, X_uv);                      //Reinitializing non-degrees-of-freedom
      }

      /******** Closing output file ********/
      fclose(file);
}


void RK4(typ tau, typ T, int output_step){

      /******** Integrates the Hamiltonian with a Runge-Kutta 4  ********/
      /******** method for a time T with a timestep tau. Outputs ********/
      /******** every output_step timestep to file pth/RK4.txt   ********/
      /******** The non-conservative error is in O(tau^4)        ********/
      
      char file_path[800];
      FILE * file;
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "RK4.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file RK4.txt in function RK4.\n");
            abort();
      }
      
      /******** Closing output file ********/
      fclose(file);
}


typ Hamiltonian(typ * X_old){

      /******** Returns the value of the Hamiltonian ********/

      typ H = 0.;
      int i, j, k, pi, pj, the_gcd, is_coorbital;
      typ p, q, l, r;
      int n, m;
      typ * Cppqk;
      typ C, angle, sq2DIsLI, sq2DJsLJ;
      typ lbd_I, lbd_J, g_I, g_J, Lbd_I, Lbd_J, D_I, D_J;
      typ mi, beta_i, mu_i, Lbd_i, mI, mJ, aJ, factor;
      typ Delta, Delta5, xi, Dvp;
      
      
      for (i = 1; i <= how_many_planet; i ++){
            /******** Keplerian part ********/
            mi     = masses[i];
            mu_i   = G*(m0 + mi);
            beta_i = m0*mi/(m0 + mi);
            Lbd_i  = X_old[4*i - 1];
            H     -= beta_i*beta_i*beta_i*mu_i*mu_i/(2.*Lbd_i*Lbd_i);
            /******** Perturbative part ********/
            for (j = i + 1; j <= how_many_planet; j ++){
                  Cppqk        = Cppq[i][j];
                  aJ           = sma[j];
                  mI           = masses[i];
                  mJ           = masses[j];
                  factor       = G*mI*mJ/aJ;
                  lbd_I        = X_old[4*i - 3];
                  lbd_J        = X_old[4*j - 3];
                  g_I          = X_old[4*i - 2]; // -vrp_I
                  g_J          = X_old[4*j - 2]; // -vrp_J
                  Lbd_I        = Lbd_0[i];       // Using nominal value since the perturbation is expanded at order 0 in semi-major axes
                  Lbd_J        = Lbd_0[j];       // Using nominal value since the perturbation is expanded at order 0 in semi-major axes
                  D_I          = X_old[4*i];
                  D_J          = X_old[4*j];
                  sq2DIsLI     = sqrt(2.0*D_I/Lbd_I);
                  sq2DJsLJ     = sqrt(2.0*D_J/Lbd_J);
                  pi           = p_i[i];
                  pj           = p_i[j];
                  is_coorbital = (pi == pj && pi != 0) ? 1 : 0; 
                  if (pi == 0 || pj == 0){
                        p = 0.;
                  }
                  else{
                        the_gcd = gcd(pi, pj);
                        p       = (typ) (pi/the_gcd);
                  }
                  
                  if (!is_coorbital){
                        for (k = 1; k < 32; k ++){
                              C = Cppqk[k]; //Coefficient of the corresponding term in the Hamiltonian
                              if (C != 0.){ //If the term actually exists
                                    q         = qnmlr[k][0];  n = (int) qnmlr[k][1];  m = (int) qnmlr[k][2];  l = qnmlr[k][3];  r = qnmlr[k][4];
                                    angle     = l*p*lbd_I - l*(p + q)*lbd_J - r*g_I - (l*q - r)*g_J; //The angle inside the cosine of that term.
                                    H        += C*fast_pow(sq2DIsLI, m)*fast_pow(sq2DJsLJ, n - m)*cos(angle)*factor;
                              }
                        }
                  }
                  else{
                        xi     = lbd_I - lbd_J;
                        Delta  = sqrt(2. - 2.*cos(xi));
                        Dvp    = g_J - g_I;
                        Delta5 = Delta*Delta*Delta*Delta*Delta;
                        H     += factor*(cos(xi) - 1./Delta);
                        if (max_deg >= 2){
                              H += factor*((D_I/Lbd_I + D_J/Lbd_J)*(5.*cos(2.*xi) + 8.*cos(xi) - 13.)/(4.*Delta5));
                              H -= factor*((D_I/Lbd_I + D_J/Lbd_J)*cos(xi));
                              H += factor*(sq2DIsLI*sq2DJsLJ*cos(2.*xi - Dvp));
                              H -= factor*(sq2DIsLI*sq2DJsLJ/(8.*Delta5)*(cos(3.*xi - Dvp) + 16.*cos(2.*xi - Dvp) - 26.*cos(xi - Dvp) + 9.*cos(xi + Dvp)));
                        }
                  }
            }
      }
      return H;
}


void PointPrint(typ * X_old, int iter){

      /******** Verifies the conservation of the first integrals and displays the point in terminal ********/

      int N = how_many_planet;
      int i, k;
      int spaces = (int) floor(log10((typ) iter));
      typ a, e, parenthesis;
      typ X_new[4*how_many_planet + 1];
      typ X_uv [4*how_many_planet + 1];
      
      /******** Verifying conservation of the first integrals. To be removed when the code is robust ********/
      old2new(X_old, X_new, X_uv);
      for (i = 1; i <= how_many_nondof; i ++){
            k = nondof[i];
            if (fabs(X_uv[4*k - 1] - X_uv_t0[4*k - 1]) > 1.e-9){
                  fprintf(stderr, "Error : Phi_%d is not conserved in function PointPrint.\n", k);
            }
      }

      printf("Iteration n° %d : ", iter);
      
      /******** Printing the lbd_i ********/
      printf("(");
      for (i = 1; i <= N - 1; i ++){
            printf("lbd_%d, ", i);
      }
      printf("lbd_%d) = (", N);
      for (i = 1; i <= N - 1; i ++){
            printf("%.12lf, ", fmod(X_old[4*i - 3], 2.*M_PI));
      }
      printf("%.12lf)\n                 ", fmod(X_old[4*N - 3], 2.*M_PI));
      for (i = 0; i < spaces; i ++){printf(" ");}
      
      /******** Printing the vrp_i ********/
      printf("(");
      for (i = 1; i <= N - 1; i ++){
            printf("vrp_%d, ", i);
      }
      printf("vrp_%d) = (", N);
      for (i = 1; i <= N - 1; i ++){
            printf("%.12lf, ", fmod(-X_old[4*i - 2], 2.*M_PI));
      }
      printf("%.12lf)\n                 ", fmod(-X_old[4*N - 2], 2.*M_PI));
      for (i = 0; i < spaces; i ++){printf(" ");}
      
      /******** Printing the a_i ********/
      printf("(");
      for (i = 1; i <= N - 1; i ++){
            printf("a_%d, ", i);
      }
      printf("a_%d) ", N); for (i = 1; i <= N; i ++){printf("  ");} printf("= (");
      for (i = 1; i <= N - 1; i ++){
            a = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            printf("%.12lf, ", a);
      }
      a = X_old[4*N - 1]*X_old[4*N - 1]*(m0 + masses[N])/(G*m0*m0*masses[N]*masses[N]);
      printf("%.12lf)\n                 ", a);
      for (i = 0; i < spaces; i ++){printf(" ");}
      
      /******** Printing the e_i ********/
      printf("(");
      for (i = 1; i <= N - 1; i ++){
            printf("e_%d, ", i);
      }
      printf("e_%d) ", N); for (i = 1; i <= N; i ++){printf("  ");} printf("= (");
      for (i = 1; i <= N - 1; i ++){
            parenthesis = 1. - X_old[4*i]/X_old[4*i - 1];
            e           = sqrt(1. - parenthesis*parenthesis);
            printf("%.12lf, ", e);
      }
      parenthesis = 1. - X_old[4*N]/X_old[4*N - 1];
      e           = sqrt(1. - parenthesis*parenthesis);
      printf("%.12lf)\n\n", e);
}


void EquilibriumFind(typ * X_old, int verbose){

      /******** Finds a fixed point of the Hamiltonian by iteratively integrating ********/
      /******** and keeping only the average of the Fourier decomposition.        ********/
      /******** verbose is a boolean that specifies if the function should talk   ********/

      int n_iter = 1;
      int i;
      typ precision = 1.;
      typ X_new[4*how_many_planet + 1];
      typ X_uv [4*how_many_planet + 1];
      typ xvXu [4*how_many_planet + 1];
      typ T   = 8000.;
      typ tau = 0.25;
      old2new(X_old, X_new, X_uv);
      
      if (verbose){
            printf("-------------------------------------------------------------------------------------------\n\n");
            printf("Starting the search for a fixed point.\n\n");
            PointPrint(X_old, 0);
      }
      
      /******** I first try to converge with a low precision using a small integration time and a large timestep ********/
      while(precision > 1.e-2 && n_iter <= 8){
            SABA1_average(2.*tau, T/2., 1, xvXu, X_old);
            nonDofReinit(X_new, xvXu);
            precision = 0.;
            for (i = 1; i <= how_many_planet; i ++){
                  precision += fabs(xvXu[4*i - 3] - X_uv[4*i - 3]) + fabs(xvXu[4*i - 2] - X_uv[4*i - 2]) + fabs(xvXu[4*i - 1] - X_uv[4*i - 1]) + fabs(xvXu[4*i] - X_uv[4*i]);
            }
            precision /= (typ) how_many_planet;   
            for (i = 1; i <= how_many_planet; i ++){
                  X_uv[4*i - 3] = xvXu[4*i - 3];
                  X_uv[4*i - 2] = xvXu[4*i - 2];
                  X_uv[4*i - 1] = xvXu[4*i - 1];
                  X_uv[4*i]     = xvXu[4*i];
            } 
            new2old(X_old, X_new, X_uv);
            PointPrint(X_old, n_iter);
            n_iter ++;
      }
      
      if (n_iter > 8 && precision > 1.e-2){ //This initial condition is hopeless
            if (verbose){
                  printf("This initial condition is hopeless.\n");
            }
            return;
      }
      if (verbose){
            printf("A precision of 10^-2 was reached.\n");
      }
      
      /******** I now refine to a moderate precision using a moderate integration time and a moderate timestep ********/
      while(precision > 1.e-6){
            SABA1_average(tau, T, 5, xvXu, X_old);
            nonDofReinit(X_new, xvXu);
            precision = 0.;
            for (i = 1; i <= how_many_planet; i ++){
                  precision += fabs(xvXu[4*i - 3] - X_uv[4*i - 3]) + fabs(xvXu[4*i - 2] - X_uv[4*i - 2]) + fabs(xvXu[4*i - 1] - X_uv[4*i - 1]) + fabs(xvXu[4*i] - X_uv[4*i]);
            }
            precision /= (typ) how_many_planet;   
            for (i = 1; i <= how_many_planet; i ++){
                  X_uv[4*i - 3] = xvXu[4*i - 3];
                  X_uv[4*i - 2] = xvXu[4*i - 2];
                  X_uv[4*i - 1] = xvXu[4*i - 1];
                  X_uv[4*i]     = xvXu[4*i];
            } 
            new2old(X_old, X_new, X_uv);
            PointPrint(X_old, n_iter);
            n_iter ++;
            if (n_iter > 20 && precision > 1.e-6){
                  fprintf(stderr, "Error : In function EquilibriumFind, the precision cannot reach 10^-6 even though it reached 10^-3. Try decreasing tau and increasing T.\n");
                  abort();
            }
      }
      if (verbose){
            printf("A precision of 10^-6 was reached.\n");
      }
      
      /******** I now refine to a high precision using a large integration time and a small timestep ********/
      while(precision > 1.e-10){
            SABA1_average(tau/4., 4*T, 5, xvXu, X_old);
            nonDofReinit(X_new, xvXu);
            precision = 0.;
            for (i = 1; i <= how_many_planet; i ++){
                  precision += fabs(xvXu[4*i - 3] - X_uv[4*i - 3]) + fabs(xvXu[4*i - 2] - X_uv[4*i - 2]) + fabs(xvXu[4*i - 1] - X_uv[4*i - 1]) + fabs(xvXu[4*i] - X_uv[4*i]);
            }
            precision /= (typ) how_many_planet;   
            for (i = 1; i <= how_many_planet; i ++){
                  X_uv[4*i - 3] = xvXu[4*i - 3];
                  X_uv[4*i - 2] = xvXu[4*i - 2];
                  X_uv[4*i - 1] = xvXu[4*i - 1];
                  X_uv[4*i]     = xvXu[4*i];
            } 
            new2old(X_old, X_new, X_uv);
            PointPrint(X_old, n_iter);
            n_iter ++;
            if (n_iter > 30 && precision > 1.e-10){
                  fprintf(stderr, "Error : In function EquilibriumFind, the precision cannot reach 10^-8 even though it reached 10^-6. Try decreasing tau and increasing T.\n");
                  abort();
            }
      }
      if (verbose){
            printf("A precision of 10^-10 was reached.\n");
      }
      
      SABA1(tau, 2.*T, 1 + (int) (2.*T/tau/32768.), X_old);
}
