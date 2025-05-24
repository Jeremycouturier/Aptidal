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
      typ X_old   [4*how_many_planet + 1];
      typ X_new   [4*how_many_planet + 1];
      typ X_uv    [4*how_many_planet + 1];
      typ X_cart  [4*how_many_planet + 1];
      typ dH_old  [4*how_many_planet + 1];
      typ dH_polar[4*how_many_planet + 1];
      typ dH_rect [4*how_many_planet + 1];
      typ epsilon = 0.;
      
      for (i = 1; i <= how_many_planet; i ++){
            epsilon += masses[i]/m0;
      }
      
      X_old_init(X_old);
      //SABA1(0.25, 2000., 1, X_old);
      //SABAn(0.25, 4000., 2, X_old, 4);
      //RK2(0.25, 10000., 1);
      EquilibriumFind(X_old);
      //LibrationCenterFind(X_old, 0);
      LibrationCenterFollow(X_old, -0.000005*epsilon, 2000, 0);
      //LibrationCenterFollow(X_old, -0.0000014*epsilon, 1000, 0);
      //PeriodicOrbitFind(X_old);
      //X_old_init(X_old);
      //SABAn(0.25, 400000., 20, X_old, 4);
      //UnaveragedSABAn(0.25/32., 400000., 100, X_old, 4);
      
      /*typ a, e, nu, mu, M, vp, l;
      a = 1.29;
      e = 0.203;
      M = 2.091;
      vp = 1.203;
      mu = G*(1. + 0.00002033);
      
      printf("a, e, M, vp = %.16lf, %.16lf, %.16lf, %.16lf\n", a, e, M, vp);
      
      nu = mean2true(M, mu, a, e);
      
      typ cart[4];
      typ alkhqp[6];
      
      ell2cart(a, e, 0., nu, vp, 0., mu, cart);
      cart2ell(cart, alkhqp, mu);
      
      a  = alkhqp[0];
      l  = alkhqp[1];
      vp = atan2(alkhqp[3],alkhqp[2]);
      M  = l - vp;
      e  = sqrt(alkhqp[3]*alkhqp[3] + alkhqp[2]*alkhqp[2]);
      
      printf("a, e, M, vp = %.16lf, %.16lf, %.16lf, %.16lf\n", a, e, M + 2.*M_PI, vp);*/
      
      /*typ a, e, E, mu, M, vp, l;
      a = 1.29;
      e = 0.203;
      M = 2.091;
      vp = 1.203;
      mu = G*(1. + 0.00002033);
      
      printf("a, e, M, vp = %.16lf, %.16lf, %.16lf, %.16lf\n", a, e, M, vp);
      
      E = mean2eccentric(M + vp, e*cos(vp), e*sin(vp));
      
      typ cart[4];
      typ alkhqp[6];
      
      ell2cart(a, e, 0., E, vp, 0., mu, cart);
      cart2ell(cart, alkhqp, mu);
      
      a  = alkhqp[0];
      l  = alkhqp[1];
      vp = atan2(alkhqp[3],alkhqp[2]);
      M  = l - vp;
      e  = sqrt(alkhqp[3]*alkhqp[3] + alkhqp[2]*alkhqp[2]);
      
      printf("a, e, M, vp = %.16lf, %.16lf, %.16lf, %.16lf\n", a, e, M + 2.*M_PI, vp);*/
      
      /*X_old_init(X_old);
      typ mu, a, e, vp, M, nu;
      
      for (i = 1; i <= how_many_planet; i ++){
            mu = G*(m0 + masses[i]);
            a  = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            e  = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
            vp = -X_old[4*i - 2];
            M  =  X_old[4*i - 3] - vp;
            nu = mean2true(M, mu, a, e);
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
      }
      
      for (i = 1; i <= how_many_planet; i ++){
            printf("(x_%d, y_%d, vx_%d, vy_%d) = (%.14lf, %.14lf, %.14lf, %.14lf)\n", i, i, i, i, X_cart[4*i - 3], X_cart[4*i - 2], X_cart[4*i - 1], X_cart[4*i]);
      }*/
      
      /*int i;
      typ a, e, sig, vp, M, mu, nu, beta, H;
      typ X_cart[4*how_many_planet + 1];
      typ alkhqp[6];

      for (i = 1; i <= how_many_planet; i ++){
            mu = G*(m0 + masses[i]);
            a  = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            e  = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
            vp = -X_old[4*i - 2];
            M  =  X_old[4*i - 3] - vp;
            nu = mean2true(M, mu, a, e);
            printf("M = %.12lf, nu = %.12lf\n", M, nu);
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
            //printf("X_old = (%.12lf, %.12lf, %.12lf, %.12lf)\n", X_old[4*i - 3], X_old[4*i - 2], X_old[4*i - 1], X_old[4*i]);
      }
      
      for (i = 1; i <= how_many_planet; i ++){
            mu = G*(m0 + masses[i]);
            kepsaut(X_cart + 4*i - 4, mu, 0.00035672);
      }
      
      for (i = 1; i <= how_many_planet; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
            beta = m0*masses[i]/(m0 + masses[i]);
            mu   = G*(m0 + masses[i]);
            cart2ell(X_cart + 4*i - 4, alkhqp, mu);
            e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
            X_old[4*i - 3] = alkhqp[1];
            printf("lbd_%d = %.12lf\n", i, alkhqp[1]);
            X_old[4*i - 2] = -atan2(alkhqp[3], alkhqp[2]);
            X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
            X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
            //printf("X_old = (%.12lf, %.12lf, %.12lf, %.12lf)\n", X_old[4*i - 3], X_old[4*i - 2], X_old[4*i - 1], X_old[4*i]);
      }
      old2new(X_old, X_new, X_uv);             // (lbd_j, -vrp_j; Lbd_j, D_j) -> (phi_j, v_j; Phi_j, u_j)
      H = 0.;
      printf("%.12lf %.12lf", 0., H);
      for (i = 1; i <= how_many_planet; i ++){
            sig = atan2(X_uv[4*i - 2], X_uv[4*i]);
            e   = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
            a   = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            printf(" %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", X_uv[4*i - 3], X_uv[4*i - 2], X_uv[4*i - 1], X_uv[4*i], a, e, sig);
      }
      printf("\n");*/
      
      
      /*X_old_init(X_old);
      dHdold(dH_old, X_old, 0);
      printf("dHkdLbd1 = %.14lf\n", dH_old[4*1 - 1]);
      printf("dHkdLbd2 = %.14lf\n", dH_old[4*2 - 1]);
      printf("dHkdLbd3 = %.14lf\n", dH_old[4*3 - 1]);
      dHdold(dH_old, X_old, 1);
      printf("dHpdlbd1 = %.23lf\n", dH_old[4*1 - 3]);
      printf("dHpdlbd2 = %.23lf\n", dH_old[4*2 - 3]);
      printf("dHpdlbd3 = %.23lf\n", dH_old[4*3 - 3]);
      printf("dHpdvrp1 = %.23lf\n", dH_old[4*1 - 2]);
      printf("dHpdvrp2 = %.23lf\n", dH_old[4*2 - 2]);
      printf("dHpdvrp3 = %.23lf\n", dH_old[4*3 - 2]);
      printf("dHpdsq2D1= %.23lf\n", dH_old[4*1]);
      printf("dHpdsq2D2= %.23lf\n", dH_old[4*2]);
      printf("dHpdsq2D3= %.23lf\n\n", dH_old[4*3]);
      
      dHdold(dH_old, X_old, 0);
      old2new(X_old, X_new, X_uv);
      dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv);
      printf("dHkdPhi1 = %.18lf\n", dH_rect[4*1 - 1]);
      printf("dHkdPhi2 = %.18lf\n", dH_rect[4*2 - 1]);
      printf("dHkdPhi3 = %.18lf\n", dH_rect[4*3 - 1]);
      printf("dHkdu1   = %.23lf\n", dH_rect[4*1]);
      printf("dHkdu2   = %.23lf\n", dH_rect[4*2]);
      printf("dHkdu3   = %.23lf\n", dH_rect[4*3]);
      printf("dHkdv1   = %.23lf\n", dH_rect[4*1 - 2]);
      printf("dHkdv2   = %.23lf\n", dH_rect[4*2 - 2]);
      printf("dHkdv3   = %.23lf\n", dH_rect[4*3 - 2]);
      printf("dHkdphi1 = %.18lf\n", dH_rect[4*1 - 3]);
      printf("dHkdphi2 = %.18lf\n", dH_rect[4*2 - 3]);
      printf("dHkdphi3 = %.18lf\n", dH_rect[4*3 - 3]);
      dHdold(dH_old, X_old, 1);
      dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv);
      printf("dHpdphi1 = %.18lf\n", dH_rect[4*1 - 3]);
      printf("dHpdu1   = %.23lf\n", dH_rect[4*1]);
      printf("dHpdu2   = %.23lf\n", dH_rect[4*2]);
      printf("dHpdu3   = %.23lf\n", dH_rect[4*3]);
      printf("dHpdv1   = %.23lf\n", dH_rect[4*1 - 2]);
      printf("dHpdv2   = %.23lf\n", dH_rect[4*2 - 2]);
      printf("dHpdv3   = %.23lf\n", dH_rect[4*3 - 2]);
      printf("dHpdphi2 = %.18lf\n", dH_rect[4*2 - 3]);
      printf("dHpdphi3 = %.18lf\n", dH_rect[4*3 - 3]);
      printf("dHpdPhi1 = %.18lf\n", dH_rect[4*1 - 1]);
      printf("dHpdPhi2 = %.18lf\n", dH_rect[4*2 - 1]);
      printf("dHpdPhi3 = %.18lf\n", dH_rect[4*3 - 1]);*/
      
      /*X_old_init(X_old);
      dHdold(dH_old, X_old, 0);
      printf("dHkdLbd1 = %.14lf\n", dH_old[4*1 - 1]);
      printf("dHkdLbd2 = %.14lf\n", dH_old[4*2 - 1]);
      dHdold(dH_old, X_old, 1);
      printf("dHpdlbd1 = %.14lf\n", dH_old[4*1 - 3]);
      printf("dHpdlbd2 = %.14lf\n", dH_old[4*2 - 3]);
      printf("dHpdvrp1 = %.14lf\n", dH_old[4*1 - 2]);
      printf("dHpdvrp2 = %.14lf\n", dH_old[4*2 - 2]);
      printf("dHpdsq2D1= %.14lf\n", dH_old[4*1]);
      printf("dHpdsq2D2= %.14lf\n", dH_old[4*2]);
      
      dHdold(dH_old, X_old, 0);
      old2new(X_old, X_new, X_uv);
      dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv);
      printf("dHkdPhi1 = %.18lf\n", dH_rect[4*1 - 1]);
      printf("dHkdPhi2 = %.18lf\n", dH_rect[4*2 - 1]);
      printf("dHkdu1   = %.26lf\n", dH_rect[4*1]);
      printf("dHkdu2   = %.26lf\n", dH_rect[4*2]);
      printf("dHkdv1   = %.26lf\n", dH_rect[4*1 - 2]);
      printf("dHkdv2   = %.26lf\n", dH_rect[4*2 - 2]);
      printf("dHkdphi1 = %.26lf\n", dH_rect[4*1 - 3]);
      printf("dHkdphi2 = %.26lf\n", dH_rect[4*2 - 3]);
      dHdold(dH_old, X_old, 1);
      dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv);
      printf("dHpdphi1 = %.18lf\n", dH_rect[4*1 - 3]);
      printf("dHpdu1   = %.26lf\n", dH_rect[4*1]);
      printf("dHpdu2   = %.26lf\n", dH_rect[4*2]);
      printf("dHpdv1   = %.26lf\n", dH_rect[4*1 - 2]);
      printf("dHpdv2   = %.26lf\n", dH_rect[4*2 - 2]);
      printf("dHpdphi2 = %.18lf\n", dH_rect[4*2 - 3]);
      printf("dHpdPhi1 = %.18lf\n", dH_rect[4*1 - 1]);
      printf("dHpdPhi2 = %.18lf\n", dH_rect[4*2 - 1]);*/
      
      /*
      typ X_old[4*how_many_planet + 1];
      typ X_new[4*how_many_planet + 1];
      typ X_uv [4*how_many_planet + 1];
      X_old_init(X_old);
      typ dH_K_old  [4*how_many_planet + 1];
      typ dH_P_old  [4*how_many_planet + 1];
      typ dH_K_polar[4*how_many_planet + 1];
      typ dH_P_polar[4*how_many_planet + 1];
      typ dH_K_rect [4*how_many_planet + 1];
      typ dH_P_rect [4*how_many_planet + 1];
      dHdold(dH_K_old, X_old, 0);
      dHdold(dH_P_old, X_old, 1);
      old2new(X_old, X_new, X_uv);
      dHdnew(dH_K_polar, dH_K_rect, dH_K_old, X_new, X_uv);
      dHdnew(dH_P_polar, dH_P_rect, dH_P_old, X_new, X_uv);

      int j;
      printf("************* Keplerian part *************\n\n");
      for (j = 1; j <= how_many_planet; j ++){
            printf("dHK/dphi_%d = %.12lf\n",   j, dH_K_rect[4*j - 3]);
            printf("dHK/dv_%d   = %.12lf\n",   j, dH_K_rect[4*j - 2]);
            printf("dHK/dPhi_%d = %.12lf\n",   j, dH_K_rect[4*j - 1]);
            printf("dHK/du_%d   = %.12lf\n\n", j, dH_K_rect[4*j]);
      }
      printf("************* Perturbative part *************\n\n");
      for (j = 1; j <= how_many_planet; j ++){
            printf("dHP/dphi_%d = %.12lf\n",   j, dH_P_rect[4*j - 3]);
            printf("dHP/dv_%d   = %.12lf\n",   j, dH_P_rect[4*j - 2]);
            printf("dHP/dPhi_%d = %.12lf\n",   j, dH_P_rect[4*j - 1]);
            printf("dHP/du_%d   = %.12lf\n\n", j, dH_P_rect[4*j]);
      }
      
      typ alpha = 0.43;
      struct pairOfReal b_12_01 = b_12(alpha, 1.0e-15);
      struct pairOfReal b_32_01 = b_k2(alpha, b_12_01, 3);
      struct pairOfReal b_52_01 = b_k2(alpha, b_32_01, 5);
      struct pairOfReal b_72_01 = b_k2(alpha, b_52_01, 7);
      
      typ C1 = -0.5*b_12_01.fst;                                               //Constant term in the inequality (0,0)
      typ C2 = -3.0/8.0*alpha*b_32_01.fst+(0.25+0.25*alpha*alpha)*b_32_01.snd; //Term propto XXb'+XbX' in the inequality (0,0)
      typ C3 = -0.125*alpha*b_32_01.snd;                                       //Term propto XXb+X'Xb' in the inequality (0,0)
      typ C13 = -9.0/128.0*alpha*alpha*b_52_01.fst+3.0/64.0*alpha*(1.0-alpha*alpha)*b_52_01.snd; //Term propto X^2Xb^2
      typ Cbeaucoup = 5.0/384.0*alpha*(1.0-alpha*alpha*alpha*alpha)*b_72_01.snd+55.0/512.0*alpha*alpha*alpha*b_72_01.snd
      -35.0/768.0*alpha*alpha*b_72_01.fst-85.0/768.0*alpha*alpha*alpha*alpha*b_72_01.fst; //Term propto X^3Xb^3
      
      printf("C1 = %.15lf\nC2 = %.15lf\nC3 = %.15lf\nC13 = %.15lf\nCbeaucoup = %.15lf\n", C1, C2, C3, C13, Cbeaucoup);*/
      
      /*typ alpha = 0.894;
      resonances[1][2](alpha,masses[1]);
      resonance_00(alpha);
      printf("C_00_1 = %.14lf\nC_00_2 = %.14lf\nC_00_3 = %.14lf\nC_00_4 = %.14lf\nC_00_5 = %.14lf\nC_00_6 = %.14lf\nC_00_7 = %.14lf\nC_00_8 = %.14lf\nC_00_9 = %.14lf\nC_00_10 = %.14lf\n"
      ,C_00_1, C_00_2, C_00_3, C_00_4, C_00_5, C_00_6, C_00_7, C_00_8, C_00_9, C_00_10);*/
      /*resonance_57(alpha,masses[1]);
      printf("C_ppp2_1 = %.14lf\nC_ppp2_2 = %.14lf\nC_ppp2_3 = %.14lf\n" ,C_ppp2_1, C_ppp2_2, C_ppp2_3);*/
      /*resonance_58(alpha,masses[1]);
      printf("C_ppp3_1 = %.14lf\nC_ppp3_2 = %.14lf\nC_ppp3_3 = %.14lf\nC_ppp3_4 = %.14lf\n" ,C_ppp3_1, C_ppp3_2, C_ppp3_3, C_ppp3_4);*/
      
      /*int i,j;
      for (i = 1; i < how_many_resonant + 1; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  if (i != j){
                        printf("k_%d%d = %d\n", i, j, k_ij[i][j]);
                  }
            }
      }
      for (i = 1; i < how_many_resonant; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  if (j > i){
                        printf("q_%d%d = %d\n", i, j, q_ij[i][j]);
                  }
            }
      }
      for (i = 1; i < how_many_resonant; i++){
            printf("k_%d = %d\n", i, k_i[i]);
            printf("k_star_%d = %d\n", i, k_star_i[i]);
            printf("q_%d = %d\n", i, q_i[i]);
      }
      
      for (i = 1; i < how_many_resonant + 1; i++){
            printf("subchain[%d] = %d\n", i, subchain[i]);
      }
      for (i = 1; i < how_many_missed + 1; i++){
            printf("chain_missed[%d] = %d\n", i, chain_missed[i]);
      }
      for (i = 1; i < how_many_resonant + 1; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  if (i != j){
                        printf("k_%d%d = ",i,j);
                        ratprint(rat_k_ij[i][j]);
                        printf("       ");
                  }
            }
            printf("\n");
      }
      for (i = 1; i < how_many_resonant - 1; i++){
            for (j = 1; j < how_many_resonant + 1; j++){
                  printf("l_%d%d = ",i,j);
                  ratprint(rat_l_ij[i][j]);
                  printf("       ");
            }
            printf("\n");
      }
      int r,s;
      for (i = 1; i < how_many_resonant - 1; i++){
            for (r = 1; r <= i; r++){
                  for (s = r + 1; s < how_many_resonant + 1; s++){
                        printf("NoverD_%d%d^(%d) = ", r, s, i);
                        ratprint(NoverD[i][r][s]);
                        printf("     ");
                  }
            }
            printf("\n");
      }
      for (i = 1; i < how_many_resonant - 1; i++){
            printf("c_%d = ", i);
            ratprint(rat_c_i[i]);
            printf("     ");
      }
      printf("\n");
      for (i = 1; i < 2*how_many_planet + 1; i++){
            for (j = 1; j < 2*how_many_planet + 1; j++){
                  printf("M_%d,%d=", i, j);
                  ratprint(transformation[i][j]);
                  printf("    ");
            }
            printf("\n");
      }
      printf("\n");
      for (i = 1; i < 2*how_many_planet + 1; i++){
            for (j = 1; j < 2*how_many_planet + 1; j++){
                  printf("Q_%d,%d=", i, j);
                  ratprint(transpose_inv[i][j]);
                  printf("    ");
            }
            printf("\n");
      }*/
      
      
      /*printf("\nPert = \n");
      printf("       %.14lf * ", C_ppp1_1); printf(string1);printf("\n");
      printf("       %.14lf * ", C_ppp1_2); printf(string2);printf("\n");
      printf("       %.14lf * ", C_ppp1_3); printf(string3);printf("\n");
      printf("       %.14lf * ", C_ppp1_4); printf(string4);printf("\n");
      printf("       %.14lf * ", C_ppp1_5); printf(string5);printf("\n");
      printf("       %.14lf * ", C_ppp1_6); printf(string6);printf("\n");
      printf("       %.14lf * ", C_ppp1_7); printf(string7);printf("\n");
      printf("       %.14lf * ", C_ppp1_8); printf(string8);printf("\n");
      printf("       %.14lf * ", C_ppp1_9); printf(string9);printf("\n");
      printf("       %.14lf * ", C_ppp1_10); printf(string10);printf("\n");
      printf("       %.14lf * ", C_ppp1_11); printf(string11);printf("\n");
      printf("       %.14lf * ", C_ppp1_12); printf(string12);printf("\n");
      printf("       %.14lf * ", C_ppp1_13); printf(string13);printf("\n");
      printf("       %.14lf * ", C_ppp1_14); printf(string14);printf("\n");
      printf("       %.14lf * ", C_ppp1_15); printf(string15);printf("\n");*/
      
      
      /******** GSL test ********/
      /*double a_data[] = { 1.0, 0.6, 0.0,
                         0.0, 1.5, 1.0,
                         0.0, 1.0, 1.0 };

      double inva[9];
     
      gsl_matrix_view m
          = gsl_matrix_view_array(a_data, 3, 3);
      gsl_matrix_view inv
          = gsl_matrix_view_array(inva,3,3);
      gsl_permutation * p = gsl_permutation_alloc (3);
 
      printf("The matrix is\n");
      for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                  printf(j==2?"%6.3f\n":"%6.3f ", gsl_matrix_get(&m.matrix,i,j));
 
      gsl_linalg_LU_decomp (&m.matrix, p, &s);    
      gsl_linalg_LU_invert (&m.matrix, p, &inv.matrix);
 
      printf("The inverse is\n");
      for (i = 0; i < 3; ++i)
            for (j = 0; j < 3; ++j)
                  printf(j==2?"%6.3f\n":"%6.3f ",gsl_matrix_get(&inv.matrix,i,j));
 
      gsl_permutation_free (p);*/


      deallocation();
      return 0;

}
