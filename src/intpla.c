#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include "parameters.h"
#include "structure.h"
#include "coefficients.h"
#include "transformation.h"
#include "calculus.h"
#include "intpla.h"

typ nu_fast = 0.;
typ nu_reso = 0.;
typ avgs[1 + 4*how_many_planet];

void ell2cart(typ a, typ e, typ i, typ nu, typ varpi, typ Omega, typ mu, typ * cart){

      /******** Updates the array [X,Y,vX,vY] of the cartesian coordinates. mu is G*M, a is the       ********/
      /******** semi-major axis, e is the eccentricity, i is the inclination, nu is the true anomaly, ********/
      /******** omega is the longitude of periapsis and Omega is the longitude of the ascending node. ********/
      
      typ X,Y,vX,vY;                     //Cartesian coordinates
      typ X_buff,Y_buff,vX_buff,vY_buff; //Buffer for cartesian coordinates
      typ r;                             //Body's distance to Earth's center
      typ g;                             //Angular momentum per unit mass
      typ dnudt;
      typ drdt;
      typ cosnu    = cos(nu);
      typ sinnu    = sin(nu);
      typ cosvarpi = cos(varpi);
      typ sinvarpi = sin(varpi);
      typ q        = sin(i/2.)*cos(Omega);
      typ p        = sin(i/2.)*sin(Omega);
      typ chi      = cos(i/2.);
      typ pp       = 1. - 2.*p*p;
      typ qq       = 1. - 2.*q*q;
      typ dpq      = 2.*p*q;

      /******** In the orbital plane (see e.g. Laskar & Robutel 1995) ********/
      r       = a*(1. - e*e)/(1. + e*cosnu);
      g       = sqrt(mu*a*(1. - e*e));
      dnudt   = g/(r*r);
      drdt    = a*e*dnudt*sinnu*(1. - e*e)/((1. + e*cosnu)*(1. + e*cosnu));
      X_buff  = r*cosnu;
      Y_buff  = r*sinnu;
      vX_buff = drdt*cosnu - r*dnudt*sinnu;
      vY_buff = drdt*sinnu + r*dnudt*cosnu;
      
      /******** Rotations to convert to reference plane (see e.g. Laskar & Robutel 1995) ********/      
      X  =  X_buff*(pp*cosvarpi + dpq*sinvarpi) +  Y_buff*(dpq*cosvarpi - pp* sinvarpi);
      vX = vX_buff*(pp*cosvarpi + dpq*sinvarpi) + vY_buff*(dpq*cosvarpi - pp* sinvarpi);
      Y  =  X_buff*(qq*sinvarpi + dpq*cosvarpi) +  Y_buff*(qq*cosvarpi  - dpq*sinvarpi);
      vY = vX_buff*(qq*sinvarpi + dpq*cosvarpi) + vY_buff*(qq*cosvarpi  - dpq*sinvarpi);
      
      /******** Writing the cartesian coordinates ********/
      *(cart + 1) = X;
      *(cart + 2) = Y;
      *(cart + 3) = vX;
      *(cart + 4) = vY;
}


void cart2ell(typ * cart, typ * alkhqp, typ mu){

      /******** Computes the elliptic elements a, l, k, h, q & p where a is the semi-major axis and l is the mean      ********/
      /******** longitude. k, h, q & p are defined as k + ih = e*exp(i*varpi) and q + ip = sin(I/2)*exp(i*Omega). The  ********/
      /******** eccentricity e, inclination I, longitude of the pericentre varpi and longitude of the ascending node   ********/
      /******** Omega are straightforward to obtain from k, h, q & p, but these variables are better since they are    ********/
      /******** regular even et zero eccentricity and inclination. We have e = sqrt(k^2+h^2), sin(I/2) = sqrt(q^2+p^2),********/
      /******** varpi = atan2(h,k) and Omega = atan2(p,q). The cartesian coordinates are in the array cart.            ********/
      /******** The elliptic elements a, l, k, h, q & p are written in the vector alkhqp. This vector is overwritten   ********/
      /******** This function originally comes from the IMCCE lab (Paris Observatory) and was adapted for Aptidal      ********/

      typ X, Y, Z, vX, vY, vZ, R, V2, RV, C1, C2, C3, DC, CC, AA, aux0;
      typ a11, a12, a21, a22, c11, c12, c21, c22, K1, H1, K2, H2, K, H;
      typ SMU, FAC1, FAC2, b12, b22, sinF, cosF, F, USQA, aux1;

      /******** Getting the cartesian coordinates ********/
      X  = cart[1];
      Y  = cart[2];
      Z  = 0.;
      vX = cart[3];
      vY = cart[4];
      vZ = 0.;

      /******** Computing the semi-major axis ********/
      R   = sqrt(X*X + Y*Y + Z*Z);
      V2  = vX*vX + vY*vY + vZ*vZ;
      AA  = R*mu/(2.0*mu - R*V2); //Division by zero if the trajectory is perfectly parabolic.
      *alkhqp = AA;

      /******** Normalizing the velocities (Adopting the convention of J. Laskar's 2004 lectures notes) ********/
      SMU = sqrt(mu);
      vX /= SMU;
      vY /= SMU;
      vZ /= SMU;

      /******** Computing the angular momentum ********/
      V2 = vX*vX + vY*vY + vZ*vZ;
      RV = X*vX  + Y*vY  + Z*vZ;
      C1 = Y*vZ  - Z*vY;
      C2 = Z*vX  - X*vZ;
      C3 = X*vY  - Y*vX;
      CC = C1*C1 + C2*C2 + C3*C3;
      DC = sqrt(CC);

      /******** Computing (q, p) ********/
      aux0          = sqrt(2.0*(CC + DC*C3));
      *(alkhqp + 4) = (aux0 == 0. ? 0. : -C2/aux0);
      *(alkhqp + 5) = (aux0 == 0. ? 0. :  C1/aux0);

      if (R == 0. || V2 == 0.){
            K = 0.;
            H = 0.;
      }
      else{
            /******** Computing the matrix coefficients needed for (k, h) ********/
            a11 = V2 - 1.0/R;
            a12 = RV/(R*DC);
            a21 = -RV;
            a22 = DC - R/DC;

            /******** Computing (k, h) ********/
            c11  =  X*a11  + vX*a21;
            c12  =  X*a12  + vX*a22;
            c21  =  Y*a11  + vY*a21;
            c22  =  Y*a12  + vY*a22;
            FAC1 =  C1/(DC + C3);
            FAC2 =  C2/(DC + C3);
            K1   =  c11 - FAC1*(Z*a11 + vZ*a21);
            H1   = -c12 + FAC1*(Z*a12 + vZ*a22);
            H2   =  c21 - FAC2*(Z*a11 + vZ*a21); //Should be equal to H1
            K2   =  c22 - FAC2*(Z*a12 + vZ*a22); //Should be equal to K1
            if (fabs(H1 - H2) + fabs(K1 - K2) > 1.e-6){
                  printf("Warning : Bad computation of (k,h) in function cart2ell. (K1 - K2, H1 - H2) = (%.13lf, %.13lf)\n", K1 - K2, H1 - H2);
            }
            K    =  0.5*(K1 + K2);
            H    =  0.5*(H1 + H2);
      }
      *(alkhqp + 2) = K;
      *(alkhqp + 3) = H;

      if (R == 0. || V2 == 0.){
            *(alkhqp + 1) = 0.;
      }
      else{
            /******** Computing the mean longitude l = M + varpi ********/
            USQA = sqrt(2.0/R - V2);
            if ((USQA) >= 0.0){ //Elliptic case
                  b12  = vX - FAC1*vZ;
                  b22  = vY - FAC2*vZ;
                  aux1 = (R*V2 - 1.0)/(1.0 + DC*USQA);
                  sinF = -b12*R*USQA + H*aux1;
                  cosF =  b22*R*USQA + K*aux1;
                  F    =  atan2(sinF, cosF);
                  *(alkhqp + 1) = F - RV*USQA;
            }
            else{ //Hyperbolic case
                  USQA  = sqrt(-(2.0/R - V2));
                  typ E = atanh(RV*USQA/(R*V2 - 1.0));
                  typ M = sqrt(K*K + H*H) * sinh(E) - E;
                  *(alkhqp + 1) = M + atan2(H, K);
            }
      }
}


typ mean2true(typ M, typ mu, typ a, typ e){

      /******** Computes the true anomaly from the mean anomaly using the differential equation ********/
      /******** dnu/dt = sqrt(mu)*(1 + e cos nu)^2/(a(1-e^2))^(3/2) and a Runge-Kutta 4 method  ********/
      /******** This is equivalent to solving Kepler's equation.                                ********/

      M = fmod(M, 2.*M_PI);
      if (M <= -M_PI){ //Reducing to the range ]-pi, pi]
            M += 2.*M_PI;
      }
      else if (M > M_PI){
            M -= 2.*M_PI;
      }
      
      int j, N_step;
      typ K1, K2, K3, K4;
      typ period, t, dt, n_step, partial_tra, num;
      typ denom = sqrt(a*(1.0 - e*e));  denom *= denom*denom;
      typ sq_mu = sqrt(mu);
      typ n     = sq_mu/sqrt(a*a*a);
      typ time  = M/n;
      typ previous_tra = 0.;

      /******** Establishing the integration time and the timestep to be used ********/
      period       = 2.*M_PI/sq_mu*sqrt(a*a*a);
      n_step       = e < 1. ? floor(256.*fabs(time)/period) + 1. : 256.;
      dt           = time/n_step;
      n_step      *= e < 1. && e > 0.8 ? 4. : 1.; //Decreasing the timestep in case of highly eccentric elliptic trajectory
      dt          /= e < 1. && e > 0.8 ? 4. : 1.;
      N_step       = (int) n_step;
      if (fabs(time - n_step*dt) > 1.e-12){
            fprintf(stderr, "\nError: Wrong computation of the integration time in function mean2true.\n");
            abort();
      }
      
      /******** Integrating ********/
      for (j = 0; j < N_step; j ++){
            partial_tra   = previous_tra;
            num           = 1. + e*cos(partial_tra);
            K1            = sq_mu*num*num/denom;
            partial_tra   = previous_tra + 0.5*K1*dt;
            num           = 1. + e*cos(partial_tra);
            K2            = sq_mu*num*num/denom;
            partial_tra   = previous_tra + 0.5*K2*dt;
            num           = 1. + e*cos(partial_tra);
            K3            = sq_mu*num*num/denom;
            partial_tra   = previous_tra + K3*dt;
            num           = 1. + e*cos(partial_tra);
            K4            = sq_mu*num*num/denom;
            previous_tra += dt*(K1 + 2.*K2 + 2.*K3 + K4)/6.;
      }
      
      return previous_tra;
}


typ mean2eccentric(typ l, typ k, typ h){
      /******** Solves the Kepler equation without a numerical integration. Much faster than the above function ********/
      /******** Returns F = E + varpi                                                                           ********/
      /******** (l, k, h) = (Mean longitude, e*cos(varpi), e*sin(varpi))                                        ********/
      /******** Written by ASD Team LTE lab (former IMCCE)                                                      ********/

      int i    = 0;
      int imax = 20;
      typ eps, a, ca, sa, se, ce, fa, f1a, f2a, f3a, f4a, f5a, d1, d2, d3, d4, d5;

      eps  = 2.*2.26e-16;
      //Order 3 method
      a   = l;
      ca  = cos(a);
      sa  = sin(a);
      se  = k*sa - h*ca;
      ce  = k*ca + h*sa;
      fa  = a - se - l;
      f1a = 1. - ce;
      f2a = se/2.;
      f3a = ce/6.;
      d1  = -fa/f1a;
      d2  = -fa/(f1a - d1*f2a);
      d3  = -fa/(f1a + d2*(f2a + d2*f3a));
      a   = a + d3;
      //Order 6 method
      ca  = cos(a);
      sa  = sin(a);
      se  = k*sa - h*ca;
      ce  = k*ca + h*sa;
      fa  = a - se - l;
      f1a = 1. - ce;
      f2a = se/2.;
      f3a = ce/6.;
      f4a = -se/24.;
      f5a = -ce/120.;
      d1  = -fa/f1a;
      d2  = -fa/(f1a - d1*f2a);
      d3  = -fa/(f1a + d2*(f2a + d2*f3a));
      d4  = -fa/(f1a + d3*(f2a + d3*(f3a + d3*f4a)));
      d5  = -fa/(f1a + d4*(f2a + d4*(f3a + d4*(f4a + d4*f5a))));
      a   = a + d5;
      //Checking current precision
      while (1){
            i ++;
            ca  = cos(a);
            sa  = sin(a);
            se  = k*sa - h*ca;
            fa  = a - se - l;
            ce  = k*ca + h*sa;
            f1a = 1. - ce;
            d1  = -fa/f1a;
            //If precision is too low, computations are continued by iterating order 1 method
            if (fabs(d1)/max(1., fabs(a)) > eps){
                  if (i > imax){
                        fprintf(stderr, "\nWarning: Cannot reach a satisfactory precision in function mean2eccentric.\n");
                        return a;
                  }
                  a = a + d1;
            }
            else{
                  return a;
            }
      }
}


void prods(typ a, typ c, const typ * X, typ b, typ d, const typ * Y, typ * R, typ * S){

      /******** Auxiliary function to kepsaut. Implemented by M. Gastineau & F. Joutel (IMCCE) ********/
      
      int i, j;
      typ a1, a2, c1, c2, b1, b2, d1, d2;
      typ X1[3], X2[3], Y1[3], Y2[3], R1[3], R2[3], RR[3];
      typ u, v, temp, DR[3];

      typ base = (1L<<((DBL_MANT_DIG + 1)/2)) + 1L;
    
      a1 = (a - a*base) + a*base;
      a2 = a - a1;
      b1 = (b - b*base) + b*base;
      b2 = b - b1;
      c1 = (c - c*base) + c*base;
      c2 = c - c1;
      d1 = (d - d*base) + d*base;
      d2 = d - d1;
    
      for (j = 0; j < 3; j ++){
            X1[j]  = (X[j] - X[j]*base) + X[j]*base;
            X2[j]  =  X[j] - X1[j];
            Y1[j]  = (Y[j] - Y[j]*base) + Y[j]*base;
            Y2[j]  =  Y[j] - Y1[j];
            R1[j]  = a*X[j];
            RR[j]  = (a1*X1[j] - R1[j]) + (a1*X2[j]) + (a2*X1[j]) + (a2*X2[j]);
            R2[j]  = b*Y[j];
            RR[j] += b1*Y1[j] - R2[j] + b1*Y2[j] + b2*Y[j] + b2*Y2[j];
      }
      for (i = 0; i < 3; i ++){
            u = R1[i];
            v = R2[i];
            if (fabs(u) < fabs(v)){
                  temp = u;
                  u    = v;
                  v    = temp;
            }
            R [i] = u + v;
            DR[i] = v + (u - R[i]);
      }
      for (j = 0; j < 3; j ++){
            DR[j] += RR[j];
      }
      for (i = 0; i < 3; i ++){
            u = X[i];
            v = R[i];
            if (fabs(u) < fabs(v)){
                  temp = u;
                  u    = v;
                  v    = temp;
            }
            R [i] = u + v;
            RR[i] = v + (u - (R[i]));
      } 
      for (i = 0; i < 3; i ++){
            R[i] += DR[i] + RR[i];
      }
      for (i = 0; i < 3; i ++){
            R1[i]  = c*X[i];
            RR[i]  = c1*X1[i] - R1[i] + c1*X2[i] + c2*X1[i] + c2*X2[i];
            R2[i]  = d*Y[i];
            RR[i] += d1*Y1[i] - R2[i] + d1*Y2[i] + d2*Y1[i] + d2*Y2[i];
      }
      for (i = 0; i < 3; i ++){
            u = R1[i];
            v = R2[i];
            if (fabs(u) < fabs(v)){
                  temp = u;
                  u    = v;
                  v    = temp;
            }
            S [i] = u + v;
            DR[i] = v + (u - (S[i]));
      }  
      for (j = 0; j < 3; j ++){
            DR[j] += RR[j];
      }
      for (i = 0; i < 3; i ++){
            u = Y[i];
            v = S[i];
            if (fabs(u) < fabs(v)){
                  temp = u;
                  u    = v;
                  v    = temp;
            }
            S [i] = u + v;
            RR[i] = v + (u - (S[i]));
      }
      for (i = 0; i < 3; i ++){
            S[i] += DR[i] + RR[i];  
      }
}


void newt(typ DM, typ A, typ B, typ * const p_X, typ * const p_C, typ * const p_S, typ * const EXP, typ * const CM1, typ * const SMX, int bexitonerror){

      /******** Auxiliary function to kepsaut. Implemented by M. Gastineau, J. Couetdic & F. Joutel (IMCCE) ********/
      /******** Solves the equation F(X) = X - A*sin(X) - B*cos(X) + B - DM = 0                             ********/

      typ B0, B1, F, F1, X, C, S, D1, DIFF, ATEMP;
      int i, NMAX, niter;
      typ EPS = DBL_EPSILON;
      NMAX = 25;

      F  = DM;
      F1 = 1. - A;
      D1 = F/F1;
      X  = D1;
      
      i = 1;
      do{
            C    = cos(X);
            S    = sin(X);
            B0   = A*S + B*C;
            B1   = A*C - B*S;
            F    = B - DM + X - B0;
            F1   = 1. - B1;
            DIFF = F/F1;
            X   -= DIFF;
            i ++;
      } while (i <= NMAX && fabs(DIFF)/max(1., fabs(X)) >= 1.5*EPS);
      niter = i-1;
      if (i>NMAX){
            typ XQ  = X;
            typ AQ  = A;
            typ BQ  = B;
            typ DMQ = DM;
            typ DIFFQ;
            typ EPSQ = 2.5*EPS; //Modified from EPSQ = 1.5*EPS
            i = 1;
            do{
                  typ CQ  = cos(XQ);
                  typ SQ  = sin(XQ);
                  typ B0Q = AQ*SQ + BQ*CQ;
                  typ B1Q = AQ*CQ - BQ*SQ;
                  typ FQ  = BQ - DMQ + XQ - B0Q;
                  typ F1Q = 1. - B1Q;
                  DIFFQ   = FQ/F1Q;
                  XQ     -= DIFFQ;
                  i ++;
            } while (i <= NMAX && fabs(DIFFQ)/max(1., fabs(XQ)) >= EPSQ);
            niter += i;
        
            if (i > NMAX){
                  #ifndef __CUDACC__
                  if (bexitonerror){
                        fprintf(stderr, "\nError: Fatal error in function newt: X = %lf, DIFF = %lf, i = %d\n", X, DIFF, i);
                        abort();
                  }
                  else{
                        printf("\nError: Fatal error in function newt: X = %lf, DIFF = %lf, i = %d\n", X, DIFF, i);
                  }
                  #endif
                  niter = -1;
            }
            else{
                  X = XQ;
            }
      }
        
      C     = cos(X);
      S     = sin(X);
      *EXP  = F1 - DIFF*B0;
      ATEMP = sin(X*0.5);
      *CM1  = -2.*ATEMP*ATEMP;    
      #if REALTYPE==0
      if (fabs(X) < 0.2){               
         *SMX = ((((2623./4929724800.*X*X
                 -97./1053360.)*X*X
                 +59./9880.)*X*X
                 -1./6.)*X*X*X)
                /
                ((23./326040.*X*X
                  +7./494)*X*X
                 +1.
                );
      }   
      else 
      #endif
      {
            *SMX = S - X;
      }
      *p_C = C;
      *p_S = S;
      *p_X = X;
}


void kepsaut(typ * cart, typ mu, typ dt){

      /******** Updates the cartesian coordinates cart = (x, y, z, vx, vy, vz) along ********/
      /******** a Keplerian arc for a time dt. mu is the gravitational parameter.    ********/
      /******** This function was provided by Mickaël Gastineau (IMCCE, ASD team)    ********/

      typ X [3] = {cart[1], cart[2], 0.};
      typ XD[3] = {cart[3], cart[4], 0.};
      typ smu;
      typ R, V2, XV, a, ce, se, dm, c, s, e;
      typ a1, a2, a3, a4, Xn[3], Xdn[3], delta, sqa;    
      typ cm1, smx;
      int j;
      
      smu = sqrt(mu);
      R   = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
      V2  = (XD[0]*XD[0] + XD[1]*XD[1] + XD[2]*XD[2])/mu;
      XV  = (X [0]*XD[0] + X [1]*XD[1] + X [2]*XD[2])/smu;
      a   = R/(2. - R*V2);
      sqa = sqrt(a);
      ce  = R*V2 - 1.;
      se  = XV/sqa;
      dm  = dt*smu/(sqa*a);
      newt(dm, ce, se, &delta, &c, &s, &e, &cm1, &smx, 0);

      a1  = cm1/(2. - R*V2); 
      a2  = dt + smx*sqa*a/smu; 
      a3  = -s/((2. - R*V2)*e*a*sqa)*smu;
      a4  = cm1/e;
      prods(a1, a3, X, a2, a4, XD, Xn, Xdn);

      for(j = 0; j < 2; j ++){
            cart[j + 1] = Xn [j];
            cart[j + 3] = Xdn[j];
      }
}


void exp_tau_LB(typ tau, typ * X_cart){

      /******** Performs one step of a SABA1 on the perturbation B = B_1(p) + B_2(q) ********/

      int k, j;
      typ sum_rtilde_x = 0.;  typ sum_rtilde_y = 0.;
      typ norm, norm3, Xk, Xj, Yk, Yj, dX, dY;

      /******** Computing the sum of the \tilde{r}_k/m_0 ********/
      for (k = 1; k <= how_many_planet; k ++){
            sum_rtilde_x += masses[k]*X_cart[4*k - 1]/m0;
            sum_rtilde_y += masses[k]*X_cart[4*k]    /m0;
      }

      /******** Step exp(1/2*tau*L_B1) ********/
      for (k = 1; k <= how_many_planet; k ++){
            X_cart[4*k - 3] += tau/2.*(sum_rtilde_x - masses[k]*X_cart[4*k - 1]/m0);
            X_cart[4*k - 2] += tau/2.*(sum_rtilde_y - masses[k]*X_cart[4*k]    /m0);
      }

      /******** Step exp(tau*L_B2) ********/
      for (k = 1; k <= how_many_planet; k ++){
            Xk = X_cart[4*k - 3];  Yk = X_cart[4*k - 2];
            for (j = 1; j <= how_many_planet; j ++){
                  if (j != k){
                        Xj   = X_cart[4*j - 3];  Yj = X_cart[4*j - 2];
                        dX   = Xk - Xj;  dY = Yk - Yj;
                        norm = sqrt(dX*dX + dY*dY);  norm3 = norm*norm*norm;
                        X_cart[4*k - 1] -= tau*G*masses[j]*dX/norm3;
                        X_cart[4*k]     -= tau*G*masses[j]*dY/norm3;
                  }
            }
      }

      /******** Step exp(1/2*tau*L_B1) ********/
      for (k = 1; k <= how_many_planet; k ++){
            X_cart[4*k - 3] += tau/2.*(sum_rtilde_x - masses[k]*X_cart[4*k - 1]/m0);
            X_cart[4*k - 2] += tau/2.*(sum_rtilde_y - masses[k]*X_cart[4*k]    /m0);
      }
}


void UnaveragedSABAn(typ tau, typ T, int output_step, typ * X_old, int n){

      /******** Integrates the complete Hamiltonian with a SABAn where      ********/
      /******** 1 <= n <= 6 for a time T with a timestep tau. Outputs every ********/
      /******** output_step timestep to file pth/UnaveragedSABAn_chain.txt. ********/
      /******** The Keplerian part A is integrated exactly, but the         ********/
      /******** perturbative part takes the form B = B_1(p) + B_2(q) and is ********/
      /******** integrated approximetaly but symplectically with a SABA1.   ********/


      char file_path[800];
      char nn[10];
      FILE * file;
      int N_step, iter, i;
      typ X_cart[4*how_many_planet + 1];
      typ X_init[4*how_many_planet + 1];
      typ X_buff[4*how_many_planet + 1];
      typ X_new [4*how_many_planet + 1];
      typ X_uv  [4*how_many_planet + 1];
      typ lbdOld[  how_many_planet + 1];
      typ   gOld[  how_many_planet + 1];
      typ a, e, vp, M, mu, nu, beta, lbd, g, dist2init;
      typ c1, c2, c3, c4, d1, d2, d3;
      typ alkhqp[6];
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "UnaveragedSABA");
      sprintf(nn, "%d", n);
      strcat(file_path, nn);
      strcat(file_path, "_");
      for (i = 1; i <= how_many_planet; i ++){
            sprintf(nn, "%d", p_i[i]);
            strcat(file_path, nn);
      }
      
      /******** To be removed ********/
      /*time_t t;
      time(&t);
      char tt[20];
      sprintf(tt, "%d", ((int) t) - 1744291130);
      strcat(file_path, "_");
      strcat(file_path, tt);*/
      
      strcat(file_path, ".txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "\nError: Cannot create or open file UnaveragedSABAn_chain.txt in function UnaveragedSABAn.\n");
            abort();
      }
      
      /******** Initializing constants (Laskar & Robutel, 2001) ********/
      c1 = 0.;  c2 = 0.;  c3 = 0.;  c4 = 0.;  d1 = 0.;  d2 = 0.;  d3 = 0.;
      if      (n == 1){c1 = 0.5;                   d1 = 1.;}
      else if (n == 2){c1 = 0.5 - sqrt(3.)/6.;     d1 = 0.5;                   c2 = sqrt(3.)/3.;}
      else if (n == 3){c1 = 0.5 - sqrt(15.)/10.;   d1 = 5./18.;                c2 = sqrt(15.)/10.;         d2 = 4./9.;}
      else if (n == 4){c1 = 0.0694318442029737124; d1 = 0.1739274225687269287; c2 = 0.2605776340045981552; d2 = 0.3260725774312730713; c3 = 0.3399810435848562648;}
      else if (n == 5){c1 = 0.0469100770306680036; d1 = 0.1184634425280945438; c2 = 0.1838552679164904509; d2 = 0.2393143352496832340; c3 = 0.2692346550528415455; d3 = 64./225.;}
      else if (n == 6){c1 = 0.0337652428984239861; d1 = 0.0856622461895851725; c2 = 0.1356300638684437571; d2 = 0.1803807865240693038; c3 = 0.2112951001915338025;
                       d3 = 0.2339569672863455237; c4 = 0.2386191860831969086;}
      else{fprintf(stderr, "\nError: n must be between 1 and 6 in function UnaveragedSABAn.\n");  abort();}

      /******** Checking validity of coefficients. To be removed when robust ********/
      if (fabs(2.*c1 + (n == 2 ? 1. : 2.)*c2 + (n == 4 ? 1. : 2.)*c3 + c4 - 1.) > 1.e-15 || fabs((n == 1 ? 1. : 2.)*d1 + (n == 3 ? 1. : 2.)*d2 + (n == 5 ? 1. : 2.)*d3 - 1.) > 1.e-15){
            fprintf(stderr, "\nError: The sum of the coefficients does not seem to be 1 in function UnaveragedSABAn.\n");
            abort();
      }
      
      /******** Initializing the cartesian coordinates ********/
      for (i = 1; i <= how_many_planet; i ++){
            mu = G*(m0 + masses[i]);
            a  = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            e  = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
            vp = -X_old[4*i - 2];
            M  =  X_old[4*i - 3] - vp;
            nu = mean2true(M, mu, a, e);
            //E  = mean2eccentric(l, e*cos(vp), e*sin(vp));
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
            lbdOld[i] = X_old[4*i - 3];  gOld[i] = X_old[4*i - 2];
      }
      for (i = 1; i <= 4*how_many_planet; i ++){ //X_init remembers the initial position in order to compute the distance with respect to it in the phase space
            X_init[i] = X_cart[i];
      }

      /******** Integrating ********/
      N_step = (int) ceil(T/tau);
      for (iter = 0; iter < N_step; iter ++){
      
            /******** Writing to file ********/
            if (iter%output_step == 0){
                  for (i = 1; i <= how_many_planet; i ++){
                        X_buff[4*i - 3] = X_cart[4*i - 3]; X_buff[4*i - 2] = X_cart[4*i - 2]; X_buff[4*i - 1] = X_cart[4*i - 1]; X_buff[4*i] = X_cart[4*i];
                  }
                  if (iter){ //Need to perform a step exp(-c1*tau*L_A) on the buffer before outputting
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_buff + 4*i - 4, mu, -c1*tau);
                        }
                  }
                  else{
                        fprintf(file, "Numerical integration of the unaveraged Hamiltonian with a SABA%d integrator in Poincaré's cartesian canonical coordinates.\n", n);
                        fprintf(file, "The Keplerian part A is integrated exactly using function kepsaut. The perturbative part B is not integrable but can be\n");
                        fprintf(file, "written B = B1 + B2 with B1 and B2 both integrable. Therefore, B is integrated approximately but symplectically with a SABA1.\n");
                        fprintf(file, "\n");
                        fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 2 + 4*how_many_planet, how_many_planet);
                        fprintf(file, "Time, distance to initial point in the phase space, a_j, e_j, phi_j, sig_j\n");
                        fprintf(file, "\n");
                  }
                  for (i = 1; i <= how_many_planet; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
                        beta = m0*masses[i]/(m0 + masses[i]);
                        mu   = G*(m0 + masses[i]);
                        cart2ell(X_buff + 4*i - 4, alkhqp, mu);
                        e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
                        lbd  = continuousAngle(alkhqp[1],                  lbdOld[i]);
                        g    = continuousAngle(-atan2(alkhqp[3], alkhqp[2]), gOld[i]);
                        X_old[4*i - 3] = lbd;
                        X_old[4*i - 2] = g;
                        X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
                        X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
                        lbdOld[i]      = lbd;  gOld[i] = g;
                  }
                  dist2init = 0.;
                  for (i = 1; i <= 4*how_many_planet; i ++){
                        dist2init += (X_buff[i] - X_init[i])*(X_buff[i] - X_init[i]);
                  }
                  dist2init = sqrt(dist2init);
                  old2new(X_old, X_new, X_uv);             // (lbd_j, -vrp_j; Lbd_j, D_j) -> (phi_j, v_j; Phi_j, u_j)
                  fprintf(file, "%.9lf %.14lf", tau*(typ) iter, dist2init);
                  for (i = 1; i <= how_many_planet; i ++){
                        e   = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
                        a   = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                        fprintf(file, " %.14lf %.14lf %.14lf %.14lf", a, e, X_new[4*i - 3], X_new[4*i - 2]);
                  }
                  fprintf(file, "\n");
            }
            
            /******** Step exp(c1*tau*L_A). For the first iteration only ********/
            if (!iter){
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c1*tau);
                  }
            }
            
            /******** Step exp(d1*tau*L_B) ********/
            exp_tau_LB(d1*tau, X_cart);
            
            if (n >= 2){
                  /******** Step exp(c2*tau*L_A) ********/
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c2*tau);
                  }

                  if (n >= 3){
                        /******** Step exp(d2*tau*L_B) ********/
                        exp_tau_LB(d2*tau, X_cart);

                        if (n >= 4){
                              /******** Step exp(c3*tau*L_A) ********/
                              for (i = 1; i <= how_many_planet; i ++){
                                    mu = G*(m0 + masses[i]);
                                    kepsaut(X_cart + 4*i - 4, mu, c3*tau);
                              }

                              if (n >= 5){
                                    /******** Step exp(d3*tau*L_B) ********/
                                    exp_tau_LB(d3*tau, X_cart);
                                    
                                    if (n >= 6){
                                          /******** Step exp(c4*tau*L_A) ********/
                                          for (i = 1; i <= how_many_planet; i ++){
                                                mu = G*(m0 + masses[i]);
                                                kepsaut(X_cart + 4*i - 4, mu, c4*tau);
                                          }
                                          
                                          /******** Step exp(d3*tau*L_B) ********/
                                          exp_tau_LB(d3*tau, X_cart);
                                    }
                                    
                                    /******** Step exp(c3*tau*L_A) ********/
                                    for (i = 1; i <= how_many_planet; i ++){
                                          mu = G*(m0 + masses[i]);
                                          kepsaut(X_cart + 4*i - 4, mu, c3*tau);
                                    }
                              }
                              /******** Step exp(d2*tau*L_B) ********/
                              exp_tau_LB(d2*tau, X_cart);
                        }
                        
                        /******** Step exp(c2*tau*L_A) ********/
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_cart + 4*i - 4, mu, c2*tau);
                        }
                  }
                  
                  /******** Step exp(d1*tau*L_B) ********/
                  exp_tau_LB(d1*tau, X_cart);
            }
            
            /******** Step exp(2*c1*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, 2.*c1*tau);
            }
      }

      /******** Closing output file ********/
      fclose(file);
}


void get_frequencies(typ tau, typ T, typ * X_old, int n){

      /******** Same as above but calculates the fast frequencies with a NAFF    ********/
      /******** method (Laskar, Froeschlé, Celletti, 1992). If nu_fast is zero,  ********/
      /******** then a first approximation of it is found and a Newtow-Raphson   ********/
      /******** method is applied on d< f(t), e^i*nu*t >/dnu in order to obtain  ********/
      /******** the exact value. Otherwise, a Newton-Raphson is applied directly ********/
      /******** I consider f(t) = sqrt(u_1**2 + v_1**1 + ... + u_N**2 + v_N**2)  ********/
      
      /* To be generalized to more complicated chains later */
      
      int N_step, iter, i, p;
      typ X_cart  [4*how_many_planet + 1];
      typ X_buff  [4*how_many_planet + 1];
      typ X_new   [4*how_many_planet + 1];
      typ X_uv    [4*how_many_planet + 1];
      typ dH_old  [4*how_many_planet + 1];
      typ dH_rect [4*how_many_planet + 1];
      typ dH_polar[4*how_many_planet + 1];
      typ lbdOld[    how_many_planet + 1];
      typ   gOld[    how_many_planet + 1];
      typ a, e, vp, M, mu, nu, beta, lbd, g, delta;
      typ ak_fast, bk_fast, ak_reso, bk_reso, t, S, dSda, dSdb, par, learning_rate, oldS, N;
      typ c1, c2, c3, c4, d1, d2, d3;
      typ alkhqp[6];
      
      /******** Initializing constants ********/
      c1 = 0.;  c2 = 0.;  c3 = 0.;  c4 = 0.;  d1 = 0.;  d2 = 0.;  d3 = 0.;
      if      (n == 1){c1 = 0.5;                   d1 = 1.;}
      else if (n == 2){c1 = 0.5 - sqrt(3.)/6.;     d1 = 0.5;                   c2 = sqrt(3.)/3.;}
      else if (n == 3){c1 = 0.5 - sqrt(15.)/10.;   d1 = 5./18.;                c2 = sqrt(15.)/10.;         d2 = 4./9.;}
      else if (n == 4){c1 = 0.0694318442029737124; d1 = 0.1739274225687269287; c2 = 0.2605776340045981552; d2 = 0.3260725774312730713; c3 = 0.3399810435848562648;}
      else if (n == 5){c1 = 0.0469100770306680036; d1 = 0.1184634425280945438; c2 = 0.1838552679164904509; d2 = 0.2393143352496832340; c3 = 0.2692346550528415455; d3 = 64./225.;}
      else if (n == 6){c1 = 0.0337652428984239861; d1 = 0.0856622461895851725; c2 = 0.1356300638684437571; d2 = 0.1803807865240693038; c3 = 0.2112951001915338025;
                       d3 = 0.2339569672863455237; c4 = 0.2386191860831969086;}
      else{fprintf(stderr, "\nError: n must be between 1 and 6 in function get_frequencies.\n");  abort();}
      
      /******** Initializing the cartesian coordinates ********/
      for (i = 1; i <= how_many_planet; i ++){
            mu = G*(m0 + masses[i]);
            a  = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            e  = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
            vp = -X_old[4*i - 2];
            M  =  X_old[4*i - 3] - vp;
            nu = mean2true(M, mu, a, e);
            //E  = mean2eccentric(l, e*cos(vp), e*sin(vp));
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
            lbdOld[i] = X_old[4*i - 3];  gOld[i] = X_old[4*i - 2];
      }
      
      /******** Integrating ********/
      tau       /= 2.;
      N_step     = (int) ceil(T/tau);
      N_step    += N_step % 2 == 0 ? 0 : 1;
      tau        = T/ (typ) N_step;
      typ * phi2 = (typ *)malloc((2*N_step + 1)*sizeof(typ));
      typ * phi3 = (typ *)malloc((2*N_step + 1)*sizeof(typ));
      if (phi2 == NULL || phi3 == NULL){
            fprintf(stderr, "\nError: Could not allocate memory for arrays in function get_frequencies.\n");
            abort();
      }
      for (iter = -N_step + 1; iter <= N_step; iter ++){
            
            /******** For the first iteration only ********/
            if (iter == -N_step + 1){
                  /******** Initializing phi2(-T) and phi3(-T) ********/
                  old2new(X_old, X_new, X_uv);
                  phi2[0] = X_uv[4*subchain[how_many_resonant - 1] - 3];
                  phi3[0] = X_uv[4*subchain[how_many_resonant]     - 3];
                  /******** Step exp(c1*tau*L_A) ********/
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c1*tau);
                  }
            }
            
            /******** Step exp(d1*tau*L_B) ********/
            exp_tau_LB(d1*tau, X_cart);
            
            if (n >= 2){
                  /******** Step exp(c2*tau*L_A) ********/
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c2*tau);
                  }

                  if (n >= 3){
                        /******** Step exp(d2*tau*L_B) ********/
                        exp_tau_LB(d2*tau, X_cart);

                        if (n >= 4){
                              /******** Step exp(c3*tau*L_A) ********/
                              for (i = 1; i <= how_many_planet; i ++){
                                    mu = G*(m0 + masses[i]);
                                    kepsaut(X_cart + 4*i - 4, mu, c3*tau);
                              }

                              if (n >= 5){
                                    /******** Step exp(d3*tau*L_B) ********/
                                    exp_tau_LB(d3*tau, X_cart);
                                    
                                    if (n >= 6){
                                          /******** Step exp(c4*tau*L_A) ********/
                                          for (i = 1; i <= how_many_planet; i ++){
                                                mu = G*(m0 + masses[i]);
                                                kepsaut(X_cart + 4*i - 4, mu, c4*tau);
                                          }
                                          
                                          /******** Step exp(d3*tau*L_B) ********/
                                          exp_tau_LB(d3*tau, X_cart);
                                    }
                                    
                                    /******** Step exp(c3*tau*L_A) ********/
                                    for (i = 1; i <= how_many_planet; i ++){
                                          mu = G*(m0 + masses[i]);
                                          kepsaut(X_cart + 4*i - 4, mu, c3*tau);
                                    }
                              }
                              /******** Step exp(d2*tau*L_B) ********/
                              exp_tau_LB(d2*tau, X_cart);
                        }
                        
                        /******** Step exp(c2*tau*L_A) ********/
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_cart + 4*i - 4, mu, c2*tau);
                        }
                  }
                  
                  /******** Step exp(d1*tau*L_B) ********/
                  exp_tau_LB(d1*tau, X_cart);
            }
            
            /******** Step exp(2*c1*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, 2.*c1*tau);
            }
            
            /******** Step exp(-c1*tau*L_A) on buffer ********/
            for (i = 1; i <= how_many_planet; i ++){
                  X_buff[4*i - 3] = X_cart[4*i - 3];  X_buff[4*i - 2] = X_cart[4*i - 2];  X_buff[4*i - 1] = X_cart[4*i - 1];  X_buff[4*i] = X_cart[4*i];
            }
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_buff + 4*i - 4, mu, -c1*tau);
            }
            
            /******** Getting relevant coordinates from cartesian coordinates ********/
            for (i = 1; i <= how_many_planet; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
                  beta = m0*masses[i]/(m0 + masses[i]);
                  mu   = G*(m0 + masses[i]);
                  cart2ell(X_buff + 4*i - 4, alkhqp, mu);
                  e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
                  lbd  = continuousAngle(alkhqp[1],                  lbdOld[i]);
                  g    = continuousAngle(-atan2(alkhqp[3], alkhqp[2]), gOld[i]);
                  X_old[4*i - 3] = lbd;
                  X_old[4*i - 2] = g;
                  X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
                  X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
                  lbdOld[i]      = lbd;  gOld[i] = g;
            }
            old2new(X_old, X_new, X_uv);
            
            /******** Computing phi2(t) and phi3(t) ********/
            phi2[iter + N_step] = X_uv[4*subchain[how_many_resonant - 1] - 3];
            phi3[iter + N_step] = X_uv[4*subchain[how_many_resonant]     - 3];
      }

      
      /******** Determination of nu_fast and nu_reso via a linear fit phi = a*t + b ********/
      if (how_many_resonant >= 2){
            /******** I obtain a first determination of a by computing dHK/dPhi2 and of b by computing the value at the origin ********/
            dHdold(dH_old, X_old, 0);
            old2new(X_old, X_new, X_uv);
            dHdnew(dH_polar, dH_rect, dH_old, X_new, X_uv);
            ak_fast = dH_rect[4*subchain[how_many_resonant - 1] - 1];
            bk_fast = phi2[N_step];
            ak_reso = dH_rect[4*subchain[how_many_resonant]     - 1];
            bk_reso = phi3[N_step];
            //printf("ak_fast = %.13lf, bk_fast = %.13lf\n", ak_fast, bk_fast); //To be removed
            //printf("ak_reso = %.13lf, bk_reso = %.13lf\n", ak_reso, bk_reso); //To be removed
            /******** I now refine a and b for the fast frequency by minimizing S(a,b) = sum_{all points} (a*t + b - phi2(t))^2 via a gradient descent ********/
            delta  = 1.;  p = 0;  learning_rate = 0.0000000298023223876953125;  oldS = 1.e300;  N = (typ) (2*N_step + 1);
            while(delta > 1.e-11 && p < 512){
                  S = 0.;  dSda = 0.; dSdb = 0.;
                  for (iter = -N_step; iter <= N_step; iter ++){
                        t     = tau* (typ) iter;
                        par   = ak_fast*t + bk_fast - phi2[iter + N_step];
                        S    += par*par/N;
                        dSda += 2.*par*t/N;
                        dSdb += 2.*par/N;
                  }
                  //printf("p = %d, S = %.13lf\n", p, S); //To be removed
                  if (S > oldS){
                        learning_rate /= 2.;
                  }
                  delta    = fabs(S - oldS);
                  oldS     = S;
                  ak_fast -= learning_rate*dSda; //Gradient descent
                  bk_fast -= learning_rate*dSdb;
                  p ++;
                  if (p > 256 && delta > 1.e-8){
                        learning_rate /= 2.;
                  }
                  if (p == 512 && delta > 1.e-8){
                        fprintf(stderr, "\nError: The gradient descent does not seem to converge in function get_frequencies. The learning rate is probably ill-chosen.\n");
                        abort();
                  }
            }
            /******** I now refine a and b for the slow frequency by minimizing S(a,b) = sum_{all points} (a*t + b - phi3(t))^2 via a gradient descent ********/
            delta  = 1.;  p = 0;  learning_rate = 0.0000000298023223876953125;  oldS = 1.e300;  N = (typ) (2*N_step + 1);
            while(delta > 1.e-11 && p < 512){
                  S = 0.;  dSda = 0.; dSdb = 0.;
                  for (iter = -N_step; iter <= N_step; iter ++){
                        t     = tau* (typ) iter;
                        par   = ak_reso*t + bk_reso - phi3[iter + N_step];
                        S    += par*par/N;
                        dSda += 2.*par*t/N;
                        dSdb += 2.*par/N;
                  }
                  //printf("p = %d, S = %.13lf\n", p, S); //To be removed
                  if (S > oldS){
                        learning_rate /= 2.;
                  }
                  delta    = fabs(S - oldS);
                  oldS     = S;
                  ak_reso -= learning_rate*dSda; //Gradient descent
                  bk_reso -= learning_rate*dSdb;
                  p ++;
                  if (p > 256 && delta > 1.e-8){
                        learning_rate /= 2.;
                  }
                  if (p == 512 && delta > 1.e-8){
                        fprintf(stderr, "\nError: The gradient descent does not seem to converge in function get_frequencies. The learning rate is probably ill-chosen.\n");
                        abort();
                  }
            }
            nu_fast = ak_fast;
            nu_reso = ak_reso;
            //printf("nu_fast = %.13lf,  nu_reso = %.13lf\n", nu_fast, nu_reso); //To be removed
      }
      
      /******** Newton-Raphson on d<f(t), e^i*nu*t>/dnu ********/  //Probably useless. Will decide later
      /*delta = 1.; nu = nu_fast; p = 0;
      while (delta > 1.e-9){
            num = 0.;  denom = 0.;*/
            /******** Computing relevant integrals ********/
            /*for (iter = -N_step + 2; iter <= N_step; iter += 2){
                  t0       = tau* (typ) (iter - 2);
                  t1       = tau* (typ) (iter - 1);
                  t2       = tau* (typ) iter;
                  Hanning0 = Cp*fast_pow(1. + cos(M_PI*t0/T), Hanning_order);
                  Hanning1 = Cp*fast_pow(1. + cos(M_PI*t1/T), Hanning_order);
                  Hanning2 = Cp*fast_pow(1. + cos(M_PI*t2/T), Hanning_order);
                  num     += t0*   Hanning0*sin(nu*t0)*f[iter + N_step - 2] + 4.*t1*   Hanning1*sin(nu*t1)*f[iter + N_step - 1] + t2*   Hanning2*sin(nu*t2)*f[iter + N_step];
                  denom   += t0*t0*Hanning0*cos(nu*t0)*f[iter + N_step - 2] + 4.*t1*t1*Hanning1*cos(nu*t1)*f[iter + N_step - 1] + t2*t2*Hanning2*cos(nu*t2)*f[iter + N_step];
            }
            num  *= tau/(6.*T);  denom *= tau/(6.*T);
            nu   -= num/denom;
            printf("nu_fast  = %.13lf,  T  = %.13lf\n", nu, 2.*M_PI/nu); //To be removed
            delta = fabs(num/denom);
            p ++;
            if (p > 20 && delta > 1.e-9){
                  fprintf(stderr, "\nError: Newton-Raphson method does not seem to converge in function get_frequencies.\n");
                  abort();
            }
      }
      nu_fast = nu;*/
      
      free(phi2); phi2 = NULL;
      free(phi3); phi3 = NULL;
}


typ UnaveragedSABAn_NAFF(typ tau, typ T, int Hanning_order, typ * X_uv, typ * X_old, int n, int how_many_harmonics){

      /******** Same as above but computes and stores the average and the ********/
      /******** amplitude of the fast terms instead of writing to file.   ********/
      /******** Uses the NAFF method of (Laskar, Froeschlé, Celletti 1992)********/
      /******** A Hanning filter of order Hanning_order is used. The time ********/
      /******** step is tau/2 to use a Simpson 3pt method to compute      ********/
      /******** integrals. Fills X_uv with (phi_j, v_j; Phi_j, u_j) eva-  ********/
      /******** luated at t = 0 with fast terms and average only          ********/

      /* To be generalized to more complicated chains later */

      int N_step, iter, i, j;
      int N = how_many_planet;
      typ X_cart[4*N + 1];
      typ X_buff[4*N + 1];
      typ X_new [4*N + 1];
      typ X_uv_0[4*N + 1]; //Left   of the segment
      typ X_uv_1[4*N + 1]; //Middle of the segment
      typ X_uv_2[4*N + 1]; //Right  of the segment
      typ lbdOld[  N + 1];
      typ   gOld[  N + 1];
      typ As[(1 + how_many_harmonics)*(1 + 4*N)]; //Index j*(1 + 4*N) + 4*i - 3 (resp. -2, -1, -0) contains amplitude A_j of phi_i (resp. v_i, Phi_i, u_i)
      typ as[1 + 4*N]; //Index 4*i - 3 (resp. -2, -1, -0) contains amplitude A_j of phi_i (resp. v_i, Phi_i, u_i) for slow term. f(t)=A_0+\sum_j (A_j*cos(nu_fast*t)+B_j*sin(nu_fast*t))
      typ a, e, vp, M, mu, nu, beta, lbd, g, toBeReturnedNum, toBeReturnedNom;
      typ facto   [11] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3628800.};
      typ two_to_p[6]  = {1., 2., 4., 8., 16., 32.};
      typ Cp, H0, H1, H2, t0, t1, t2;
      typ c1, c2, c3, c4, d1, d2, d3;
      typ alkhqp[6];
      
      /******** Getting the fast frequencies ********/
      tau    /= 2.;
      N_step  = (int) ceil(T/tau);
      N_step += N_step % 2 == 0 ? 0 : 1;
      tau     = T/(typ) N_step;
      if (how_many_harmonics){
            for (i = 1; i <= how_many_planet; i ++){
                  X_buff[4*i - 3] = X_old[4*i - 3];  X_buff[4*i - 2] = X_old[4*i - 2];  X_buff[4*i - 1] = X_old[4*i - 1];  X_buff[4*i] = X_old[4*i];
            }
            get_frequencies(2.*tau, T, X_buff, n);
      }
      
      /******** Initializing arrays As and as *********/
      for (i = 0; i < (1 + how_many_harmonics)*(1 + 4*N); i ++){
            As[i] = 0.; //Initializing As
      }
      for (i = 0; i < 1 + 4*N; i ++){
            as[i] = 0.; //Initializing as
      }
      
      /******** Initializing constants ********/
      c1 = 0.;  c2 = 0.;  c3 = 0.;  c4 = 0.;  d1 = 0.;  d2 = 0.;  d3 = 0.;
      if      (n == 1){c1 = 0.5;                   d1 = 1.;}
      else if (n == 2){c1 = 0.5 - sqrt(3.)/6.;     d1 = 0.5;                   c2 = sqrt(3.)/3.;}
      else if (n == 3){c1 = 0.5 - sqrt(15.)/10.;   d1 = 5./18.;                c2 = sqrt(15.)/10.;         d2 = 4./9.;}
      else if (n == 4){c1 = 0.0694318442029737124; d1 = 0.1739274225687269287; c2 = 0.2605776340045981552; d2 = 0.3260725774312730713; c3 = 0.3399810435848562648;}
      else if (n == 5){c1 = 0.0469100770306680036; d1 = 0.1184634425280945438; c2 = 0.1838552679164904509; d2 = 0.2393143352496832340; c3 = 0.2692346550528415455; d3 = 64./225.;}
      else if (n == 6){c1 = 0.0337652428984239861; d1 = 0.0856622461895851725; c2 = 0.1356300638684437571; d2 = 0.1803807865240693038; c3 = 0.2112951001915338025;
                       d3 = 0.2339569672863455237; c4 = 0.2386191860831969086;}
      else{fprintf(stderr, "\nError: n must be between 1 and 6 in function UnaveragedSABAn_NAFF.\n");  abort();}
      
      /******** Initializing constant of Hanning filter ********/
      if (Hanning_order > 5){
            fprintf(stderr, "\nError: The order of the Hanning filter cannot exceed 5.\n");
            abort();
      }
      Cp = two_to_p[Hanning_order]*facto[Hanning_order]*facto[Hanning_order]/facto[2*Hanning_order];
      
      /******** Initializing the cartesian coordinates ********/
      for (i = 1; i <= how_many_planet; i ++){
            mu = G*(m0 + masses[i]);
            a  = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            e  = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
            vp = -X_old[4*i - 2];
            M  =  X_old[4*i - 3] - vp;
            nu = mean2true(M, mu, a, e);
            //E  = mean2eccentric(l, e*cos(vp), e*sin(vp));
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
            lbdOld[i] = X_old[4*i - 3];  gOld[i] = X_old[4*i - 2];
      }
      
      /******** Initializing output to 0 ********/
      for (i = 1; i <= 4*how_many_planet; i ++){
            X_uv[i] = 0.;
      }

      /******** Integrating ********/
      for (iter = -N_step + 1; iter <= N_step; iter ++){
            
            /******** For the first iteration only ********/
            if (iter == -N_step + 1){
                  /******** Initializing X_uv_0 ********/
                  old2new(X_old, X_new, X_uv_0);
                  /******** Step exp(c1*tau*L_A) ********/
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c1*tau);
                  }
            }
            
            /******** Step exp(d1*tau*L_B) ********/
            exp_tau_LB(d1*tau, X_cart);
            
            if (n >= 2){
                  /******** Step exp(c2*tau*L_A) ********/
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c2*tau);
                  }

                  if (n >= 3){
                        /******** Step exp(d2*tau*L_B) ********/
                        exp_tau_LB(d2*tau, X_cart);

                        if (n >= 4){
                              /******** Step exp(c3*tau*L_A) ********/
                              for (i = 1; i <= how_many_planet; i ++){
                                    mu = G*(m0 + masses[i]);
                                    kepsaut(X_cart + 4*i - 4, mu, c3*tau);
                              }

                              if (n >= 5){
                                    /******** Step exp(d3*tau*L_B) ********/
                                    exp_tau_LB(d3*tau, X_cart);
                                    
                                    if (n >= 6){
                                          /******** Step exp(c4*tau*L_A) ********/
                                          for (i = 1; i <= how_many_planet; i ++){
                                                mu = G*(m0 + masses[i]);
                                                kepsaut(X_cart + 4*i - 4, mu, c4*tau);
                                          }
                                          
                                          /******** Step exp(d3*tau*L_B) ********/
                                          exp_tau_LB(d3*tau, X_cart);
                                    }
                                    
                                    /******** Step exp(c3*tau*L_A) ********/
                                    for (i = 1; i <= how_many_planet; i ++){
                                          mu = G*(m0 + masses[i]);
                                          kepsaut(X_cart + 4*i - 4, mu, c3*tau);
                                    }
                              }
                              /******** Step exp(d2*tau*L_B) ********/
                              exp_tau_LB(d2*tau, X_cart);
                        }
                        
                        /******** Step exp(c2*tau*L_A) ********/
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_cart + 4*i - 4, mu, c2*tau);
                        }
                  }
                  
                  /******** Step exp(d1*tau*L_B) ********/
                  exp_tau_LB(d1*tau, X_cart);
            }
            
            /******** Step exp(2*c1*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, 2.*c1*tau);
            }
            
            /******** Averaging ********/
            for (i = 1; i <= N; i ++){
                  X_buff[4*i - 3] = X_cart[4*i - 3];  X_buff[4*i - 2] = X_cart[4*i - 2];  X_buff[4*i - 1] = X_cart[4*i - 1];  X_buff[4*i] = X_cart[4*i];
            }
            for (i = 1; i <= N; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_buff + 4*i - 4, mu, -c1*tau);
            }
            for (i = 1; i <= N; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
                  beta = m0*masses[i]/(m0 + masses[i]);
                  mu   = G*(m0 + masses[i]);
                  cart2ell(X_buff + 4*i - 4, alkhqp, mu);
                  e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
                  lbd  = continuousAngle(alkhqp[1],                  lbdOld[i]);
                  g    = continuousAngle(-atan2(alkhqp[3], alkhqp[2]), gOld[i]);
                  X_old[4*i - 3] = lbd;
                  X_old[4*i - 2] = g;
                  X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
                  X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
                  lbdOld[i]      = lbd;  gOld[i] = g;
            }
            if (iter % 2 == 0){
                  old2new(X_old, X_new, X_uv_2);
                  /******** Accumulating the Fourier coefficients ********/
                  t0 = tau* (typ) (iter - 2);
                  t1 = tau* (typ) (iter - 1);
                  t2 = tau* (typ) iter;
                  H0 = Cp*fast_pow(1. + cos(M_PI*t0/T), Hanning_order);
                  H1 = Cp*fast_pow(1. + cos(M_PI*t1/T), Hanning_order);
                  H2 = Cp*fast_pow(1. + cos(M_PI*t2/T), Hanning_order);
                  /******** Getting the amplitude of the fast terms ********/
                  for (j = 0; j <= how_many_harmonics; j ++){
                        nu = nu_fast * (typ) j;
                        for (i = 1; i <= N; i ++){
                              As[j*(1 + 4*N) + 4*i - 3] += tau*(H0*X_uv_0[4*i - 3]*cos(nu*t0) + 4.*H1*X_uv_1[4*i - 3]*cos(nu*t1) + H2*X_uv_2[4*i - 3]*cos(nu*t2))/(6.*T);
                              As[j*(1 + 4*N) + 4*i - 2] += tau*(H0*X_uv_0[4*i - 2]*cos(nu*t0) + 4.*H1*X_uv_1[4*i - 2]*cos(nu*t1) + H2*X_uv_2[4*i - 2]*cos(nu*t2))/(6.*T);
                              As[j*(1 + 4*N) + 4*i - 1] += tau*(H0*X_uv_0[4*i - 1]*cos(nu*t0) + 4.*H1*X_uv_1[4*i - 1]*cos(nu*t1) + H2*X_uv_2[4*i - 1]*cos(nu*t2))/(6.*T);
                              As[j*(1 + 4*N) + 4*i]     += tau*(H0*X_uv_0[4*i]    *cos(nu*t0) + 4.*H1*X_uv_1[4*i]    *cos(nu*t1) + H2*X_uv_2[4*i]    *cos(nu*t2))/(6.*T);
                        }
                  }
                  /******** Getting the amplitude of the slow term ********/
                  nu = nu_reso;
                  for (i = 1; i <= N; i ++){
                        as[4*i - 3] += tau*(H0*X_uv_0[4*i - 3]*cos(nu*t0) + 4.*H1*X_uv_1[4*i - 3]*cos(nu*t1) + H2*X_uv_2[4*i - 3]*cos(nu*t2))/(6.*T);
                        as[4*i - 2] += tau*(H0*X_uv_0[4*i - 2]*cos(nu*t0) + 4.*H1*X_uv_1[4*i - 2]*cos(nu*t1) + H2*X_uv_2[4*i - 2]*cos(nu*t2))/(6.*T);
                        as[4*i - 1] += tau*(H0*X_uv_0[4*i - 1]*cos(nu*t0) + 4.*H1*X_uv_1[4*i - 1]*cos(nu*t1) + H2*X_uv_2[4*i - 1]*cos(nu*t2))/(6.*T);
                        as[4*i]     += tau*(H0*X_uv_0[4*i]    *cos(nu*t0) + 4.*H1*X_uv_1[4*i]    *cos(nu*t1) + H2*X_uv_2[4*i]    *cos(nu*t2))/(6.*T);
                  }
                  for (i = 1; i <= 4*N; i ++){
                        X_uv_0[i] = X_uv_2[i];
                  }
            }
            else{
                  old2new(X_old, X_new, X_uv_1);
            }
      }
      
      /******** Evaluating at t = 0 and updating the averages (Used by LibrationCenterFollow for outputting) ********/
      for (i = 1; i <= 4*N; i ++){
            avgs[i] = As[i];
            As[i]  /= 2.;
      }
      for (j = 0; j <= how_many_harmonics; j ++){
            for (i = 1; i <= N; i ++){
                  X_uv[4*i - 3] += 2.*As[j*(1 + 4*N) + 4*i - 3];
                  X_uv[4*i - 2] += 2.*As[j*(1 + 4*N) + 4*i - 2];
                  X_uv[4*i - 1] += 2.*As[j*(1 + 4*N) + 4*i - 1];
                  X_uv[4*i]     += 2.*As[j*(1 + 4*N) + 4*i];
            }
      }
      
      /******** Evaluating amplitude of slow mode. Using first harmonic of f(t) = sqrt(sum_j (u_j^2 + v_j^2)) ********/
      toBeReturnedNum = 0.;
      for (i = 1; i <= N; i ++){
            toBeReturnedNum += as[4*i - 2]*as[4*i - 2] + as[4*i]*as[4*i];
      }
      toBeReturnedNum = sqrt(toBeReturnedNum);
      return toBeReturnedNum;
}


void Renormalization(typ * X_old){

      /******** Renormalizes a1 to 1 and lbd1 to 0 ********/

      int i;
      typ a1, lbd1, Lbd_renorm;

      a1 = X_old[4*1 - 1]*X_old[4*1 - 1]*(m0 + masses[1])/(G*m0*m0*masses[1]*masses[1]);
      Lbd_renorm = sqrt(a1);
      lbd1 = X_old[1];
      for (i = 1; i <= how_many_planet; i ++){
            //X_old[4*i - 3] -= lbd1;
            //X_old[4*i - 2] += lbd1;
            X_old[4*i - 1] /= Lbd_renorm;
            X_old[4*i]     /= Lbd_renorm;
      }
}


void get_n(typ * n){

      /******** Stores the mean motions in n ********/

      int i;
      typ mu_i, a_i, n_i;
      typ X_old[4*how_many_planet + 1];
      typ X_new[4*how_many_planet + 1];
      new2old(X_old, X_new, avgs);
      for (i = 1; i <= how_many_planet; i ++){
            mu_i = G*(m0 + masses[i]);
            a_i  = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            n_i  = sqrt(mu_i/(a_i*a_i*a_i));
            *(n + i) = n_i;
      }
}


void ConstantParameter(typ * X_new, typ * X_uv){

      /******** Similar to function nonDofReinit but for the non-averaged system ********/
      /******** The action variables conjugated to fast angles are not constant  ********/
      /******** and cannot be reinitialized after each iteration. I don't        ********/
      /******** reinitialize them but make sure that their ratio with the total  ********/
      /******** angular momentum G is constant by updating G accordingly         ********/
      /******** All angles are left untouched by this function                   ********/
      
      /* To be generalized to more complicated chains later */
      
      int i, k;
      typ AM0, Gamma0, param, AM, Gamma;
      
      /*AM0    = X_uv_t0[4*subchain[how_many_resonant]     - 1];
      Gamma0 = X_uv_t0[4*subchain[how_many_resonant - 1] - 1];
      param  = Gamma0/AM0;
      AM     = X_uv   [4*subchain[how_many_resonant]     - 1];
      Gamma  = X_uv   [4*subchain[how_many_resonant - 1] - 1];
      AM     = Gamma/param;
      
      X_uv[4*subchain[how_many_resonant] - 1] = AM;*/
      
      X_uv [4*subchain[how_many_resonant] - 1] = X_uv_t0 [4*subchain[how_many_resonant] - 1];
      X_uv [4*subchain[how_many_resonant] - 3] = X_uv_t0 [4*subchain[how_many_resonant] - 3];
      X_new[4*subchain[how_many_resonant] - 1] = X_new_t0[4*subchain[how_many_resonant] - 1];
      X_new[4*subchain[how_many_resonant] - 3] = X_new_t0[4*subchain[how_many_resonant] - 3];
      
      if      (how_many_resonant >= 2){
            
      }
      else if (how_many_resonant == 1){
      
      }
}


//void LibrationCenterFind(typ * X_old, int verbose){

      /******** Finds a libration center of the complete Hamiltonian by iteratively ********/
      /******** integrating and keeping only the average and fast frequency.        ********/
      /******** verbose is a boolean that specifies if the function should talk     ********/

      /*int i, j;
      int fast = subchain[how_many_resonant - 1];
      int slow = subchain[how_many_resonant];
      typ amplitudeRatio;
      typ X_new[4*how_many_planet + 1];
      typ X_uv [4*how_many_planet + 1];
      typ xvXu [4*how_many_planet + 1];
      typ tau = 0.0078125;
      typ dt[4] = {2., 1., 0.5, 0.5};            //Timestep in units of tau
      typ T [4] = {3000., 4500., 4500., 5500.};  //Integration time
      int Hf[4] = {2, 2, 5, 5};                  //order of Hanning filter
      int Hr[4] = {15, 25, 50, 80};              //Number of harmonics
      int Sn[4] = {1, 1, 1, 1};                  //Order of the SABA integrator
      
      if (verbose){
            printf("-------------------------------------------------------------------------------------------\n\n");
            printf("Starting the search for a libration center.\n\n");
            PointPrint(X_old, 0);
      }*/

      /******** Converging towards a libration center via the NAFF method ********/
      /*for (j = 0; j < 4; j ++){
            amplitudeRatio = UnaveragedSABAn_NAFF(tau*dt[j], T[j], Hf[j], xvXu, X_old, Sn[j], Hr[j]);
            ConstantParameter(X_new, xvXu);
            for (i = 1; i <= how_many_planet; i ++){
                  X_uv[4*i - 3] = xvXu[4*i - 3];
                  X_uv[4*i - 2] = xvXu[4*i - 2];
                  X_uv[4*i - 1] = xvXu[4*i - 1];
                  X_uv[4*i]     = xvXu[4*i];
            }
            new2old(X_old, X_new, X_uv);
            if (verbose){
                  PointPrint(X_old, j + 1);
                  printf("nu_%d = dphi_%d/dt = %.8lf, nu_%d = dphi_%d/dt = %.8lf, nu_%d/nu_%d = %.8lf\n", fast, fast, nu_fast, slow, slow, nu_reso, fast, slow, nu_fast/nu_reso);
                  printf("Amplitude ratio = %.10lf\n\n", amplitudeRatio);
                  //UnaveragedSABAn(tau/2., 4000., 2, X_old, 4); //To be removed
                  //new2old(X_old, X_new, X_uv);                 //To be removed
            }
      }

      UnaveragedSABAn(tau/2., 4000., 2, X_old, 5);
      new2old(X_old, X_new, X_uv);
}*/


void LibrationCenterFind(typ * X_old, int precision){

      /******** Finds a libration center of the complete Hamiltonian by iteratively ********/
      /******** integrating and keeping only the average and fast frequency.        ********/
      /******** precision = {0, 1, 2} -> {low, medium, high} precision.             ********/

      int i, j, went2fixedPoint, noProgress;
      int fast = subchain[how_many_resonant - 1];
      int slow = subchain[how_many_resonant];
      typ amplitude = 1.;
      typ oldAmplitude, mean_e, P;
      typ X_new[4*how_many_planet + 1];
      typ X_uv [4*how_many_planet + 1];
      typ xvXu [4*how_many_planet + 1];
      typ n    [  how_many_planet + 1];
      typ tau = 0.0078125;
      
      typ dt[3] = {2., 0.75, 0.5};         //Timestep in units of tau
      typ T [3] = {3000., 4500., 4500.};   //Integration time
      int Hf[3] = {2, 5, 5};               //order of Hanning filter
      int Hr[3] = {15, 45, 75};            //Number of harmonics
      int Sn[3] = {1, 1, 4};               //Order of the SABA integrator
      typ AR[3] = {0.5e-2, 0.5e-4, 2.e-6}; //Required value for the amplitude
      
      /******** Obtaining the average value of the eccentricity in order to adapt the required amplitude for convergence ********/
      mean_e = 0.;
      for (i = 1; i <= how_many_planet; i ++){
            mean_e += sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
      }
      mean_e /= (typ) how_many_planet;
      mean_e  = max(1.e-5, mean_e);
      AR[0] *= mean_e;  AR[1] *= mean_e;  AR[2] *= mean_e;
      
      /******** Adapting the length of the integration according to the largest fundamental period ********/
      if (nu_reso != 0. && nu_fast != 0.){
            P    = max(fabs(2.*M_PI/nu_reso), fabs(2.*M_PI/nu_fast));
            T[0] = 4.*P;  T[1] = 8.*P;  T[2] = 16.*P;
      }
      
      if (precision < 0 || precision > 2){
            fprintf(stderr, "\nError: The integer precision must be 0, 1 or 2 in function LibrationCenterFind\n");
            abort();
      }
      
      /******** To be removed ********/
      /*typ dt[4] = {2., 1., 0.5, 0.5};            //Timestep in units of tau
      typ T [4] = {3000., 4500., 4500., 5500.};  //Integration time
      int Hf[4] = {2, 2, 5, 5};                  //order of Hanning filter
      int Hr[4] = {15, 25, 50, 80};              //Number of harmonics
      int Sn[4] = {1, 1, 1, 1};                  //Order of the SABA integrator*/
      
      
      printf("-------------------------------------------------------------------------------------------\n\n");
      printf("Starting the search for a libration center.\n\n");
      PointPrint(X_old, 0);

      /******** Converging towards a libration center ********/
      j = 1;  oldAmplitude = 1.e300;  went2fixedPoint = 0;  noProgress = 0;
      while (amplitude > AR[precision]){
            amplitude = UnaveragedSABAn_NAFF(tau*dt[precision], T[precision], Hf[precision], xvXu, X_old, Sn[precision], Hr[precision] + 10*(how_many_planet - 2));
            ConstantParameter(X_new, xvXu);
            for (i = 1; i <= how_many_planet; i ++){
                  X_uv[4*i - 3] = xvXu[4*i - 3];
                  X_uv[4*i - 2] = xvXu[4*i - 2];
                  X_uv[4*i - 1] = xvXu[4*i - 1];
                  X_uv[4*i]     = xvXu[4*i];
            }
            new2old(X_old, X_new, X_uv);
            if (amplitude > AR[precision]){
                  PointPrint(X_old, j);
                  get_n(n);
                  printf("nu_%d = dphi_%d/dt = %.14lf, nu_%d = dphi_%d/dt = %.14lf, nu_%d/nu_%d = %.14lf\n", fast, fast, nu_fast, slow, slow, nu_reso, fast, slow, nu_fast/nu_reso);
                  for (i = 1; i < how_many_planet; i ++){
                        printf("n_%d/n_%d = %.8lf", i, i + 1, n[i]/n[i + 1]);
                        if (i < how_many_planet - 1){printf(", ");} else{printf("\n");}
                  }
                  printf("Amplitude = %.20lf, required = %.13lf\n\n", amplitude, AR[precision]);
            }
            j ++;
            if (amplitude > 0.9*oldAmplitude && amplitude < 1.4*oldAmplitude){
                  if (noProgress < 3){
                        noProgress ++;
                  }
                  else{
                        Hr[precision] += 10;
                        T [precision] *= 2.;
                        noProgress     = 0;
                  }
            }
            if (j > 32 && !went2fixedPoint ){
                  printf("Cannot converge. Starting back from a fixed point.\n");
                  EquilibriumFind(X_old, precision);
                  went2fixedPoint = 1;
            }
            if (j > 64){
                  fprintf(stderr, "\nError: Cannot converge in function LibrationCenterFind.\n");
                  abort();
            }
            oldAmplitude = amplitude;
      }
      Renormalization(X_old);
      PointPrint(X_old, j - 1);
      get_n(n);
      printf("nu_%d = dphi_%d/dt = %.14lf, nu_%d = dphi_%d/dt = %.14lf, nu_%d/nu_%d = %.14lf\n", fast, fast, nu_fast, slow, slow, nu_reso, fast, slow, nu_fast/nu_reso);
      for (i = 1; i < how_many_planet; i ++){
            printf("n_%d/n_%d = %.8lf", i, i + 1, n[i]/n[i + 1]);
            if (i < how_many_planet - 1){printf(", ");} else{printf("\n");}
      }
      printf("Amplitude = %.20lf, required = %.13lf\n\n", amplitude, AR[precision]);
      /******** To be removed ********/
      //UnaveragedSABAn(tau/2., 4000., 2, X_old, 5);
      //new2old(X_old, X_new, X_uv);
}


void LibrationCenterFollow(typ * X_old, typ dG, int Npoints, int precision){

      /******** Finds a libration center and follows its family by updating the total    ********/
      /******** angular momentum G by an amount dG before each search. Npoints libration ********/
      /******** centers are found and stored to the file pth/LibrationCenters_chain.txt  ********/
      /******** precision = {0, 1, 2} -> {low, medium, high} precision.                  ********/

      int i, j, dG_inc_count, dG_dec_count;
      int fast = subchain[how_many_resonant - 1];
      int slow = subchain[how_many_resonant];
      typ a, e, sig, lbd, e_av, a_av, oldRatio, newRatio;
      char file_path[800];
      char nn[10];
      FILE * file;
      typ X_new    [4*how_many_planet + 1];
      typ X_uv     [4*how_many_planet + 1];
      typ X_old_av [4*how_many_planet + 1];
      typ X_new_av [4*how_many_planet + 1];
      typ sigOld   [  how_many_planet + 1];

      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "LibrationCenters");
      strcat(file_path, "_");
      for (i = 1; i <= how_many_planet; i ++){
            sprintf(nn, "%d", p_i[i]);
            strcat(file_path, nn);
      }
      strcat(file_path, ".txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "\nError: Cannot create or open file LibrationCenters_chain.txt in function LibrationCenterFollow.\n");
            abort();
      }
      
      /******** Writing to file ********/
      fprintf(file, "This file contains a family of libration centers of the resonance chain ");
      for (i = 1; i < how_many_planet; i ++){
            fprintf(file, "%d:", p_i[i]);
      }
      fprintf(file, "%d", p_i[how_many_planet]);
      fprintf(file, " in the unaveraged problem.\n");
      fprintf(file, "The family is parameterized by Phi_%d.\n", slow);
      fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 4 + 8*how_many_planet, how_many_planet);
      fprintf(file, "nu_%d = dphi_%d/dt, nu_%d = dphi_%d/dt, Phi_%d, <Phi_%d>, a_j, e_j, phi_j, sig_j, <a_j>, <e_j>, <phi_j>, <sig_j>\n", fast, fast, slow, slow, slow, slow);
      fprintf(file, "\n");

      /******** Finding the first libration center ********/
      LibrationCenterFind(X_old, precision);
      /******** Writing to file and initializing sigOld ********/
      old2new(X_old,    X_new,    X_uv);
      printf("Phi_%d = %.20lf, <Phi_%d> = %.20lf\n", slow, X_uv[4*slow - 1], slow, avgs[4*slow - 1]);
      new2old(X_old_av, X_new_av, avgs);
      fprintf(file, "%.8lf %.8lf %.20lf %.20lf", nu_fast, nu_reso, X_uv[4*slow - 1], avgs[4*slow - 1]);
      for (i = 1; i <= how_many_planet; i ++){
            sigOld[i] = X_new[4*i - 2];
            sig       = X_new[4*i - 2];
            e         = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
            a         = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            e_av      = sqrt(1. - (1. - X_old_av[4*i]/X_old_av[4*i - 1])*(1. - X_old_av[4*i]/X_old_av[4*i - 1]));
            a_av      = X_old_av[4*i - 1]*X_old_av[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
            fprintf(file, " %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf", a, e, X_uv[4*i - 3], sig, a_av, e_av, avgs[4*i - 3], X_new_av[4*i - 2]);
      }
      fprintf(file, "\n");
      
      /******** Updating the total angular momentum ********/
      X_uv    [4*slow - 1] += dG;
      X_uv_t0 [4*slow - 1]  = X_uv [4*slow - 1];
      X_new_t0[4*slow - 1]  = X_new[4*slow - 1];
      new2old(X_old, X_new, X_uv);

      /******** Displaying progress ********/
      printf("progress = %d/%d\n", 1, Npoints);
      
      /******** Following the family of libration centers ********/
      oldRatio = nu_fast/nu_reso;  newRatio = oldRatio;  dG_inc_count = 0;  dG_dec_count = 0;
      for (j = 1; j < Npoints; j ++){
            LibrationCenterFind(X_old, precision);
            old2new(X_old,    X_new,    X_uv);
            printf("Phi_%d = %.20lf, <Phi_%d> = %.20lf\n", slow, X_uv[4*slow - 1], slow, avgs[4*slow - 1]);
            new2old(X_old_av, X_new_av, avgs);
            /******** Increasing dG if travelling too slowly. Decreasing it if travelling too fast ********/
            newRatio = nu_fast/nu_reso;
            /*if      (fabs(oldRatio - newRatio) < 0.2 && dG_inc_count < 10){
                  dG *= 1.05;
                  printf("dG = %.18lf\n", dG);
                  dG_inc_count ++;
            }
            else if (fabs(oldRatio - newRatio) > 2. && dG_dec_count < 10){
                  dG *= 0.95;
                  printf("dG = %.18lf\n", dG);
                  dG_dec_count ++;
            }*/
            oldRatio = newRatio;
            /******** Writing to file ********/
            fprintf(file, "%.8lf %.8lf %.20lf %.20lf", nu_fast, nu_reso, X_uv[4*slow - 1], avgs[4*slow - 1]);
            for (i = 1; i <= how_many_planet; i ++){
                  sig  = continuousAngle(X_new[4*i - 2], sigOld[i]);
                  e    = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
                  a    = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                  e_av = sqrt(1. - (1. - X_old_av[4*i]/X_old_av[4*i - 1])*(1. - X_old_av[4*i]/X_old_av[4*i - 1]));
                  a_av = X_old_av[4*i - 1]*X_old_av[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                  fprintf(file, " %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf", a, e, X_uv[4*i - 3], sig, a_av, e_av, avgs[4*i - 3], X_new_av[4*i - 2]);
                  sigOld[i] = sig;
            }
            fprintf(file, "\n");
            /******** Updating the total angular momentum ********/
            if (j + 1 != Npoints){
                  X_uv    [4*slow - 1] += dG;
                  X_uv_t0 [4*slow - 1]  = X_uv [4*slow - 1];
                  X_new_t0[4*slow - 1]  = X_new[4*slow - 1];
                  new2old(X_old, X_new, X_uv);
            }
            /******** Displaying progress ********/
            printf("progress = %d/%d\n", j + 1, Npoints);
      }

      /******** Closing output file ********/
      fclose(file);
}


void PeriodicOrbitFind(typ * X_old){

      /******** Looks for a periodic orbit of the N-body problem.                               ********/
      /******** Finds a libration center, that is, a quasi-periodic orbit with two frequencies. ********/
      /******** Then changes the relevant first integral until the ratio between the two        ********/
      /******** frequencies is a rational number a/b with b and a - b not too large             ********/
      
      int i, j, targetDenom, precision, count;
      int fast = subchain[how_many_resonant - 1];
      int slow = subchain[how_many_resonant];
      typ X_new [4*how_many_planet + 1];
      typ X_uv  [4*how_many_planet + 1];
      typ dG, current_ratio, old_ratio, target, t, T1, T2, T, dt;
      struct rational target_ratio;
      typ epsilon = 0.;
      //typ toTarget[3] = {0.01, 0.0001, 0.0000001}; //The maximal distance between the current frequency ratio and the target frequency ratio once converged.
      typ toTarget[3] = {0.01, 0.0001, 1.e-11}; //The maximal distance between the current frequency ratio and the target frequency ratio once converged.
      
      for (i = 1; i <= how_many_planet; i ++){
            epsilon += masses[i]/m0;
      }

      /******** Finding the first libration center ********/
      LibrationCenterFind(X_old, 0);
      old2new(X_old, X_new, X_uv);
      
      /******** Trying first value for dG ********/
      dG            = 0.00002*epsilon;
      current_ratio = nu_fast/nu_reso;
      old_ratio     = current_ratio;
      
      /******** Obtaining the target frequency ratio ********/
      target_ratio = int2rat(1000.);
      targetDenom  = 1;
      while (fabs(rat2real(target_ratio) - current_ratio) > 0.05){
            target_ratio = real2rat(current_ratio, targetDenom);
            targetDenom ++;
      }
      printf("\nTarget value for nu_%d/nu_%d : ", fast, slow);
      ratprint(target_ratio);
      printf("\n");
      
      /******** Updating the total angular momentum ********/
      X_uv    [4*slow - 1] += dG;
      X_uv_t0 [4*slow - 1]  = X_uv [4*slow - 1];
      X_new_t0[4*slow - 1]  = X_new[4*slow - 1];
      new2old(X_old, X_new, X_uv);
      
      /******** Following the family of libration centers ********/
      target = rat2real(target_ratio);
      for (precision = 0; precision < 3; precision ++){
            if (precision > 0){
                  dG /= 2.;
            }
            count = 0;
            while (fabs(current_ratio - target) > toTarget[precision]){ //Obtaining a small precision
                  LibrationCenterFind(X_old, precision);
                  current_ratio = nu_fast/nu_reso;
                  if (fabs(current_ratio - target) > fabs(old_ratio - target) && (current_ratio - target)*(old_ratio - target) > 0. && count){ //Going the wrong direction
                        dG *= -0.95;
                  }
                  else if ((current_ratio - target)*(old_ratio - target) < 0.){ //Went too far
                        t   = (target - old_ratio)/(current_ratio - old_ratio);
                        dG *= -(1. - t);
                  }
                  else if (fabs(current_ratio - target)/fabs(old_ratio - target) > 0.9 && fabs(current_ratio - target)/fabs(old_ratio - target) < 1. && count){ //Going too slowly
                        dG *= 1.2;
                  }
                  /******** Updating the total angular momentum ********/
                  if (fabs(current_ratio - target) > toTarget[precision]){
                        X_uv    [4*slow - 1] += dG;
                        X_uv_t0 [4*slow - 1]  = X_uv [4*slow - 1];
                        X_new_t0[4*slow - 1]  = X_new[4*slow - 1];
                        new2old(X_old, X_new, X_uv);
                  }
                  old_ratio = current_ratio;
                  count ++;
            }
      }
      /******** Numerically integrating the periodic orbit found ********/
      T1 = min(target_ratio.numerator, target_ratio.denominator)*max(2.*M_PI/nu_fast, 2.*M_PI/nu_reso);
      T2 = max(target_ratio.numerator, target_ratio.denominator)*min(2.*M_PI/nu_fast, 2.*M_PI/nu_reso);
      T  = 0.5*(fabs(T1) + fabs(T2));
      dt = T/(256.*round(T));
      UnaveragedSABAn(dt, 3.2*T, 1, X_old, 4);
}







