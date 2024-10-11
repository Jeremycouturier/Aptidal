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



void ell2cart(typ a, typ e, typ i, typ nu, typ varpi, typ Omega, typ mu, typ * cart){

      /******** Updates the array [X,Y,Z,vX,vY,vZ] of the cartesian coordinates. mu is G*M                    ********/
      /******** a is the semi-major axis, e is the eccentricity, i is the inclination, nu is the true anomaly,********/
      /******** omega is the longitude of periapsis and Omega is the longitude of the ascending node          ********/
      
      typ X,Y,vX,vY;                //Cartesian coordinates
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
      X  =  X_buff*(pp*cosvarpi + dpq*sinvarpi) +  Y_buff*(dpq*cosvarpi - pp*sinvarpi);
      vX = vX_buff*(pp*cosvarpi + dpq*sinvarpi) + vY_buff*(dpq*cosvarpi - pp*sinvarpi);
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
            if (fabs(H1 - H2) + fabs(K1 - K2) > 1.0e-6){
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
            fprintf(stderr, "Error : Wrong computation of the integration time in function mean2true.\n");
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
      NMAX = 15;

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
            typ EPSQ = 1.5*EPS;
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
                        fprintf(stderr, "Error: Fatal error in function newt: X = %lf, DIFF = %lf, i = %d\n", X, DIFF, i);
                        abort();
                  }
                  else{
                        printf("Error: Fatal error in function newt: X = %lf, DIFF = %lf, i = %d\n", X, DIFF, i);
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
      newt(dm, ce, se, &delta, &c, &s, &e, &cm1, &smx, 1);

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


void UnaveragedSABA1(typ tau, typ T, int output_step, typ * X_old){

      /******** Integrates the unaveraged Hamiltonian with a SABA1 method   ********/
      /******** (Leapfrog) for a time T with a timestep tau. Outputs every  ********/
      /******** output_step timestep to file pth/UnaveragedSABA1.txt.       ********/
      /******** The Keplerian part A is integrated exactly, but the         ********/
      /******** perturbative part takes the form B = B_1(p) + B_2(q) and is ********/
      /******** integrated approximetaly but symplectically with a SABA1.   ********/


      char file_path[800];
      FILE * file;
      int N_step, iter, i;
      typ X_cart[4*how_many_planet + 1];
      typ X_buff[4*how_many_planet + 1];
      typ X_new [4*how_many_planet + 1];
      typ X_uv  [4*how_many_planet + 1];
      typ a, e, sig, vp, M, mu, nu, beta, H;
      typ alkhqp[6];
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "UnaveragedSABA1.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file UnaveragedSABA1.txt in function UnaveragedSABA1.\n");
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
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
      }

      /******** Integrating ********/
      N_step = (int) ceil(T/tau);
      for (iter = 0; iter < N_step; iter ++){
      
            /******** Writing to file ********/
            if (iter%output_step == 0){
                  for (i = 1; i <= how_many_planet; i ++){
                        X_buff[4*i - 3] = X_cart[4*i - 3]; X_buff[4*i - 2] = X_cart[4*i - 2]; X_buff[4*i - 1] = X_cart[4*i - 1]; X_buff[4*i] = X_cart[4*i];
                  }
                  if (iter){ //Need to perform a step exp(-tau/2*L_A) on the buffer before outputting
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_buff + 4*i - 4, mu, -tau/2.);
                        }
                  }
                  else{
                        fprintf(file, "Numerical integration of the unaveraged Hamiltonian with a SABA1 integrator in the cartesian coordinates.\n");
                        fprintf(file, "The Keplerian part A is integrated exactly using function kepsaut. The perturbative part B is not integrable but can be\n");
                        fprintf(file, "written B = B1 + B2 with B1 and B2 both integrable. Therefore, B is integrated approximately but symplectically with a SABA1.\n");
                        fprintf(file, "\n");
                        fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 2 + 7*how_many_planet, how_many_planet);
                        fprintf(file, "Time, Hamiltonian, phi_j, v_j, Phi_j, u_j, a_j, e_j, sig_j.\n");
                        fprintf(file, "\n");
                  }
                  for (i = 1; i <= how_many_planet; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
                        beta = m0*masses[i]/(m0 + masses[i]);
                        mu   = G*(m0 + masses[i]);
                        cart2ell(X_buff + 4*i - 4, alkhqp, mu);
                        e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
                        X_old[4*i - 3] = alkhqp[1];
                        X_old[4*i - 2] = -atan2(alkhqp[3], alkhqp[2]);
                        X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
                        X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
                  }
                  old2new(X_old, X_new, X_uv);             // (lbd_j, -vrp_j; Lbd_j, D_j) -> (phi_j, v_j; Phi_j, u_j)
                  H = 0;
                  fprintf(file, "%.12lf %.12lf", tau*(typ) iter, H);
                  for (i = 1; i <= how_many_planet; i ++){
                        sig = atan2(X_uv[4*i - 2], X_uv[4*i]);
                        e   = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
                        a   = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                        fprintf(file, " %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", X_uv[4*i - 3], X_uv[4*i - 2], X_uv[4*i - 1], X_uv[4*i], a, e, sig);
                  }
                  fprintf(file, "\n");
            }
            
            /******** Step exp(tau/2*L_A). For the first iteration only ********/
            if (!iter){
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, tau/2.);
                  }
            }
            
            /******** Step exp(tau*L_B) ********/
            exp_tau_LB(tau, X_cart);
            
            /******** Step exp(tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, tau);
            }
      }

      /******** Closing output file ********/
      fclose(file);
}


void UnaveragedSABA2(typ tau, typ T, int output_step, typ * X_old){

      /******** Integrates the unaveraged Hamiltonian with a SABA2 method   ********/
      /******** for a time T with a timestep tau. Outputs every output_step ********/
      /******** timestep to file pth/UnaveragedSABA2.txt.                   ********/
      /******** The Keplerian part A is integrated exactly, but the         ********/
      /******** perturbative part takes the form B = B_1(p) + B_2(q) and is ********/
      /******** integrated approximetaly but symplectically with a SABA1.   ********/


      char file_path[800];
      FILE * file;
      int N_step, iter, i;
      typ X_cart[4*how_many_planet + 1];
      typ X_buff[4*how_many_planet + 1];
      typ X_new [4*how_many_planet + 1];
      typ X_uv  [4*how_many_planet + 1];
      typ a, e, sig, vp, M, mu, nu, beta, H;
      typ c_2 = 1./sqrt(3.);
      typ c_1 = 0.5*(1. - c_2);
      typ d_1 = 0.5;
      typ alkhqp[6];
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "UnaveragedSABA2.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file UnaveragedSABA2.txt in function UnaveragedSABA2.\n");
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
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
      }

      /******** Integrating ********/
      N_step = (int) ceil(T/tau);
      for (iter = 0; iter < N_step; iter ++){
      
            /******** Writing to file ********/
            if (iter%output_step == 0){
                  for (i = 1; i <= how_many_planet; i ++){
                        X_buff[4*i - 3] = X_cart[4*i - 3]; X_buff[4*i - 2] = X_cart[4*i - 2]; X_buff[4*i - 1] = X_cart[4*i - 1]; X_buff[4*i] = X_cart[4*i];
                  }
                  if (iter){ //Need to perform a step exp(-c_1*tau*L_A) on the buffer before outputting
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_buff + 4*i - 4, mu, -c_1*tau);
                        }
                  }
                  else{
                        fprintf(file, "Numerical integration of the unaveraged Hamiltonian with a SABA2 integrator in the cartesian coordinates.\n");
                        fprintf(file, "The Keplerian part A is integrated exactly using function kepsaut. The perturbative part B is not integrable but can be\n");
                        fprintf(file, "written B = B1 + B2 with B1 and B2 both integrable. Therefore, B is integrated approximately but symplectically with a SABA1.\n");
                        fprintf(file, "\n");
                        fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 2 + 7*how_many_planet, how_many_planet);
                        fprintf(file, "Time, Hamiltonian, phi_j, v_j, Phi_j, u_j, a_j, e_j, sig_j.\n");
                        fprintf(file, "\n");
                  }
                  for (i = 1; i <= how_many_planet; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
                        beta = m0*masses[i]/(m0 + masses[i]);
                        mu   = G*(m0 + masses[i]);
                        cart2ell(X_buff + 4*i - 4, alkhqp, mu);
                        e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
                        X_old[4*i - 3] = alkhqp[1];
                        X_old[4*i - 2] = -atan2(alkhqp[3], alkhqp[2]);
                        X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
                        X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
                  }
                  old2new(X_old, X_new, X_uv);             // (lbd_j, -vrp_j; Lbd_j, D_j) -> (phi_j, v_j; Phi_j, u_j)
                  H = 0;
                  fprintf(file, "%.12lf %.12lf", tau*(typ) iter, H);
                  for (i = 1; i <= how_many_planet; i ++){
                        sig = atan2(X_uv[4*i - 2], X_uv[4*i]);
                        e   = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
                        a   = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                        fprintf(file, " %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", X_uv[4*i - 3], X_uv[4*i - 2], X_uv[4*i - 1], X_uv[4*i], a, e, sig);
                  }
                  fprintf(file, "\n");
            }
            
            /******** Step exp(c_1*tau*L_A). For the first iteration only ********/
            if (!iter){
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c_1*tau);
                  }
            }
            
            /******** Step exp(d_1*tau*L_B) ********/
            exp_tau_LB(d_1*tau, X_cart);
            
            /******** Step exp(c_2*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, c_2*tau);
            }

            /******** Step exp(d_1*tau*L_B) ********/
            exp_tau_LB(d_1*tau, X_cart);
            
            /******** Step exp(2*c_1*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, 2.*c_1*tau);
            }
      }

      /******** Closing output file ********/
      fclose(file);
}


void UnaveragedSABA3(typ tau, typ T, int output_step, typ * X_old){

      /******** Integrates the unaveraged Hamiltonian with a SABA3 method   ********/
      /******** for a time T with a timestep tau. Outputs every output_step ********/
      /******** timestep to file pth/UnaveragedSABA3.txt.                   ********/
      /******** The Keplerian part A is integrated exactly, but the         ********/
      /******** perturbative part takes the form B = B_1(p) + B_2(q) and is ********/
      /******** integrated approximetaly but symplectically with a SABA1.   ********/


      char file_path[800];
      FILE * file;
      int N_step, iter, i;
      typ X_cart[4*how_many_planet + 1];
      typ X_buff[4*how_many_planet + 1];
      typ X_new [4*how_many_planet + 1];
      typ X_uv  [4*how_many_planet + 1];
      typ a, e, sig, vp, M, mu, nu, beta, H;
      typ c_2 = sqrt(15.)/10.;
      typ c_1 = 0.5 - c_2;
      typ d_1 = 5./18.;
      typ d_2 = 4./9.;
      typ alkhqp[6];
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "UnaveragedSABA3.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file UnaveragedSABA3.txt in function UnaveragedSABA3.\n");
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
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
      }

      /******** Integrating ********/
      N_step = (int) ceil(T/tau);
      for (iter = 0; iter < N_step; iter ++){
      
            /******** Writing to file ********/
            if (iter%output_step == 0){
                  for (i = 1; i <= how_many_planet; i ++){
                        X_buff[4*i - 3] = X_cart[4*i - 3]; X_buff[4*i - 2] = X_cart[4*i - 2]; X_buff[4*i - 1] = X_cart[4*i - 1]; X_buff[4*i] = X_cart[4*i];
                  }
                  if (iter){ //Need to perform a step exp(-c_1*tau*L_A) on the buffer before outputting
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_buff + 4*i - 4, mu, -c_1*tau);
                        }
                  }
                  else{
                        fprintf(file, "Numerical integration of the unaveraged Hamiltonian with a SABA3 integrator in the cartesian coordinates.\n");
                        fprintf(file, "The Keplerian part A is integrated exactly using function kepsaut. The perturbative part B is not integrable but can be\n");
                        fprintf(file, "written B = B1 + B2 with B1 and B2 both integrable. Therefore, B is integrated approximately but symplectically with a SABA1.\n");
                        fprintf(file, "\n");
                        fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 2 + 7*how_many_planet, how_many_planet);
                        fprintf(file, "Time, Hamiltonian, phi_j, v_j, Phi_j, u_j, a_j, e_j, sig_j.\n");
                        fprintf(file, "\n");
                  }
                  for (i = 1; i <= how_many_planet; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
                        beta = m0*masses[i]/(m0 + masses[i]);
                        mu   = G*(m0 + masses[i]);
                        cart2ell(X_buff + 4*i - 4, alkhqp, mu);
                        e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
                        X_old[4*i - 3] = alkhqp[1];
                        X_old[4*i - 2] = -atan2(alkhqp[3], alkhqp[2]);
                        X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
                        X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
                  }
                  old2new(X_old, X_new, X_uv);             // (lbd_j, -vrp_j; Lbd_j, D_j) -> (phi_j, v_j; Phi_j, u_j)
                  H = 0;
                  fprintf(file, "%.12lf %.12lf", tau*(typ) iter, H);
                  for (i = 1; i <= how_many_planet; i ++){
                        sig = atan2(X_uv[4*i - 2], X_uv[4*i]);
                        e   = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
                        a   = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                        fprintf(file, " %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", X_uv[4*i - 3], X_uv[4*i - 2], X_uv[4*i - 1], X_uv[4*i], a, e, sig);
                  }
                  fprintf(file, "\n");
            }
            
            /******** Step exp(c_1*tau*L_A). For the first iteration only ********/
            if (!iter){
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c_1*tau);
                  }
            }
            
            /******** Step exp(d_1*tau*L_B) ********/
            exp_tau_LB(d_1*tau, X_cart);
            
            /******** Step exp(c_2*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, c_2*tau);
            }

            /******** Step exp(d_2*tau*L_B) ********/
            exp_tau_LB(d_2*tau, X_cart);
            
            /******** Step exp(c_2*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, c_2*tau);
            }

            /******** Step exp(d_1*tau*L_B) ********/
            exp_tau_LB(d_1*tau, X_cart);

            /******** Step exp(2*c_1*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, 2.*c_1*tau);
            }
      }

      /******** Closing output file ********/
      fclose(file);
}


void UnaveragedSABA4(typ tau, typ T, int output_step, typ * X_old){

      /******** Integrates the unaveraged Hamiltonian with a SABA4 method   ********/
      /******** for a time T with a timestep tau. Outputs every output_step ********/
      /******** timestep to file pth/UnaveragedSABA4.txt.                   ********/
      /******** The Keplerian part A is integrated exactly, but the         ********/
      /******** perturbative part takes the form B = B_1(p) + B_2(q) and is ********/
      /******** integrated approximetaly but symplectically with a SABA1.   ********/


      char file_path[800];
      FILE * file;
      int N_step, iter, i;
      typ X_cart[4*how_many_planet + 1];
      typ X_buff[4*how_many_planet + 1];
      typ X_new [4*how_many_planet + 1];
      typ X_uv  [4*how_many_planet + 1];
      typ a, e, sig, vp, M, mu, nu, beta, H;
      typ c_1 = 0.5 - sqrt(525. + 70.*sqrt(30.))/70.;
      typ c_2 = (sqrt(525. + 70.*sqrt(30.)) - sqrt(525. - 70.*sqrt(30.)))/70.;
      typ c_3 = sqrt(525. - 70.*sqrt(30.))/35.;
      typ d_1 = 0.25 - sqrt(30.)/72.;
      typ d_2 = 0.25 + sqrt(30.)/72.;
      typ alkhqp[6];
      
      /******** Opening output file ********/
      strcpy(file_path, pth);
      strcat(file_path, "UnaveragedSABA4.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file UnaveragedSABA4.txt in function UnaveragedSABA4.\n");
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
            ell2cart(a, e, 0., nu, vp, 0., mu, X_cart + 4*i - 4);
      }

      /******** Integrating ********/
      N_step = (int) ceil(T/tau);
      for (iter = 0; iter < N_step; iter ++){
      
            /******** Writing to file ********/
            if (iter%output_step == 0){
                  for (i = 1; i <= how_many_planet; i ++){
                        X_buff[4*i - 3] = X_cart[4*i - 3]; X_buff[4*i - 2] = X_cart[4*i - 2]; X_buff[4*i - 1] = X_cart[4*i - 1]; X_buff[4*i] = X_cart[4*i];
                  }
                  if (iter){ //Need to perform a step exp(-c_1*tau*L_A) on the buffer before outputting
                        for (i = 1; i <= how_many_planet; i ++){
                              mu = G*(m0 + masses[i]);
                              kepsaut(X_buff + 4*i - 4, mu, -c_1*tau);
                        }
                  }
                  else{
                        fprintf(file, "Numerical integration of the unaveraged Hamiltonian with a SABA4 integrator in the cartesian coordinates.\n");
                        fprintf(file, "The Keplerian part A is integrated exactly using function kepsaut. The perturbative part B is not integrable but can be\n");
                        fprintf(file, "written B = B1 + B2 with B1 and B2 both integrable. Therefore, B is integrated approximately but symplectically with a SABA1.\n");
                        fprintf(file, "\n");
                        fprintf(file, "This file has %d columns that are (for 1 <= j <= %d):\n", 2 + 7*how_many_planet, how_many_planet);
                        fprintf(file, "Time, Hamiltonian, phi_j, v_j, Phi_j, u_j, a_j, e_j, sig_j.\n");
                        fprintf(file, "\n");
                  }
                  for (i = 1; i <= how_many_planet; i ++){ // (x, y; vx, vy) -> (lbd_j, -vrp_j; Lbd_j, D_j)
                        beta = m0*masses[i]/(m0 + masses[i]);
                        mu   = G*(m0 + masses[i]);
                        cart2ell(X_buff + 4*i - 4, alkhqp, mu);
                        e    = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]);
                        X_old[4*i - 3] = alkhqp[1];
                        X_old[4*i - 2] = -atan2(alkhqp[3], alkhqp[2]);
                        X_old[4*i - 1] = beta*sqrt(mu*alkhqp[0]);
                        X_old[4*i]     = X_old[4*i - 1]*(1. - sqrt(1. - e*e));
                  }
                  old2new(X_old, X_new, X_uv);             // (lbd_j, -vrp_j; Lbd_j, D_j) -> (phi_j, v_j; Phi_j, u_j)
                  H = 0;
                  fprintf(file, "%.12lf %.12lf", tau*(typ) iter, H);
                  for (i = 1; i <= how_many_planet; i ++){
                        sig = atan2(X_uv[4*i - 2], X_uv[4*i]);
                        e   = sqrt(1. - (1. - X_old[4*i]/X_old[4*i - 1])*(1. - X_old[4*i]/X_old[4*i - 1]));
                        a   = X_old[4*i - 1]*X_old[4*i - 1]*(m0 + masses[i])/(G*m0*m0*masses[i]*masses[i]);
                        fprintf(file, " %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf", X_uv[4*i - 3], X_uv[4*i - 2], X_uv[4*i - 1], X_uv[4*i], a, e, sig);
                  }
                  fprintf(file, "\n");
            }
            
            /******** Step exp(c_1*tau*L_A). For the first iteration only ********/
            if (!iter){
                  for (i = 1; i <= how_many_planet; i ++){
                        mu = G*(m0 + masses[i]);
                        kepsaut(X_cart + 4*i - 4, mu, c_1*tau);
                  }
            }
            
            /******** Step exp(d_1*tau*L_B) ********/
            exp_tau_LB(d_1*tau, X_cart);
            
            /******** Step exp(c_2*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, c_2*tau);
            }

            /******** Step exp(d_2*tau*L_B) ********/
            exp_tau_LB(d_2*tau, X_cart);

            /******** Step exp(c_3*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, c_3*tau);
            }

            /******** Step exp(d_2*tau*L_B) ********/
            exp_tau_LB(d_2*tau, X_cart);
            
            /******** Step exp(c_2*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, c_2*tau);
            }

            /******** Step exp(d_1*tau*L_B) ********/
            exp_tau_LB(d_1*tau, X_cart);

            /******** Step exp(2*c_1*tau*L_A) ********/
            for (i = 1; i <= how_many_planet; i ++){
                  mu = G*(m0 + masses[i]);
                  kepsaut(X_cart + 4*i - 4, mu, 2.*c_1*tau);
            }
      }

      /******** Closing output file ********/
      fclose(file);
}
