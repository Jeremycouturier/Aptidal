#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include "parameters.h"
#include "transformation.h"
#include "coefficients.h"
#include "structure.h"
#include "calculus.h"


/******** Defining some global variables ********/
typ masses  [how_many_planet + 1];
typ sma     [how_many_planet + 1];
typ ecc     [how_many_planet + 1];
typ lbd     [how_many_planet + 1];
typ vrp     [how_many_planet + 1];
typ Lbd_0   [how_many_planet + 1];
int p_i     [how_many_planet + 1];
int k_ij    [how_many_planet + 1][how_many_planet + 1];
int q_ij    [how_many_planet]    [how_many_planet + 1];
int k_i     [how_many_planet];
int k_star_i[how_many_planet];
int q_i     [how_many_planet];
struct rational rat_k_ij    [how_many_planet + 1][how_many_planet + 1];
struct rational rat_q_ij    [how_many_planet]    [how_many_planet + 1];
struct rational rat_k_i     [how_many_planet];
struct rational rat_k_star_i[how_many_planet];
struct rational rat_q_i     [how_many_planet];
int * subchain;
int * chain_missed;
typ m0;
int how_many_resonant;
int how_many_missed;



/******** A generic term of the Hamiltonian for any resonance (except 1:1) takes the form                              ********/
/******** C_p,p+q^(k) * (2*D_i/Lbd_i)^m/2 * (2*D_j/Lbd_j)^(n-m)/2 * cos(l*p*lbd_i-l*(p+q)*lbd_j+r*vrp_i+(l*q-r)*vrp_j) ********/
/******** The array qnmlr is filled consequently. The first 10 terms correspond to q = 0 (resonance 0:0), the 15 next  ********/
/******** to q = 1 (resonance p:p+1), the 3 next to q = 2 (resonance p:p+2) and the 4 last to q = 3 (resonance p:p+3)  ********/
/******** For example, the coefficients q,n,m,l,r corresponding to C_p,p+2^(3) are given by qnmlr[10+15+3-1]           ********/
typ qnmlr[32][5] = {
      {0.,0.,0.,0.,0.},
      {0.,2.,2.,0.,0.},
      {0.,2.,0.,0.,0.},
      {0.,2.,1.,0.,1.},
      {0.,4.,4.,0.,0.},
      {0.,4.,0.,0.,0.},
      {0.,4.,3.,0.,1.},
      {0.,4.,2.,0.,0.},
      {0.,4.,2.,0.,2.},
      {0.,4.,1.,0.,1.},
      {1.,1.,1.,1.,1.},
      {1.,1.,0.,1.,0.},
      {1.,2.,2.,2.,2.},
      {1.,2.,1.,2.,1.},
      {1.,2.,0.,2.,0.},
      {1.,3.,3.,1.,1.},
      {1.,3.,2.,1.,0.},
      {1.,3.,2.,1.,2.},
      {1.,3.,1.,1.,1.},
      {1.,3.,1.,1.,-1.},
      {1.,3.,0.,1.,0.},
      {1.,3.,3.,3.,3.},
      {1.,3.,2.,3.,2.},
      {1.,3.,1.,3.,1.},
      {1.,3.,0.,3.,0.},
      {2.,2.,2.,1.,2.},
      {2.,2.,1.,1.,1.},
      {2.,2.,0.,1.,0.},
      {3.,3.,3.,1.,3.},
      {3.,3.,2.,1.,2.},
      {3.,3.,1.,1.,1.},
      {3.,3.,0.,1.,0.}
};


void init(){

      /******** Calls all the initialization-related functions ********/
      
      time_t t;
      time(&t);
      srand((unsigned) t);
      initialization();
      chain_validity();
      array_init();
      array2rational();
      resonance_init();
      if (how_many_resonant > 2){
            l_ij_init();
            NoverD_init();
            rat_c_i_init();
      }
      transformation_init();
      transpose_inv_init();
      verification();
      matrix_fill();
      transformation_display();
      Cppq_init();
      Hamiltonian_display();
      X_old_init(X_old_t0);
      old2new(X_old_t0, X_new_t0, X_uv_t0);
}


void initialization(){

      /******** Initializes the code ********/
      
      int i, j, k;
      
      /******** Initializing the planetary masses and the resonance chain ********/
      how_many_resonant = 0;
      typ buffer_mass  [how_many_planet + 1] = body_masses;
      typ buffer_sma   [how_many_planet]     = body_sma;
      typ buffer_ecc   [how_many_planet]     = body_ecc;
      typ buffer_lbd   [how_many_planet]     = body_lambda;
      typ buffer_vrp   [how_many_planet]     = body_varpi;
      typ buffer_chain [how_many_planet]     = resonance_chain;
      for (i = 1; i < how_many_planet + 1; i ++){
            p_i[i] = buffer_chain[i - 1];
            sma[i] = buffer_sma  [i - 1];
            ecc[i] = buffer_ecc  [i - 1];
            lbd[i] = buffer_lbd  [i - 1];
            vrp[i] = buffer_vrp  [i - 1];
      }
      for (i = 0; i < how_many_planet + 1; i ++){
            masses[i] = buffer_mass[i];
      }
      m0 = *masses;
      for (i = 1; i < how_many_planet + 1; i ++){
            Lbd_0[i] = masses[i]*sqrt(G*m0*sma[i]);
      }
      
      for (i = 1; i < how_many_planet + 1; i ++){
            if (p_i[i] != 0 && (i == how_many_planet || p_i[i + 1] != p_i[i])){ //If the considered planet is resonant and is not the first member of a co-orbital pair
                  how_many_resonant ++;
            }
      }
      how_many_missed = how_many_planet - how_many_resonant;
      subchain        = (int *)malloc((how_many_resonant + 1)*sizeof(int));
      chain_missed    = (int *)malloc((how_many_missed   + 1)*sizeof(int));
      j               = 1;
      k               = 1;
      for (i = 1; i < how_many_planet + 1; i ++){
            if (p_i[i] != 0 && (i == how_many_planet || p_i[i + 1] != p_i[i])){ //If the considered planet is resonant and is not the first member of a co-orbital pair
                  *(subchain + j) = i;
                  j ++;
            }
            else {
                  *(chain_missed + k) = i;
                  k ++;
            }
      }
}


void chain_validity(){

      /******** This function makes sure that the chain given by the user is valid ********/

      int i;

      /******** The system cannot have less than two planets *********/
      if (how_many_planet < 2){
            fprintf(stderr, "\nError: The system must have at least two planets.\n");
            abort();
      }
      
      /******** The chain cannot have only one resonant planet or negative nominal periods ********/
      int how_many_non_secular = 0;
      for (i = 1; i <= how_many_planet; i ++){
            if (p_i[i] < 0){
                  fprintf(stderr, "\nError: An invalid resonance chain was given. The nominal orbital periods cannot be negative.\n");
                  abort();
            }
            if (p_i[i] != 0){
                  how_many_non_secular ++;
            }
      }
      if (how_many_non_secular == 1){
            fprintf(stderr, "\nError: An invalid resonance chain was given. The chain cannot contain exactly one resonant planet.\n");
            abort();
      }
      
      /******** The nominal periods must increase ********/
      int largest_orbital_period_so_far = 0;
      for (i = 1; i < how_many_planet + 1; i ++){
            if (p_i[i] != 0 && p_i[i] < largest_orbital_period_so_far){
                  fprintf(stderr,
                  "\nError: An invalid resonance chain was given. The nominal orbital periods must be given in increasing order (planets are ordered from innermost to outermost).\n"
                  );
                  abort();
            }
            if (p_i[i] != 0){
                  largest_orbital_period_so_far = p_i[i];
            }
      }
      
      /******** There can't be secular planets between co-orbital planets ********/
      int latest_orbital_period = 0;
      int secular_in_between    = 0;
      for (i = 1; i < how_many_planet + 1; i ++){
            if (p_i[i] != 0 && p_i[i] == latest_orbital_period && secular_in_between){
                  fprintf(stderr,
                  "\nError: An invalid resonance chain was given. A non-resonant planet cannot be between co-orbital planets (planets are ordered from innermost to outermost).\n");
                  abort();
            }
            if (p_i[i] == 0){
                  secular_in_between = 1;
            }
            else {
                  secular_in_between = 0;
                  latest_orbital_period = p_i[i];
            }
      }
}


void array_init(){

      /******** Initializes the arrays k_ij, q_ij, k_i, k_star_i and q_i ********/
      
      int i, j;
      
      /******** k_ij ********/
      for (i = 1; i < how_many_resonant + 1; i ++){
            for (j = 1; j < how_many_resonant + 1; j ++){
                  if (i != j){
                        k_ij[i][j] = p_i[subchain[i]] / gcd(p_i[subchain[i]],p_i[subchain[j]]);
                  }
            }
      }
      
      /******** q_ij ********/
      for (i = 1; i < how_many_resonant; i ++){
            for (j = i + 1; j < how_many_resonant + 1; j ++){
                  q_ij[i][j] = k_ij[j][i] - k_ij[i][j];
            }
      }
      
      /******** k_i, k_star_i and q_i ********/
      for (i = 1; i < how_many_resonant; i ++){
            k_i     [i] = k_ij[i][i + 1];
            k_star_i[i] = k_ij[i + 1][i];
            q_i     [i] = q_ij[i][i + 1];
      }

}


void array2rational(){

      /******** Converts the arrays k_ij, q_ij, k_i, k_star_i and q_i from integer to rational ********/
      /******** Updates  the arrays rat_ accordingly                                           ********/

      int i, j;
      
      /******** k_ij ********/
      for (i = 1; i < how_many_resonant + 1; i ++){
            for (j = 1; j < how_many_resonant + 1; j ++){
                  if (i != j){
                        rat_k_ij[i][j] = int2rat(k_ij[i][j]);
                  }
            }
      }
      
      /******** q_ij ********/
      for (i = 1; i < how_many_resonant; i ++){
            for (j = i + 1; j < how_many_resonant + 1; j ++){
                  rat_q_ij[i][j] = int2rat(q_ij[i][j]);
            }
      }
      
      /******** k_i, k_star_i and q_i ********/
      for (i = 1; i < how_many_resonant; i ++){
            rat_k_i     [i] = int2rat(k_i[i]);
            rat_k_star_i[i] = int2rat(k_star_i[i]);
            rat_q_i     [i] = int2rat(q_i[i]);
      }

}


void deallocation(){

      /******** Deallocates the global memory ********/

      free(subchain);
      subchain = NULL;
      free(chain_missed);
      chain_missed = NULL;

}


struct pairOfReal b_12(typ alpha, typ precision){

      /******** Returns the Laplace coefficients b_1/2^0(alpha) and b_1/2^1(alpha) ********/
      /******** The relative error on b_1/2^0(alpha) is better than precision      ********/
      /******** Based on Eqs. (59) & (60) p. 500 of Brouwer & Clemence             ********/
      
      struct pairOfReal toBeReturned;
      typ m_buffer;
      typ m                 = 1.;
      typ n                 = sqrt(1. - alpha*alpha);
      typ current_precision = fabs(m - n)/m;
      typ b_12_1            = alpha;
      typ loop              = 2.;
      
      while(current_precision > precision){
            m_buffer          = m;
            m                 = 0.5*(m + n);
            n                 = sqrt(m_buffer*n);
            current_precision = fabs(m - n)/m;
            b_12_1           += loop*(m*m - n*n)/alpha;
            loop             *= 2.;
      }
      
      b_12_1          /= m;
      toBeReturned.fst = 2./m;
      toBeReturned.snd = b_12_1;
      return toBeReturned;
}


struct pairOfReal b_k2(typ alpha, struct pairOfReal b_km12, int k){

      /******** Returns the Laplace coefficients b_k/2^0(alpha) and b_k/2^1(alpha) ********/
      /******** from the knowlegde of b_(k-1)/2^0(alpha) and b_(k-1)/2^1(alpha),   ********/
      /******** using Eqs. (63) & (64) p. 501 of Brouwer & Clemence                ********/
      
      typ kk = (typ) k;
      typ coef1, coef2, coef3, coef4;
      if      (k == 3){
            coef1 = 0.5;  coef2 = -1.0;  coef3 = -0.5; coef4 = 1.0;
      }
      else if (k == 5){
            coef1 = 1.5;  coef2 = 1.0;   coef3 = 0.5;  coef4 = 3.0;
      }
      else if (k == 7){
            coef1 = 2.5;  coef2 = 3.0;   coef3 = 1.5;  coef4 = 5.0;
      }
      else{
            fprintf(stderr, "\nError: The argument k must be 3, 5 or 7 in function b_k2.\n");
            abort();
      }
      
      struct pairOfReal toBeReturned;
      typ b_km12_0     = b_km12.fst;
      typ b_km12_1     = b_km12.snd;
      typ malp2        = 1.0 - alpha*alpha;
      typ palp2        = 1.0 + alpha*alpha;
      typ denominator  = coef1*malp2*malp2;
      typ b_k2_0       = (coef1*palp2*b_km12_0 + coef2*alpha*b_km12_1)/denominator;
      typ b_k2_1       = (coef3*palp2*b_km12_1 + coef4*alpha*b_km12_0)/denominator;
      toBeReturned.fst = b_k2_0;
      toBeReturned.snd = b_k2_1;
      return toBeReturned;
}


int gcd(int a, int b){

      /******** Computes the greatest common divisor of a and b with Euclide algorithm ********/
      /******** a and b are both positive                                              ********/

      
      if (a < b){
            return gcd(b, a);
      }
      
      if (b == 0){
            if (a == 0){
                  fprintf(stderr, "\nError: Trying to compute gcd(0,0) in function gcd.\n");
                  abort();
            }
            return a;
      }
      
      int remainder;
      
      while (a > 1 && b > 1){
            remainder = a % b;
            if (remainder == 0){
                  return b;
            }
            a = b;
            b = remainder;
      }
      return 1;
}


int lcm(int a, int b){

      /******** Computes the least common multiple of a and b ********/
      /******** a and b are both positive                     ********/
      
      return a*b/gcd(a, b);
}


struct rational ratmul(struct rational r1, struct rational r2){

      /******** Computes the product of the two rational numbers r1 and r2 ********/
      
      struct rational r12;
      int p1 = r1.numerator;
      int p2 = r2.numerator;
      int q1 = r1.denominator;
      int q2 = r2.denominator;
      int r12_num, r12_denom, irrec;
      
      r12_num         = p1*p2;
      r12_denom       = q1*q2;
      irrec           = gcd(abs(r12_num),r12_denom);
      r12_num         = r12_num  /  irrec;
      r12_denom       = r12_denom / irrec;
      r12.numerator   = r12_num;
      r12.denominator = r12_denom;
      
      return r12;
}


struct rational ratadd(struct rational r1, struct rational r2){

      /******** Computes the sum of the two rational numbers r1 and r2 ********/
      
      struct rational r12;
      int p1 = r1.numerator;
      int p2 = r2.numerator;
      int q1 = r1.denominator;
      int q2 = r2.denominator;
      int r12_num, r12_denom, irrec;
      
      r12_num         = p1*q2 + p2*q1;
      r12_denom       = q1*q2;
      irrec           = gcd(abs(r12_num), r12_denom);
      r12_num         = r12_num  /  irrec;
      r12_denom       = r12_denom / irrec;
      r12.numerator   = r12_num;
      r12.denominator = r12_denom;
      
      return r12;
}


struct rational ratopp(struct rational r1){

      /******** Returns -r1 where r1 is a rational number ********/
      
      struct rational mr1;
      mr1.numerator   = -r1.numerator;
      mr1.denominator =  r1.denominator;
      return mr1;
}


struct rational ratinv(struct rational r1){

      /******** Returns 1/r1 where r1 is a rational number ********/

      struct rational mr1;
      if (r1.numerator == 0){
            fprintf(stderr, "\nError: Trying to divide by 0 in function ratinv.\n");
            abort();
      }
      if (r1.numerator > 0){
            mr1.numerator   = r1.denominator;
            mr1.denominator = r1.numerator;
      }
      else {
            mr1.numerator   = -r1.denominator;
            mr1.denominator = -r1.numerator;
      }
      
      return mr1;
}


struct rational ratdiff(struct rational r1, struct rational r2){

      /******** Computes the difference r1 - r2 ********/
      
      return ratadd(r1, ratopp(r2));
}


struct rational ratover(struct rational r1, struct rational r2){

      /******** Computes the ratio r1/r2 ********/
      
      return ratmul(r1, ratinv(r2));
}


struct rational int2rat(int a){

      /******** Converts the integer a into a rational number ********/
      
      struct rational r;
      r.numerator   = a;
      r.denominator = 1;

      return r;
}


typ rat2real(struct rational r1){

      /******** Converts the rational r1 to floating point number ********/

      typ r;
      int num   = r1.numerator;
      int denom = r1.denominator;
      r = (typ) (((typ) num)/((typ) denom));
      return r;
}


struct rational real2rat(typ r, int targetDenom){

      /******** Converts the real number r to a rational approximation whose denominator is close to targetDenom in order of magnitude ********/

      typ a, a0, R;
      struct rational output;
      
      if (r < 0.){
            return ratopp(real2rat(-r, targetDenom));
      }
      if (targetDenom <= 1){
            return int2rat((int) floor(r));
      }
      a0     = floor(r);
      R      = 1./(r - a0);
      a      = floor(R);
      output = ratadd(int2rat((int) a0), ratinv(int2rat((int) a)));
      if (output.denominator >= targetDenom/2){
            return output;
      }
      else{
            if (r >= rat2real(output)){
                  return ratadd(output, real2rat(r - rat2real(output), targetDenom));
            }
            else{
                  return ratdiff(output, real2rat(rat2real(output) - r, targetDenom));
            }
      }
}


struct rational ratabs(struct rational r1){

      /******** Returns the absolute value of the rational r1 ********/

      struct rational r;
      r.numerator   = abs(r1.numerator);
      r.denominator = r1.denominator;

      return r;
}


void ratprint(struct rational r){

      /******** Prints the rational r ********/
      int num   = r.numerator;
      int denom = r.denominator;
      if (denom == 1){
            printf("%d", num);
      }
      else{
            printf("%d/%d", num, denom);
      }   
}


typ rdm(typ min, typ max){
      
      /******** Returns a random number between min and max according to a uniform distribution ********/

      typ MyRand = ((typ) rand())/((typ) RAND_MAX); //between 0 and 1
      MyRand = min + (max - min)*MyRand;            //between min and max

      return MyRand;
}


typ continuousAngle(typ newAngle, typ oldAngle){

      /******** Maintain continuity of the angle ********/

      typ dAngle = newAngle - oldAngle;
      typ qtient = round(dAngle/(2.*M_PI));
      return newAngle - qtient*2.*M_PI;
}
