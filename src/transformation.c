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


/******** Defining some global variables ********/
struct rational transformation[2*how_many_planet + 1][2*how_many_planet + 1];
struct rational transpose_inv [2*how_many_planet + 1][2*how_many_planet + 1];
typ             Transformation[2*how_many_planet + 1][2*how_many_planet + 1];
typ             Transpose_inv [2*how_many_planet + 1][2*how_many_planet + 1];
struct rational rat_l_ij      [  how_many_planet - 1]  [how_many_planet + 1];
struct rational NoverD        [  how_many_planet - 1]  [how_many_planet - 1][how_many_planet + 1];
struct rational rat_c_i       [  how_many_planet - 1];
int dof                       [  how_many_planet + 1];
int how_many_dof;
int nondof             [  how_many_planet + 1];
int how_many_nondof;


void l_ij_init(){

      /******** Computes the coefficients l_ij and stores them into the array rat_l_ij ********/

      int i,j;
      
      for (i = 1; i < how_many_resonant - 1; i ++){
            for (j = how_many_resonant; j > 0; j --){
                  if (j > i){
                        rat_l_ij[i][j] = int2rat(0);
                  }
                  else if (j == i){
                        rat_l_ij[i][j] = ratover(rat_q_i[i],rat_k_i[i]);
                  }
                  else{
                        rat_l_ij[i][j] = ratadd(rat_l_ij[i][j + 1], ratmul(ratover(rat_q_i[j],rat_k_i[j]),ratover(rat_k_ij[i + 1][j + 1],rat_k_ij[j + 1][i + 1])));
                  }
            }
      }
}


void NoverD_init(){

      /******** Computes the coefficients N_r,s^(i) and D_r,s^(i) under the form of a rational number ********/
      
      int i, r, s;
      
      for (i = 1; i < how_many_resonant - 1; i ++){
            for (r = 1; r <= i; r ++){
                  for (s = r + 1; s < how_many_resonant + 1; s ++){
                        NoverD[i][r][s] = ratabs(ratdiff(ratmul(rat_k_ij[s][r], rat_l_ij[i][s]), ratmul(rat_k_ij[r][s], rat_l_ij[i][r])));
                  }
            }
      }
}


void rat_c_i_init(){

      /******** Computes the coefficients c_i and stores them into rat_c_i ********/

      int i, j, r, s;
      int * N;
      int * D;
      int k; //Number of coefficients N_r,s^(i)/D_r,s^(i) for a fixed i
      struct rational R, c;
      int num, denom, irrec;
      
      for (i = 1; i < how_many_resonant - 1; i ++){
            k = how_many_resonant*i - (i*(i + 1))/2;
            N = (int *)malloc(k * sizeof(int));
            D = (int *)malloc(k * sizeof(int));
            j = 0;
            for (r = 1; r <= i; r ++){
                  for (s = r + 1; s < how_many_resonant + 1; s ++){
                        R     = NoverD[i][r][s];
                        num   = R.numerator;
                        denom = R.denominator;
                        if (q_ij[r][s] <= max_deg && k_ij[s][r] <= max_res){ //If the pair (r,s) is really resonant despite the approximations
                              *(N + j) = num;
                              *(D + j) = denom;
                              j ++;
                        }
                  }
            }
            if (j == 0){
                  num   = 1;
                  denom = 1;
            }
            else{
                  num   = LCM(D, j);
                  denom = GCD(N, j);
                  irrec = gcd(num, denom);
                  num   = num / irrec;
                  denom = denom/irrec;
            }
            c.numerator   = num;
            c.denominator = denom;
            
            rat_c_i[i] = c;
            
            free(N);  N = NULL;
            free(D);  D = NULL;
      }
}


void transformation_init(){

      /******** Fills the matrix M (array transformation) of the transformation ********/

      int i, j;
      int is_coorbital;

      for (i = 1; i < 2*how_many_planet + 1; i ++){
            for (j = 1; j < 2*how_many_planet + 1; j ++){
                  transformation[i][j] = int2rat(0);
            }
      }

      if (how_many_resonant == 0){ //All planets are secular
            for (i = 1; i < 2*how_many_planet + 1; i ++){
                  transformation[i][i] = int2rat(1); //The transformation is identity
            }
      }
      else if (how_many_resonant == 1){ //All planets are secular except for a pair or a triplet of co-orbitals
            for (i = 1; i < how_many_planet + 1; i ++){
                  transformation[i][i] = int2rat(1);
                  is_coorbital         = p_i[i];
                  if (is_coorbital && i != subchain[how_many_resonant]){
                        transformation[i][i + 1] = int2rat(-1);
                  }
            }
            for (i = how_many_planet + 1; i < 2*how_many_planet + 1; i ++){
                  transformation[i][i] = int2rat(1);
            }
      }
      else{ //There are at least two resonant non-coorbital planets
      
            /******** Upper-left block of the transformation matrix ********/
            /******** Resonant planets and trailing co-orbitals ********/
            for (i = 1; i < how_many_resonant - 1; i ++){
                  transformation[subchain[i]][subchain[i]]     = ratmul(ratinv(rat_c_i[i]), ratover(rat_k_i[i], rat_q_i[i]));
                  transformation[subchain[i]][subchain[i + 1]] = ratopp(ratmul(ratinv(rat_c_i[i]),ratadd(ratover(rat_k_star_i[i],rat_q_i[i]),ratover(rat_k_i[i + 1],rat_q_i[i + 1]))));
                  transformation[subchain[i]][subchain[i + 2]] = ratmul(ratinv(rat_c_i[i]),ratover(rat_k_star_i[i + 1],rat_q_i[i + 1]));
            }
            transformation[subchain[how_many_resonant - 1]][subchain[how_many_resonant - 1]] = int2rat(1);
            transformation[subchain[how_many_resonant - 1]][subchain[how_many_resonant]]     = int2rat(-1);
            transformation[subchain[how_many_resonant]]    [subchain[how_many_resonant]]     = ratover(  rat_k_star_i[how_many_resonant - 1],rat_q_i[how_many_resonant - 1]);
            transformation[subchain[how_many_resonant]]    [subchain[how_many_resonant - 1]] = ratopp(ratover(rat_k_i[how_many_resonant - 1],rat_q_i[how_many_resonant - 1]));
            /******** Secular planets and leading co-orbitals ********/
            for (i = 1; i < how_many_missed + 1; i ++){
                  transformation[chain_missed[i]][chain_missed[i]] = int2rat(1);
                  is_coorbital                                     = p_i[chain_missed[i]];
                  if (is_coorbital){
                        transformation[chain_missed[i]][chain_missed[i] + 1] = int2rat(-1);
                  }
            }
            
            /******** Bottom-right block of the transformation matrix ********/
            for (i = how_many_planet + 1; i < 2*how_many_planet + 1; i ++){
                  transformation[i][i] = int2rat(1); //The Bottom-right block is identity
            }
            
            /******** Bottom-left block of the transformation matrix ********/
            for (i = how_many_planet + 1; i < 2*how_many_planet + 1; i ++){
                  transformation[i][subchain[how_many_resonant - 1]] = ratopp(ratover(rat_k_i[how_many_resonant - 1],rat_q_i[how_many_resonant - 1]));
                  transformation[i][subchain[how_many_resonant]]     = ratover(  rat_k_star_i[how_many_resonant - 1],rat_q_i[how_many_resonant - 1]);
            }
      }
}


void transpose_inv_init(){

      /******** Fills the matrix tM^-1 (array transpose_inv) of the inverse transpose of the transformation ********/

      int i,j;
      int is_coorbital;
      struct rational Gamma_factor;

      for (i = 1; i < 2*how_many_planet + 1; i ++){
            for (j = 1; j < 2*how_many_planet + 1; j ++){
                  transpose_inv[i][j] = int2rat(0);
            }
      }
      
      if (how_many_resonant == 0){ //All planets are secular
            for (i = 1; i < 2*how_many_planet + 1; i ++){
                  transpose_inv[i][i] = int2rat(1); //The transformation is identity
            }
      }
      else if (how_many_resonant == 1){ //All planets are secular except for a pair or a triplet of co-orbitals
            for (i = 1; i < how_many_planet + 1; i ++){
                  transpose_inv[i][i] = int2rat(1);
                  is_coorbital        = p_i[i];
                  if (is_coorbital){
                        for (j = 1; j < i; j ++){
                              if (p_i[j]){
                                    transpose_inv[i][j] = int2rat(1);
                              }
                        }
                  }
            }
            for (i = how_many_planet + 1; i < 2*how_many_planet + 1; i ++){
                  transpose_inv[i][i] = int2rat(1);
            }
      }
      else{ //There are at least two resonant non-coorbital planets
      
            /******** Upper-left block of the inverse transpose matrix ********/
            /******** Resonant planets and trailing co-orbitals ********/
            Gamma_factor = ratover(rat_k_i[how_many_resonant - 1],rat_q_i[how_many_resonant - 1]);
            for (i = 1; i < how_many_resonant - 1; i ++){
                  for (j = 1; j <= i; j ++){
                        transpose_inv[subchain[i]][subchain[j]] = ratmul(rat_c_i[i], rat_l_ij[i][j]);
                  }
            }
            for (j = 1; j < how_many_resonant; j ++){
                  transpose_inv[subchain[how_many_resonant - 1]][subchain[j]] = ratmul(Gamma_factor,ratover(rat_k_ij[how_many_resonant][j],rat_k_ij[j][how_many_resonant]));
            }
            transpose_inv[subchain[how_many_resonant - 1]][subchain[how_many_resonant]] = Gamma_factor;
            for (j = 1; j < how_many_resonant + 1; j ++){
                  transpose_inv[subchain[how_many_resonant]][subchain[j]] = int2rat(1);
            }
            
            /******** Secular planets and leading co-orbitals ********/
            for (i = how_many_missed; i > 0; i --){
                  transpose_inv[chain_missed[i]][chain_missed[i]] = int2rat(1);
                  is_coorbital = p_i[chain_missed[i]];
                  if (is_coorbital){
                        for (j = chain_missed[i] + 1; j < how_many_planet + 1; j ++){
                              transpose_inv[j][chain_missed[i]] = transpose_inv[j][chain_missed[i] + 1];
                        }
                        transpose_inv[subchain[how_many_resonant - 1]][chain_missed[i]] = transpose_inv[subchain[how_many_resonant - 1]][chain_missed[i] + 1];
                  }
            }
            
            /******** Upper and bottom-right block of the inverse transpose matrix ********/
            for (j = how_many_planet + 1; j < 2*how_many_planet + 1; j ++){
                  transpose_inv[subchain[how_many_resonant]][j] = int2rat(-1);
                  transpose_inv[j][j]                           = int2rat(1);
            }
      }
}


void verification(){

      /******** Checks that transpose_inv is indeed t(transformation)^-1 by computing ********/
      /******** the product transformation*t(transpose_inv), that should be identity  ********/

      int i, j, k;
      struct rational coef;
      int num, denom;
      
      for (i = 1; i < 2*how_many_planet + 1; i ++){
            for (j = 1; j < 2*how_many_planet + 1; j ++){
                  coef = int2rat(0);
                  for (k = 1; k < 2*how_many_planet + 1; k ++){
                        coef = ratadd(coef,ratmul(transformation[i][k],transpose_inv[j][k]));
                  }
                  num   = coef.numerator;
                  denom = coef.denominator;
                  if (i != j && (num != 0 || denom != 1) || i == j && (num != 1 || denom != 1)){
                        fprintf(stderr, "\nError: The transformation failed verification. This is a bug.\n");
                        abort();
                  }
            }
      }
}


void matrix_fill(){

      /******** Fills the arrays Transformation and Transpose_inv ********/

      int i, j;
      for (i = 1; i < 2*how_many_planet + 1; i ++){
            for (j = 1; j < 2*how_many_planet + 1; j ++){
                  Transformation[i][j] = rat2real(transformation[i][j]);
                  Transpose_inv [i][j] = rat2real(transpose_inv [i][j]);
            }
      }
}


void transformation_display(){

      /******** Displays in the terminal the change of variables performed by Aptidal ********/
      
      int i, j, how_many_secular;
      struct rational coef;
      int num, denom;
      int first_print;
      
      how_many_secular = 0;
      for (i = 1; i < how_many_planet + 1; i ++){
            if (p_i[i] == 0){
                  how_many_secular ++;
            }
      }
      
      printf("\n-------------------------------------------------------------------------------------------\n\n");
      printf("The resonance chain to be studied is ");
      for (i = 1; i < how_many_planet; i ++){
            printf("%d:", p_i[i]);
      }
      printf("%d\n", p_i[how_many_planet]);
      printf("The system currently has %d degrees of freedom.\n\n", 2*how_many_planet);
      printf("Aptidal will perform a canonical transformation (lbd, -vrp; Lbd, D) -> (phi, sig; Phi, D) adapted to\n");
      printf("the resonance where (lbd, -vrp; Lbd, D) are the traditional Poincaré variables.\n");
      printf("lbd (lambda) is the mean longitude, vrp (varpi) is the longitude of the ascending node, Lambda (resp. Phi)\n");
      printf("is the action conjugated to lambda (resp. phi) and D (the AMD) is conjugated to both sigma and -varpi.\n\n");
      printf("The transformation reads\n\n");
      
      /******** Displaying (phi_j, sig_j) as a function of (lbd_j,-vrp_j) ********/
      for (i = 1; i < how_many_planet + 1; i ++){
            printf(" phi_%d =", i);
            first_print = 1;
            for (j = 1; j < how_many_planet + 1; j ++){
                  coef  = transformation[i][j];
                  num   = coef.numerator;
                  denom = coef.denominator;
                  if (num != 0){
                        if (num < 0){
                              printf(" - ");
                        }
                        else if (num > 0 && !first_print){
                              printf(" + ");
                        }
                        else{
                              printf(" ");
                        }
                        if (abs(num) == 1 && denom == 1){
                              printf("lbd_%d", j);
                        }
                        else{
                              ratprint(ratabs(coef));
                              printf(" lbd_%d", j);
                        }
                        first_print = 0;
                  }
            }
            printf("\n");
      }
      printf(" sig_j =");
      first_print = 1;
      for (j = 1; j < how_many_planet + 1; j ++){
            coef  = transformation[i][j];
            num   = coef.numerator;
            denom = coef.denominator;
            if (num != 0){
                  if (num < 0){
                        printf(" - ");
                  }
                  else if (num > 0 && !first_print){
                        printf(" + ");
                  }
                  else{
                        printf(" ");
                  }
                  if (abs(num) == 1 && denom == 1){
                        printf("lbd_%d", j);
                  }
                  else{
                        ratprint(ratabs(coef));
                        printf(" lbd_%d", j);
                  }
                  first_print = 0;
            }
      }
      printf(" - vrp_j,  1 <= j <= %d\n\n", how_many_planet);
      
      /******** Displaying (Phi_j, D_j) as a function of (Lbd_j, D_j) ********/
      printf("It is canonical if the actions are transformed according to (the D_j are unchanged)\n\n");
      for (i = 1; i < how_many_planet + 1; i ++){
            printf(" Phi_%d =", i);
            first_print = 1;
            for (j = 1; j < how_many_planet + 1; j ++){
                  coef  = transpose_inv[i][j];
                  num   = coef.numerator;
                  denom = coef.denominator;
                  if (num != 0){
                        if (num < 0){
                              printf(" - ");
                        }
                        else if (num > 0 && !first_print){
                              printf(" + ");
                        }
                        else{
                              printf(" ");
                        }
                        if (abs(num) == 1 && denom == 1){
                              printf("Lbd_%d", j);
                        }
                        else{
                              ratprint(ratabs(coef));
                              printf(" Lbd_%d", j);
                        }
                        first_print = 0;
                  }
            }
            for (j = how_many_planet + 1; j < 2*how_many_planet + 1; j ++){
                  coef  = transpose_inv[i][j];
                  num   = coef.numerator;
                  denom = coef.denominator;
                  if (num != 0){
                        if (num < 0){
                              printf(" - ");
                        }
                        else if (num > 0 && !first_print){
                              printf(" + ");
                        }
                        else{
                              printf(" ");
                        }
                        if (abs(num) == 1 && denom == 1){
                              printf("D_%d", j - how_many_planet);
                        }
                        else{
                              ratprint(ratabs(coef));
                              printf(" D_%d", j - how_many_planet);
                        }
                        first_print = 0;
                  }
            }
            printf("\n");
      }
      printf("\n");
      
      /******** Displaying (lbd_j, -vrp_j) as a function of (phi_j, sig_j) ********/
      printf("The inverse transformation is\n\n");
      for (i = 1; i < how_many_planet + 1; i ++){
            printf(" lbd_%d =", i);
            first_print = 1;
            for (j = 1; j < how_many_planet + 1; j ++){
                  coef  = transpose_inv[j][i];
                  num   = coef.numerator;
                  denom = coef.denominator;
                  if (num != 0){
                        if (num < 0){
                              printf(" - ");
                        }
                        else if (num > 0 && !first_print){
                              printf(" + ");
                        }
                        else{
                              printf(" ");
                        }
                        if (abs(num) == 1 && denom == 1){
                              printf("phi_%d", j);
                        }
                        else{
                              ratprint(ratabs(coef));
                              printf(" phi_%d", j);
                        }
                        first_print = 0;
                  }
            }
            printf("\n");
      }
      if (how_many_resonant < 2){
            printf(" vrp_j = - sig_j,  1 <= j <= %d\n\n", how_many_planet);
      }
      else{
            printf(" vrp_j = phi_%d - sig_j,  1 <= j <= %d\n\n", subchain[how_many_resonant], how_many_planet);
      }
      
      /******** Displaying (Lbd_j, D_j) as a function of (Phi_j, D_j) ********/
      for (i = 1; i < how_many_planet + 1; i ++){
            printf(" Lbd_%d =", i);
            first_print = 1;
            for (j = 1; j < how_many_planet + 1; j ++){
                  coef  = transformation[j][i];
                  num   = coef.numerator;
                  denom = coef.denominator;
                  if (num != 0){
                        if (num < 0){
                              printf(" - ");
                        }
                        else if (num > 0 && !first_print){
                              printf(" + ");
                        }
                        else{
                              printf(" ");
                        }
                        if (abs(num) == 1 && denom == 1){
                              printf("Phi_%d", j);
                        }
                        else{
                              ratprint(ratabs(coef));
                              printf(" Phi_%d", j);
                        }
                        first_print = 0;
                  }
            }
            for (j = how_many_planet + 1; j < 2*how_many_planet + 1; j ++){
                  coef  = transformation[j][i];
                  num   = coef.numerator;
                  denom = coef.denominator;
                  if (num != 0){
                        if (num < 0){
                              printf(" - ");
                        }
                        else if (num > 0 && !first_print){
                              printf(" + ");
                        }
                        else{
                              printf(" ");
                        }
                        if (abs(num) == 1 && denom == 1){
                              printf("D_%d", j - how_many_planet);
                        }
                        else{
                              ratprint(ratabs(coef));
                              printf(" D_%d", j - how_many_planet);
                        }
                        first_print = 0;
                  }
            }
            printf("\n");
      }
      printf("\n");
      
      /******** Displaying the degrees of freedom ********/
      if (how_many_resonant == 0){
            printf("After averaging over the fast angles, the %d remaining degrees of freedom are\n\n", how_many_planet);
      }
      else if (how_many_resonant == 1){
            if (how_many_secular == 0){
                  printf("After averaging over the fast angle, the %d remaining degrees of freedom are\n\n",  2*how_many_planet - 1);
            }
            else{
                  printf("After averaging over the fast angles, the %d remaining degrees of freedom are\n\n", 2*how_many_planet - how_many_secular - 1);
            }
      }
      else{
            if (how_many_secular == 0){
                  printf("After averaging over the fast angle, the %d remaining degrees of freedom are\n\n",  2*how_many_planet - 2);
            }
            else{
                  printf("After averaging over the fast angles, the %d remaining degrees of freedom are\n\n", 2*how_many_planet - how_many_secular - 2);
            }
      }
      
      /******** Printing the degrees of freedom and storing their indexes in array dof ********/
      how_many_dof    = 0;
      how_many_nondof = 0;
      printf(" ");
      for (i = 1; i < how_many_planet + 1; i ++){
            if (how_many_resonant >= 2){
                  if (i != subchain[how_many_resonant - 1] && i != subchain[how_many_resonant] && p_i[i] != 0){
                        printf("(phi_%d; Phi_%d), ", i, i);
                        how_many_dof ++;
                        dof[how_many_dof] = i;
                  }
                  else{
                        how_many_nondof ++;
                        nondof[how_many_nondof] = i;
                  }
            }
            else if (how_many_resonant == 1){
                  if (i != subchain[how_many_resonant] && p_i[i] != 0){
                        printf("(phi_%d; Phi_%d), ", i, i);
                        how_many_dof ++;
                        dof[how_many_dof] = i;
                  }
                  else{
                        how_many_nondof ++;
                        nondof[how_many_nondof] = i;
                  }
            }
            else if (how_many_resonant == 0){
                  how_many_nondof ++;
                  nondof[how_many_nondof] = i;
            }
      }
      for (i = 1; i < how_many_planet + 1; i ++){
            if (i < how_many_planet){
                  printf("(sig_%d; D_%d), ", i, i);
            }
            else{
                  printf("(sig_%d; D_%d)\n\n", i, i);
            }
      }
}


void Hamiltonian_display(){

      /******** Displays the Hamiltonian in the terminal ********/
      
      int i, j, k, pi, pj, the_gcd, how_many_backward;
      int p, q, n, m, l, r, s;
      
      printf("-------------------------------------------------------------------------------------------\n\n");
      printf("The Hamiltonian of the model reads, in the old variables\n\n H = ");
      
      /******** Displaying the Keplerian part ********/
      printf("- sum_{1 <= j <= %d} beta_j**3 * mu_j**2 / (2*Lbd_j**2)\n\n", how_many_planet);
      
      /******** Displaying the perturbation ********/
      for (i = 1; i <= how_many_planet; i ++){
            for (j = i + 1; j <= how_many_planet; j ++){
                  
                  /******** Retrieving the value of p for a resonance p : p + q ********/
                  pi      = p_i[i];
                  pj      = p_i[j];
                  if (pi == 0 || pj == 0){
                        p = 0;
                  }
                  else{
                        the_gcd = gcd(pi, pj);
                        p       = pi/the_gcd;
                        pi      = pi/the_gcd;
                        pj      = pj/the_gcd;
                  }
                  
                  if ((non_resonant_bool && (max_deg + one_more_deg_bool) >= 2) || ((pj - pi) <= max_deg && pj <= max_res)){
                        printf(" + G*m_%d*m_%d/a_%d * (\n", i, j, j);
                        /******** Printing the secular contribution and the non-co-orbital MMR contribution ********/
                        for (k = 1; k < 32; k ++){
                              if (Cppq[i][j][k] != 0.){ //If the corresponding term exists in the Hamiltonian
                                    q = (int) qnmlr[k][0];  n = (int) qnmlr[k][1];  m = (int) qnmlr[k][2];  l = (int) qnmlr[k][3];  r = (int) qnmlr[k][4];
                              
                                    if (Cppq[i][j][k] > 0.){
                                          if (k > 1){
                                                printf("    + %.10lf", Cppq[i][j][k]);
                                          }
                                          else{
                                                printf("    %.10lf", Cppq[i][j][k]);
                                          }
                                    }
                                    else{
                                          printf("    - %.10lf", -Cppq[i][j][k]);
                                    }
                                    how_many_backward = (int) log10(fabs(Cppq[i][j][k]));
                                    for (s = 0; s < how_many_backward; s ++){
                                          printf("\b");
                                    }
                                    for (s = 0; s < how_many_backward; s ++){
                                          printf(" ");
                                    }
                                    for (s = 0; s < how_many_backward; s ++){
                                          printf("\b");
                                    }
                                    printf(" ");
                                    if (m != 0){
                                          if (m == 1){
                                                printf("* sqrt(2*D_%d/Lbd_%d)    ", i, i);
                                          }
                                          else{
                                                printf("* sqrt(2*D_%d/Lbd_%d)**%d ", i, i, m);
                                          }
                                    }
                                    if (n - m != 0){
                                          if (n - m == 1){
                                                printf("* sqrt(2*D_%d/Lbd_%d)    ", j, j);
                                          }
                                          else{
                                                printf("* sqrt(2*D_%d/Lbd_%d)**%d ", j, j, n - m);
                                          }
                                    }
                                    if (m == 0 && max_deg >= 2){
                                          printf("                       ");
                                    }
                                    if (n - m == 0 && max_deg >= 2){
                                          printf("                       ");
                                    }
                                    if (l*p != 0 || l*(p + q) != 0 || r != 0 || l*q - r != 0){ //If there is a cosine factor
                                          printf("* cos(");
                                          if (l*p != 0){
                                                if (l*p == 1){
                                                      printf("lbd_%d", i);
                                                }
                                                else{
                                                      printf("%d lbd_%d", l*p, i);
                                                }
                                          }
                                          if (l*(p + q) != 0){
                                                if (l*(p + q) == 1){
                                                      printf(" - lbd_%d", j);
                                                }
                                                else{
                                                      printf(" - %d lbd_%d", l*(p + q), j);
                                                }
                                          }
                                          if (r != 0){
                                                if (r < -1){
                                                      printf(" - %d vrp_%d", -r, i);
                                                }
                                                else if (r == -1){
                                                      printf(" - vrp_%d", i);
                                                }
                                                else if (r == 1){
                                                      if (l*p != 0 || l*(p + q) != 0){
                                                            printf(" + vrp_%d", i);
                                                      }
                                                      else{
                                                            printf("vrp_%d", i);
                                                      }
                                                }
                                                else{
                                                      if (l*p != 0 || l*(p + q) != 0){
                                                            printf(" + %d vrp_%d", r, i);
                                                      }
                                                      else{
                                                            printf("%d vrp_%d", r, i);
                                                      }
                                                }
                                          }
                                          if (l*q - r != 0){
                                                if (l*q - r < -1){
                                                      printf(" - %d vrp_%d", r - l*q, j);
                                                }
                                                else if (l*q - r == -1){
                                                      printf(" - vrp_%d", j);
                                                }
                                                else if (l*q - r == 1){
                                                      printf(" + vrp_%d", j);
                                                }
                                                else{
                                                      printf(" + %d vrp_%d", l*q - r, j);
                                                }
                                          }
                                          printf(")");
                                    }
                                    printf("\n");
                              }
                        }
                  
                        /******** Printing the co-orbital contribution ********/
                        if (pi != 0 && pj != 0 && pi == pj){ //The pair (i,j) is in a MMR
                              printf("      cos(lbd_%d - lbd_%d) - 1/sqrt(2 - 2 cos(lbd_%d - lbd_%d))\n", i, j, i, j);
                              if (max_deg >= 2){
                                    printf("    - (D_%d/Lbd_%d + D_%d/Lbd_%d) * cos(lbd_%d - lbd_%d)\n", i, i, j, j, i, j);
                                    printf("    + (D_%d/Lbd_%d + D_%d/Lbd_%d) * [5 cos(2 lbd_%d - 2 lbd_%d) + 8 cos(lbd_%d - lbd_%d) - 13] / [4 sqrt(2 - 2 cos(lbd_%d - lbd_%d))**5]\n",
                                    i, i, j, j, i, j, i, j, i, j);
                                    printf("    + sqrt(2*D_%d/Lbd_%d) * sqrt(2*D_%d/Lbd_%d) * cos(vrp_%d - vrp_%d + 2 lbd_%d - 2 lbd_%d)\n", i, i, j, j, j, i, i, j);
                                    printf("    - sqrt(2*D_%d/Lbd_%d) * sqrt(2*D_%d/Lbd_%d) / [8 sqrt(2 - 2 cos(lbd_%d - lbd_%d))**5] * [\n", i, i, j, j, i, j);
                                    printf("            cos(vrp_%d - vrp_%d + 3 lbd_%d - 3 lbd_%d) + 16 cos(vrp_%d - vrp_%d + 2 lbd_%d - 2 lbd_%d)\n", j, i, i, j, j, i, i, j);
                                    printf("            - 26 cos(vrp_%d - vrp_%d + lbd_%d - lbd_%d) + 9 cos(vrp_%d - vrp_%d + lbd_%d - lbd_%d)\n", j, i, i, j, j, i, j, i);
                                    printf("            ]\n");
                              }
                        }

                        printf(" )\n\n");
                  }
            }
      }
}


int GCD(int * As, int k){

      /******** Returns the gcd of integers *As to *(As + k - 1) ********/

      int i, toBeReturned;
      toBeReturned = *As;
      
      for (i = 1; i < k; i ++){
            toBeReturned = gcd(toBeReturned, *(As + i));
      }
      
      return toBeReturned;
}


int LCM(int * As, int k){

      /******** Returns the lcm of integers *As to *(As + k - 1) ********/

      int i, toBeReturned;
      toBeReturned = *As;
      
      for (i = 1; i < k; i ++){
            toBeReturned = lcm(toBeReturned, *(As + i));
      }
      
      return toBeReturned;
}
