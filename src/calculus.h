#ifndef _CALCULUS_H_
#define _CALCULUS_H_


#include "parameters.h"
#include "structure.h"



extern typ grad_old[4 * how_many_planet + 1]; //The gradient of the Hamiltonian in the         old variables (lbd, -vrp; Lbd, D)
extern typ * grad_new;                        //The gradient of the Hamiltonian in the         new variables (phi, sig;  Phi, D)
extern typ * grad_uv;                         //The gradient of the Hamiltonian in the rectangular variables (phi, v;    Phi, u)
extern typ X_old   [4 * how_many_planet + 1]; //The old         variables (lbd, -vrp; Lbd, D). Indexes 4*i-3 to 4*i contain lbd_i to D_i
extern typ * X_new;                           //The new         variables (phi, sig;  Phi, D). Indexes 4*i-3 to 4*i contain phi_i to D_i
extern typ * X_uv;                            //The rectangular variables (phi, v;    Phi, u). Indexes 4*i-3 to 4*i contain phi_i to u_i



typ fast_pow(typ x, int power);


void dH_old(int K, typ * dHk);


void gradH_old();


typ average_test(typ A, typ dt, typ T);























#endif
