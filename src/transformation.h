#ifndef _TRANSFORMATION_H_
#define _TRANSFORMATION_H_

#include "parameters.h"
#include "structure.h"

extern struct rational transformation[2*how_many_planet + 1][2*how_many_planet + 1]; //matrix M such that (phi, sigma) = M(lambda, -varpi)
extern struct rational transpose_inv [2*how_many_planet + 1][2*how_many_planet + 1]; //tM^-1, that is, (Phi,D) = (tM^-1)(Lambda, D)
extern typ             Transformation[2*how_many_planet + 1][2*how_many_planet + 1]; //matrix M such that (phi, sigma) = M(lambda, -varpi)
extern typ             Transpose_inv [2*how_many_planet + 1][2*how_many_planet + 1]; //tM^-1, that is, (Phi,D) = (tM^-1)(Lambda, D)
extern struct rational rat_l_ij      [  how_many_planet - 1][  how_many_planet + 1]; //Coefficients l_ij
extern struct rational NoverD        [  how_many_planet - 1][  how_many_planet - 1][how_many_planet + 1]; //Coefficients N_r,s^(i) and D_r,s^(i)
extern struct rational rat_c_i       [  how_many_planet - 1];                        //Coefficients c_i
extern int dof                       [  how_many_planet + 1];                        //The indexes j of the degrees of freedom (phi_j; Phi_j)
extern int how_many_dof;                                                             //The number of degrees of freedom of the form (phi_j; Phi_j)
extern int nondof                    [  how_many_planet + 1];                        //The indexes j of the non-degrees of freedom (phi_j; Phi_j)
extern int how_many_nondof;                                                          //The number of non-degrees of freedom of the form (phi_j; Phi_j)


void l_ij_init();


void NoverD_init();


void rat_c_i_init();


void transformation_init();


void transpose_inv_init();


void verification();


void matrix_fill();


void transformation_display();


void Hamiltonian_display();


int GCD(int * As, int k);


int LCM(int * As, int k);
#endif
