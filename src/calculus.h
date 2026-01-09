#ifndef _CALCULUS_H_
#define _CALCULUS_H_

#include "parameters.h"
#include "structure.h"

extern typ X_old_t0[4*how_many_planet + 1]; 
extern typ X_new_t0[4*how_many_planet + 1];
extern typ X_uv_t0 [4*how_many_planet + 1];


typ fast_pow(typ x, int power);


void dHdold(typ * dH, typ * X_old, int KP);


void dHdnew(typ * dH_polar, typ * dH_rect, typ * dH_old, typ * X_new, typ * X_uv);


void old2new(typ * X_old, typ * X_new, typ * X_uv);


void new2old(typ * X_old, typ * X_new, typ * X_uv);


void canonical2nonCanonical(typ * X_cart);


void nonCanonical2canonical(typ * X_cart);


#if (toInvar_bool && _3D_bool)
void toInvar(typ * X_cart);
#endif


void X_old_init(typ * X_old);


void nonDofReinit(typ * X_new, typ * X_uv);


void exp_tau_LB_Ralston(typ tau, typ * X_old);


void AveragedSABAn(typ tau, typ T, int output_step, typ * X_old, int n);


void SABAn_average(typ tau, typ T, int Hanning_order, typ * X_uv_mean, typ * X_old, int n);


void RK2(typ tau, typ T, int output_step);


typ AveragedHamiltonian(typ * X_old);


void PointPrint(typ * X_old, int iter);


int EquilibriumFind(typ * X_old, int precision);


void EquilibriumFollow(typ * X_old, typ dG, int Npoints, int precision);


void EquilibriumFindUntil(typ * X_old, int precision);
#endif
