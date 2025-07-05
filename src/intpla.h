#ifndef _INTPLA_H_
#define _INTPLA_H_

#include "parameters.h"
#include "structure.h"

extern typ nu_fast; //The fast frequency
extern typ nu_reso; //The frequency of the periapses
extern typ avgs[1 + 4*how_many_planet]; //Index 4*i - 3 (resp. -2, -1, -0) contains the average of phi_i (resp. v_i, Phi_i, u_i)

void ell2cart(typ a, typ e, typ i, typ nu, typ varpi, typ Omega, typ mu, typ * cart);


void cart2ell(typ * cart, typ * alkhqp, typ mu);


typ mean2true(typ M, typ mu, typ a, typ e);


typ mean2eccentric(typ l, typ k, typ h);


void prods(typ a, typ c, const typ * X, typ b, typ d, const typ * Y, typ * R, typ * S);


void newt(typ DM, typ A, typ B, typ * const p_X, typ * const p_C, typ * const p_S, typ * const EXP, typ * const CM1, typ * const SMX, int bexitonerror);


void kepsaut(typ * cart, typ mu, typ dt);


void exp_tau_LB(typ tau, typ * X_cart);


void UnaveragedSABAn(typ tau, typ T, int output_step, typ * X_old, int n);


void get_frequencies(typ tau, typ T, typ * X_old, int n);


typ UnaveragedSABAn_NAFF(typ tau, typ T, int Hanning_order, typ * X_uv, typ * X_old, int n, int how_many_harmonics);


void Renormalization(typ * X_old);


void get_n(typ * n);


void ConstantParameter(typ * X_new, typ * X_uv);


void LibrationCenterFind(typ * X_old, int precision);


void LibrationCenterFollow(typ * X_old, typ dG, int Npoints, int precision);


void PeriodicOrbitFind(typ * X_old);
#endif
