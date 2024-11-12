#ifndef _INTPLA_H_
#define _INTPLA_H_


#include "parameters.h"
#include "structure.h"


void ell2cart(typ a, typ e, typ i, typ nu, typ varpi, typ Omega, typ mu, typ * cart);


void cart2ell(typ * cart, typ * alkhqp, typ mu);


typ mean2true(typ M, typ mu, typ a, typ e);


void prods(typ a, typ c, const typ * X, typ b, typ d, const typ * Y, typ * R, typ * S);


void newt(typ DM, typ A, typ B, typ * const p_X, typ * const p_C, typ * const p_S, typ * const EXP, typ * const CM1, typ * const SMX, int bexitonerror);


void kepsaut(typ * cart, typ mu, typ dt);


void exp_tau_LB(typ tau, typ * X_cart);


void UnaveragedSABAn(typ tau, typ T, int output_step, typ * X_old, int n);


void get_fast_frequency(typ tau, typ T, int Hanning_order, typ * X_old, int n, typ * nu_2);


void UnaveragedSABAn_NAFF(typ tau, typ T, int Hanning_order, typ * X_uv, typ * X_old, int n, int how_many_harmonics);


void ConstantParameter(typ * X_new, typ * X_uv);


void LibrationCenterFind(typ * X_old, int verbose);

#endif
