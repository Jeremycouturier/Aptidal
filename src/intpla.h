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

#endif
