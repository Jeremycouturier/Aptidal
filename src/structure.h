#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include "parameters.h"


struct pairOfReal {                           //Simple data structure defining a pair of real numbers
      typ fst;
      typ snd;
};


struct rational {                             //Simple data structure representing a rational number
      int numerator;
      int denominator;
};


/******** Parameters of the system ********/
extern typ masses[how_many_planet + 1];      //Masses                           of the planets
extern typ sma   [how_many_planet + 1];      //Semi-major axes                  of the planets
extern typ ecc   [how_many_planet + 1];      //Eccentricities                   of the planets
extern typ lbd   [how_many_planet + 1];      //Mean longitudes                  of the planets
extern typ vrp   [how_many_planet + 1];      //Longitudes of the ascending node of the planets
extern typ m0;


/******** Some global arrays ********/
extern int p_i     [how_many_planet + 1];
extern int k_ij    [how_many_planet + 1][how_many_planet + 1];
extern int q_ij    [how_many_planet]    [how_many_planet + 1];
extern int k_i     [how_many_planet];
extern int k_star_i[how_many_planet];
extern int q_i     [how_many_planet];
extern struct rational rat_k_ij    [how_many_planet + 1][how_many_planet + 1];
extern struct rational rat_q_ij    [how_many_planet]    [how_many_planet + 1];
extern struct rational rat_k_i     [how_many_planet];
extern struct rational rat_k_star_i[how_many_planet];
extern struct rational rat_q_i     [how_many_planet];
extern int how_many_resonant;                 //number of planets in the resonance chain that are not unresonant of the first member of a co-orbital pair
extern int how_many_missed;                   //number of planets in the resonance chain that are unresonant of that are the first member of a co-orbital pair
extern int * subchain;
extern int * chain_missed;
extern typ qnmlr[32][5];                        



void init();


void initialization();


void chain_validity();


void array_init();


void array2rational();


void deallocation();


struct pairOfReal b_12(typ alpha, typ precision);


struct pairOfReal b_k2(typ alpha, struct pairOfReal b_km12, int k);


int gcd(int a, int b);


int lcm(int a, int b);


struct rational ratmul(struct rational r1, struct rational r2);


struct rational ratadd(struct rational r1, struct rational r2);


struct rational ratopp(struct rational r1);


struct rational ratinv(struct rational r1);


struct rational ratdiff(struct rational r1, struct rational r2);


struct rational ratover(struct rational r1, struct rational r2);


struct rational int2rat(int a);


typ rat2real(struct rational r1);


struct rational ratabs(struct rational r1);


void ratprint(struct rational r);


#endif
