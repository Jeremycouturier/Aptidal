#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include "parameters.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) < (b) ? (b) : (a))

struct pairOfReal {                          //Simple data structure defining a pair of real numbers
      typ fst;
      typ snd;
};


struct rational {                            //Simple data structure representing a rational number
      int numerator;
      int denominator;
};


#if toInvar_bool
/******** Defining a quaternion structure to perform rotations from one vector to another ********/
struct quaternion {
      typ w;
      typ x;
      typ y;
      typ z;
};
#endif


/******** Parameters of the system ********/
extern typ masses[how_many_planet + 1];      //Masses                            of the planets
extern typ sma   [how_many_planet + 1];      //Semi-major axes                   of the planets
extern typ ecc   [how_many_planet + 1];      //Eccentricities                    of the planets
extern typ lbd   [how_many_planet + 1];      //Mean longitudes                   of the planets
extern typ vrp   [how_many_planet + 1];      //Longitudes of the periapses       of the planets
#if _3D_bool
extern typ inc   [how_many_planet + 1];      //inclinations                      of the planets
extern typ Om    [how_many_planet + 1];      //Longitudes of the ascending nodes of the planets
#endif
extern typ Lbd_0 [how_many_planet + 1];      //Nominal circular angular momenta  of the planets
#if tides_bool
extern typ radii [how_many_planet + 1];      //Radii                             of the planets
extern typ k2s   [how_many_planet + 1];      //Second Love numbers               of the planets
extern typ Dts   [how_many_planet + 1];      //Tidal timelags                    of the planets
extern typ alps  [how_many_planet + 1];      //Dimensionless structure constants of the planets
#if _3D_bool
extern typ Omx   [how_many_planet + 1];      //x-coordinate of sideral rotations of the planets
extern typ Omy   [how_many_planet + 1];      //y-coordinate of sideral rotations of the planets
extern typ Omz   [how_many_planet + 1];      //z-coordinate of sideral rotations of the planets
#else
extern typ Omg   [how_many_planet + 1];      //Sideral rotations                 of the planets
#endif
#endif
extern typ m0;                               //Mass of the star


/******** Some global variables and arrays ********/
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
extern int how_many_resonant;                 //number of planets in the system that are resonant     and are not the first member of a co-orbital pair
extern int how_many_missed;                   //number of planets in the system that are not resonant or that are the first member of a co-orbital pair
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


struct rational real2rat(typ r, int maxDenom);


struct rational ratabs(struct rational r1);


void ratprint(struct rational r);


typ rdm(typ min, typ max);


typ continuousAngle(typ newAngle, typ oldAngle);


#if toInvar_bool
void quaternion_norm(struct quaternion * q);


struct quaternion get_quaternion(typ ux, typ uy, typ uz, typ vx, typ vy, typ vz);


void rotate_with_quaternion(typ x, typ y, typ z, struct quaternion q, typ * xr, typ * yr, typ * zr);
#endif
#endif
