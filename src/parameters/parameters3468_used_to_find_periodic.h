#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/Aptidal_simulation/towards_right_less_precise/"


/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {3, 4, 6, 8}



/****************************************/
/******** Defining some booleans ********/
/****************************************/
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the Hamiltonian
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.
                             //For example, if max_deg is 3 and one_more_deg_bool is 1, the non-resonant terms will be expanded to 4th degree. Unimportant if non_resonant_bool is 0.


                            
/**********************************************************************************************************************************/
/********                                  Defining some physical constants                                                ********/
/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** orbiting body and the initial orbital period of the innermost orbiting body (hence G=4*pi^2). The canonical      ********/
/******** heliocentric coordinates are used (heliocentric positions and barycentric velocities). See Laskar & Robutel 1995 ********/
/**********************************************************************************************************************************/
#define how_many_planet 4                                                           //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964                                               //Gravitational constant is set to 4*pi^2, so the orbital period of the innermost planet is 1
#define body_masses {1.,   0.000019756, 0.000013616, 0.000021358, 0.000012815}      //Masses of the bodies of the system, beginning with the star.
#define body_sma          {1.,          1.215717260, 1.604470493, 1.957775389}      //Initial and nominal semi-major axes of the planets. Used to computed the Hamiltonian's coef
#define body_ecc          {0.001825403, 0.002563375, 0.003476109, 0.000907217}      //Initial eccentricities              of the planets.
#define body_lambda       {1.88022,     2.0087,      0.00,        2.00827}          //Initial mean longitudes             of the planets.
#define body_varpi        {2.71880,    -1.2602,     -2.9007,     -3.0981}           //Initial longitude of the periapsis  of the planets.
#define body_radii  {0.12, 0.006744,    0.010036}                                   //Radii               of the planets, beginning with the star. Unimportant if no tides
#define body_k2     {0.0,  0.8,         1.32}                                       //Second Love numbers of the planets, beginning with the star. Unimportant if no tides
#define body_Q  {10000.0,  234.8,       78.9}                                       //Quality factors     of the planets, beginning with the star. Unimportant if no tides



/******************************************/
/******** Defining some thresholds ********/
/******************************************/
#define max_deg 3            //Maximum degree in eccentricity. Aptidal allows up to 3.
#define max_res 9            //Maximum value of q for a resonance p : q (with gcd(p,q) = 1 and q >= p). Aptidal currently allows up to 9


#define typ double           //Renaming double as typ. To color gedit : Update field <context id="types" style-ref="type"> of /usr/share/gtksourceview/language-specs/c.lang

#endif
