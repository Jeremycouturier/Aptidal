#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_


#define pth "/home/jeremy/Documents/Aptidal_simulation/test/"


/**********************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n.  ********/
/**********************************************************************/
#define resonance_chain   {1, 1, 2}



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
#define how_many_planet 3                                                           //Number of orbiting bodies in the system
#define G 39.478417604357434475337964                                               //Gravitational constant is set to 4*pi^2, so the orbital period of the innermost planet is 1
#define body_masses {1.,   0.0000184, 0.0000108,    0.0000212}                      //Masses of the bodies of the system, beginning with the star.
#define body_sma          {1.,        1.,           1.5874010520}                   //Initial and nominal semi-major axes of the orbiting bodies. Used to computed the Hamiltonian's coef
#define body_ecc          {0.02,      0.02,         0.02}                           //Initial eccentricities              of the orbiting bodies.
#define body_lambda      {-1.2678,    0.27182,     -2.8035}                         //Initial mean longitudes             of the orbiting bodies.
#define body_varpi       {-1.0614,    0.81578,     -3.07126}                        //Initial longitude of the periapsis  of the orbiting bodies.
#define body_radii  {0.12, 0.006744,  0.010036}                                     //Radii               of the bodies of the system, beginning with the star. Unimportant if no tides
#define body_k2     {0.0,  0.8,       1.32}                                         //Second Love numbers of the bodies of the system, beginning with the star. Unimportant if no tides
#define body_Q  {10000.0,  234.8,     78.9}                                         //Quality factors     of the bodies of the system, beginning with the star. Unimportant if no tides



/******************************************/
/******** Defining some thresholds ********/
/******************************************/
#define max_deg 1            //Maximum degree in eccentricity. Aptidal allows up to 3.
#define max_res 9            //Maximum value of q for a resonance p:q (with gcd(p,q) = 1). Aptidal currently allows up to 9


#define typ double           //Renaming double as typ. To color gedit : Update field <context id="types" style-ref="type"> of /usr/share/gtksourceview/language-specs/c.lang

#endif
