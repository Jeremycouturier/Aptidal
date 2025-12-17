#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/git/Aptidal/"


/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {1, 1}



/****************************************/
/******** Defining some booleans ********/
/****************************************/
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the Hamiltonian
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.
                             //For example, if max_deg is 3 and one_more_deg_bool is 1, the non-resonant terms will be expanded to 4^th degree. Unimportant if non_resonant_bool is 0.
#define tides_bool 1         //Determines if there are tides raised by the star on the planets in the system. The tidal model of Couturier et al. 2021 is used


                            
/**********************************************************************************************************************************/
/****************************************  Defining some physical constants *******************************************************/
/**********************************************************************************************************************************/
/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** planet and the orbital period of a massless particle with semi-major axis 1 (hence G=4*pi^2). The canonical      ********/
/******** heliocentric coordinates are used (heliocentric positions and barycentric velocities). See Laskar & Robutel 1995 ********/
/**********************************************************************************************************************************/
#define how_many_planet 2                                                           //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964                                               //Gravitational constant is set to 4*pi^2, so the orbital period of the innermost planet is 1
#define body_masses {1., 0.0000729,         0.0000509}                              //Masses of the bodies of the system, beginning with the star.
#define body_sma        {1.,                1.}                                     //Initial and nominal semi-major axes of the planets. Used to computed the Hamiltonian's coef
#define body_ecc        {0.000000000000271, 0.000000000000194}                      //Initial eccentricities              of the planets.
#define body_lambda     {1.048,             0.}                                     //Initial mean longitudes             of the planets.
#define body_varpi      {0.,                0.}                                     //Initial longitudes of the periapses of the planets.
#if tides_bool
#define body_radii      {0.006,             0.006}                                  //Radii of the planets, in units of the innermost planet's semi-major axis
#define body_k2         {1.5,               1.5}                                    //Second Love number of the planets
#define body_Dt         {0.007957747154,    0.007957747154}                         //Tidal timelag of the bodies, in units of the period of a massless particle with semi-major axis 1
#define body_Omega      {6.2831853071796,   6.2831853071796}                        //Initial sideral rotation of the bodies, in radians/period of a massless particle with semi-major axis 1
#define body_alpha      {0.33,              0.33}                                   //Dimensionless structure constant of the bodies. 2/5 for an homogeneous body. 
#endif



/******************************************/
/******** Defining some thresholds ********/
/******************************************/
#define max_deg 3            //Maximum degree in eccentricity. Aptidal allows up to 3.
#define max_res 9            //Maximum value of q for a resonance p : q (with gcd(p,q) = 1 and q >= p). Aptidal currently allows up to 9


#define typ double           //Renaming double as typ. To color gedit, update the field <context id="types" style-ref="type"> of the file c.lang used by gedit
                             //Typical paths where this file is located are /usr/share/gtksourceview-5/language-specs or /usr/share/libgedit-gtksourceview-300/language-specs

#endif
