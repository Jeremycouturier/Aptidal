#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/git/K2138/K2138_in_3pla_0th_degree/Numerical_integration/data/5pla_present/"


/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {16, 24, 36, 54, 81, 243}//{16, 24, 36, 54, 81, 243}



/****************************************/
/******** Defining some booleans ********/
/****************************************/
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the Hamiltonian
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.
                             //For example, if max_deg is 3 and one_more_deg_bool is 1, the non-resonant terms will be expanded to 4^th degree. Unimportant if non_resonant_bool is 0.
#define tides_bool 0         //Determines if there are tides raised by the star on the planets in the system. The tidal model of Couturier et al. 2021 is used


                            
/**********************************************************************************************************************************/
/****************************************  Defining some physical constants *******************************************************/
/**********************************************************************************************************************************/
/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** planet and the orbital period of a massless particle with semi-major axis 1 (hence G=4*pi^2). The canonical      ********/
/******** heliocentric coordinates are used (heliocentric positions and barycentric velocities). See Laskar & Robutel 1995 ********/
/**********************************************************************************************************************************/
#define how_many_planet 6                                                           //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964                                               //Gravitational constant is set to 4*pi^2, so the orbital period of the innermost planet is 1
#define body_masses {1., 1.337103473310441e-5,1.4290562159203059e-5,2.1684355013014528e-5,3.122664547488287e-5,1.7893159976480633e-5,2.0354248e-5}   //Masses of the bodies of the system, beginning with the star.
#define body_sma        {1.000004456991713, 1.317710524299074, 1.7407583104386735, 2.3102117987653443, 3.0865771944685254, 6.826756893864246}    //Initial and nominal semi-major axes of the planets. Used to computed the Hamiltonian's coef
#define body_ecc        {0.0004834870730296,0.00049117635460117,0.00030839053219745046,0.006405488954639363,1.6004960508244145e-6,0.0005367253195170454} //Initial eccentricities              of the planets.
#define body_lambda     {0.7751686256869148,3.2533512640111515,1.5302110617377698, 5.631820995206427,  1.0831604767994083, 2.4426447224951158}  //Initial mean longitudes             of the planets.
#define body_varpi      {0.2394931999756539,2.882371561867995,-1.8105969411716547,-1.5707513360222807,-0.8699367319193446,-0.010223777925911098}  //Initial longitudes of the periapses of the planets.
#if tides_bool
#define body_radii      {0.002083879934854, 0.00281041371739,  0.0031803169951999}  //Radii of the planets, in units of the innermost planet's semi-major axis
#define body_k2         {4.5,               4.5,               4.5}                 //Second Love number of the planets
//#define body_Dt         {0.0056840250358102,0.0085286276976171,0.01279877379834141} //Tidal timelag of the bodies, in units of the period of a massless particle with semi-major axis 1  // Q=28
#define body_Dt         {0.028420125179051, 0.0426431384880855,0.06399386899170705} //Tidal timelag of the bodies, in units of the period of a massless particle with semi-major axis 1  // Q=5.6
#define body_Omega      {6.2832738225606,   4.187577061695874, 2.790445887786056}   //Initial sideral rotation of the bodies, in radians/period of a massless particle with semi-major axis 1
#define body_alpha      {0.33,              0.33,              0.33}                //Dimensionless structure constant of the bodies. 2/5 for an homogeneous body. 
#endif



/******************************************/
/******** Defining some thresholds ********/
/******************************************/
#define max_deg 3            //Maximum degree in eccentricity. Aptidal allows up to 3.
#define max_res 9            //Maximum value of q for a resonance p : q (with gcd(p,q) = 1 and q >= p). Aptidal currently allows up to 9


#define typ double           //Renaming double as typ. To color gedit, update the field <context id="types" style-ref="type"> of the file c.lang used by gedit
                             //Typical paths where this file is located are /usr/share/gtksourceview-5/language-specs or /usr/share/libgedit-gtksourceview-300/language-specs

#endif
