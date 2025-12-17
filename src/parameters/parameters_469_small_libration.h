#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/git/K2138/K2138_in_3pla_0th_degree/Numerical_integration/data/"


/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {4, 6, 9}//{12, 15, 20} //{4, 6, 9}



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
#define how_many_planet 3                                                           //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964                                               //Gravitational constant is set to 4*pi^2, so the orbital period of the innermost planet is 1
#define body_masses {1., 1.337103473310441e-5,1.4290562159203059e-5,2.1684355013014528e-5}   //Masses of the bodies of the system, beginning with the star.
#define body_sma        {1.00003707373117,  1.34656616289081,  1.84195039398514}    //Initial and nominal semi-major axes of the planets. Used to computed the Hamiltonian's coef
#define body_ecc        {0.0002366091964394,0.0004810360841520,0.00013722998705962} //Initial eccentricities              of the planets.
#define body_lambda     {2.18222340496659,  2.50181742918604,  5.85666849075520}    //Initial mean longitudes             of the planets.
#define body_varpi     {-3.07901689658590,  0.05141267463016, -3.07903197012034}    //Initial longitudes of the periapses of the planets.
/*#define body_masses {1., 0.00001396343247, 0.00001018846391, 0.00001168182951}      //Masses of the bodies of the system, beginning with the star.
#define body_sma        {1.,               1.1761317071578,  1.45803441881}         //Initial and nominal semi-major axes of the planets. Used to computed the Hamiltonian's coef
#define body_ecc        {0.00349902571220, 0.00633968349406, 0.001088416895255772}  //Initial eccentricities              of the planets.
#define body_lambda     {3.24613039776969, 2.00403931789809, 0.00767556394251}      //Initial mean longitudes             of the planets.
#define body_varpi      {0.27785699707616, 3.41403973628331, 0.26235900507790}      //Initial longitudes of the periapses of the planets.*/
#if tides_bool
#define body_radii      {0.002083879934854, 0.00281041371739,  0.0031803169951999}  //Radii of the planets, in units of the innermost planet's semi-major axis
#define body_k2         {4.5,               4.5,               4.5}                 //Second Love number of the planets
#define body_Dt         {0.005684105110,    0.008596531650,    0.013056447233}      //Tidal timelag of the bodies, in units of the period of a massless particle with semi-major axis 1
#define body_Omega      {6.2831853071796,   4.1544994152875,   2.7353754873595}     //Initial sideral rotation of the bodies, in radians/period of a massless particle with semi-major axis 1
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
