#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/Aptidal_simulation/3pla_resonance_exit/K1987/"



/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {8, 12, 18, 27}



/****************************************/
/******** Defining some booleans ********/
/****************************************/
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the Hamiltonian of the model
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.
                             //For example, if max_deg is 3 and one_more_deg_bool is 1, the non-resonant terms will be expanded to 4^th degree. Unimportant if non_resonant_bool is 0.
#define second_mass_bool 0   //Determines if the model is expanded to second order in mass. Terms at second order in mass are truncated to degree 0 in eccentricity.
#define tides_bool 1         //Determines if there are tides raised by the star on the planets in the system. The tidal model of Couturier et al. 2021 is used
#define GR_bool 0            //To be coded. Determines if the first PPN order of General Relativity is included in the Hamiltonian. Eq. (3.55) of https://jeremycouturier.com/img/PhD_manuscript.pdf
#define _3D_bool 0           //For function UnaveragedSABAn only. Determines if the problem is 3D or coplanar
#define canon_input_bool 1   //Determines if the heliocentric coordinates in input  are canonical (heliocentric position, barycentric speed) or not (heliocentric position, heliocentric speed)
#define canon_output_bool 1  //Determines if the heliocentric coordinates in output are canonical (heliocentric position, barycentric speed) or not (heliocentric position, heliocentric speed)
                             //Both canon_input_bool and canon_output_bool must be left to 1 for now
#define toInvar_bool 1       //Determines if the system is rotated so that the total angular momentum points in the z-direction. Irrelevant if _3D_bool is 0
#define close_enc_bool 0     //Determines if close encounters are checked when integrating the complete system.


                            
/**********************************************************************************************************************************/
/****************************************  Defining some physical constants *******************************************************/
/**********************************************************************************************************************************/
/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** planet and the orbital period of a massless particle with semi-major axis the unit of length. The canonical      ********/
/******** heliocentric coordinates are used (heliocentric positions and barycentric velocities). See Laskar & Robutel 1995 ********/
/**********************************************************************************************************************************/
#define how_many_planet 4              //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964  //Gravitational constant. Should be left to 4*pi^2 (the unit of time is the orbital period of a massless particle with semi-major axis the unit of length)
#define m0 1.                          //Mass of the star. Should be left to 1. (the unit of mass is the star mass)
#if GR_bool
#define c_light 1.e9                   //Speed of light in units of length per units of time. Only when GR_bool = 1
#endif
#define body_masses {1.103e-5,            1.195e-5,             2.422e-5,             1.737e-5}             //Masses                              of the planets in stellar masses
#define body_sma    {1.,                  1.31063388956530,     1.71785897200592,     2.25140923008155}     //Initial and nominal semi-major axes of the planets.
#define body_ecc    {.025969963962830640, .059382364145342155,  .037409363918984279,  .029034289643334182}  //Initial eccentricities              of the planets.
#define body_lambda {2.65725525333552,    3.78233356367910,     3.48546755071036,     .14584719428602}      //Initial mean longitudes             of the planets.
#define body_varpi {-0.25011308400757,   -3.39240540126691,    -.25098431811849,     -3.39237472136929}     //Initial longitudes of the periapses of the planets.

#if _3D_bool
#define body_inc    {}       //Initial inclinations of the planets in radians.
#define body_Omeg   {}       //Initial longitudes of the ascending node of the planets in radians.
#endif

#if tides_bool
#define body_radii  {.0016682,            .0018038,             .0036077,             .0026040}             //Radii of the planets, in units of length
#define body_k2     {1.5,                 1.5,                  1.5,                  1.5}                  //Second Love number of the planets
#define body_Dt     {.06578400722471152,  .06845556018452882,   .08839218851765285,   .17204735546952463}   //Tidal timelag of the bodies, in units of time. k2/Q = {6.2e-1, 4.3e-1, 3.7e-1, 4.8e-1}
#define body_alpha  {.33,                 .33,                  .33,                  .33}                  //Dimensionless structure constant of the bodies. 2/5 for an homogeneous body.
#if _3D_bool
#define body_Omegx  {}                     //x-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegy  {}                     //y-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegz  {}                     //z-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#else
#define body_Omega  {6.283219958851004,   4.187553542589783,    2.790640330614937,    1.8599530293662827}   //Initial sideral rotation of the bodies, in radians/unit of time
#endif
#endif



/******************************************/
/******** Defining some thresholds ********/
/******************************************/
#define max_deg 3            //Maximum degree in eccentricity. Aptidal allows up to 3.
#define max_res 9            //Maximum value of q for a resonance p : q (with gcd(p,q) = 1 and q >= p). Aptidal currently allows up to 9


#define typ double           //Renaming double as typ. To color gedit, update the field <context id="types" style-ref="type"> of the file c.lang used by gedit
                             //Typical paths where this file is located are /usr/share/gtksourceview-5/language-specs or /usr/share/libgedit-gtksourceview-300/language-specs

#if _3D_bool
#define Nd 6
#else
#define Nd 4
#endif

#endif
