#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/path/towards/output/location/with/slash/at/the/end/"



/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {1, 2, 4, 8}



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
#define body_masses  {3.105962448956341e-5, 3.802386494804825e-5, 2.75622879221497e-5, 2.463327515053729e-5}  //Masses                              of the planets.
#define body_sma     {1.,                   1.588616530008093,    2.522519965126477,   4.002895283219303}     //Initial and nominal semi-major axes of the planets.
#define body_ecc     {.2037366402560644,    .23327738738722253,   .10879789578159746,  .010608528803873652}   //Initial eccentricities              of the planets.
#define body_lambda  {1.4002712016807659,   1.2232899783393805,   .3829287109240201,  -.3338275493247094}     //Initial mean longitudes             of the planets.
#define body_varpi   {1.5685233484435790,  -.222635740123498,    -1.1583119628385794,  .6315248701576947}     //Initial longitudes of the periapses of the planets.

#if _3D_bool
#define body_inc     {}       //Initial inclinations of the planets in radians.
#define body_Omeg    {}       //Initial longitudes of the ascending node of the planets in radians.
#endif

#if tides_bool
#define body_radii   {.008613817699103963, .010503474561304208,  .007661821216636959,  .006862709794492351}   //Radii of the planets, in units of length
#define body_k2      {.15,                 .15,                  .15,                  .15}                   //Second Love number of the planets
#define body_Dt      {.0001909955049112505,.00023359884046906003,.00034024622326981984,.0011046589597291702}  //Tidal timelag of the planets, in units of time. k21/Q1 = 1.8e-4, k22/Q2 = 1.1e-4, k23/Q3 = 0.8e-4, k24/Q4 = 1.3e-4
#define body_alpha   {.33,                 .33,                  .33,                  .33}                   //Dimensionless momoent of inertia of the planets. 2/5 for an homogeneous body.
#if _3D_bool
#define body_Omegx   {}                     //x-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegy   {}                     //y-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegz   {}                     //z-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#else
#define body_Omega   {6.283282883110044,   3.13804747150582,     1.568317306533526,    .7845558658928975}     //Initial sideral rotation of the planets, in radians/unit of time
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
