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
#define resonance_chain {0, 0, 0, 0, 0, 0, 0}



/****************************************/
/******** Defining some booleans ********/
/****************************************/
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the Hamiltonian of the model
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.
                             //For example, if max_deg is 3 and one_more_deg_bool is 1, the non-resonant terms will be expanded to 4^th degree. Unimportant if non_resonant_bool is 0.
#define second_mass_bool 0   //Determines if the model is expanded to second order in mass. Terms at second order in mass are truncated to degree 0 in eccentricity.
#define tides_bool 1         //Determines if there are tides raised by the star on the planets in the system. The tidal model of Couturier et al. 2021 is used
#define GR_bool 0            //To be coded. Determines if the first PPN order of General Relativity is included in the Hamiltonian. Eq. (3.55) of https://jeremycouturier.com/img/PhD_manuscript.pdf
#define _3D_bool 1           //For function UnaveragedSABAn only. Determines if the problem is 3D or coplanar
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
#define how_many_planet 7              //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964  //Gravitational constant. Should be left to 4*pi^2 (the unit of time is the orbital period of a massless particle with semi-major axis the unit of length)
#define m0 1.                          //Mass of the star. Should be left to 1. (the unit of mass is the star mass)
#if GR_bool
#define c_light 1.e9                   //Speed of light in units of length per units of time. Only when GR_bool = 1
#endif
#define body_masses  {3.4441111046788264e-6, 5.006385307140297e-6, 5.098523511875069e-6, 1.2100080582565424e-5, 1.8110574767183746e-5, 1.0163623683153116e-5, 1.1231554910743344e-5} //Masses                              of the planets.
#define body_sma     {1.,                    1.405849712250732,    1.841685550334578,    2.4166699708090813,    3.1649937095083525,    3.8278412697379145,    5.396333153416171}     //Initial and nominal semi-major axes of the planets.
#define body_ecc     {.024811735891496517,   .07766009377916062,   .12162524097221987,   .04966311478765057,    .033377725032406276,   .07445790010177279,    .022427897517380543}   //Initial eccentricities              of the planets.
#define body_lambda {-2.2528309377327713,   -2.1339982888798703,   2.3878744905464764,  -1.7676471348488083,   -1.4659673999363996,   -.24797123764044776,    2.2339884324328354}    //Initial mean longitudes             of the planets.
#define body_varpi   {2.528221268504456,    -.8843142085925679,   -3.0174410561598943,   .659786403715292,      2.889181913716233,     .10876356574395767,   -.541762350471001}      //Initial longitudes of the periapses of the planets.

#if _3D_bool
#define body_inc     {.17547906999312193,    .11172037206374609,   .006722428448706527,  .03143059339382074,    .013591145845092552,   .022464498425872937,   .01208994631651904}    //Initial inclinations of the planets in radians.
#define body_Omeg   {-0.42887660077686823,   2.2985517298239966,   .5735373770827941,   -.874560860550921,      1.6569697571984747,   -2.956609751405752,    -1.2753958772284055}    //Initial longitudes of the ascending node of the planets in radians.
#endif

#if tides_bool
#define body_radii   {.0018494378572032591,  .002668711479582546,  .002716855133060715,  .006339436191818241,   .009413701223667024,   .005343133689986821,   .005893001327721096}   //Radii of the planets, in units of length
#define body_k2      {1.5,                   1.5,                  1.5,                  1.5,                   1.5                    1.5                    1.5}                   //Second Love number of the planets
#define body_Dt      {.004138021394492695,   .014499986495027591,  .007695873853896611,  .02031627039579169,    .028079028306872273,   .010330013725938593,   .11571610772595842}    //Tidal timelag of the planets, in units of time. k21/Q1 = 3.9e-2, k22/Q2 = 8.2e-2, k23/Q3 = 2.9e-2, k24/Q4 = 5.1e-2, k25/Q5 = 4.7e-2, k26/Q6 = 1.3e-2, , k27/Q7 = 8.7e-2
#define body_alpha   {.33,                   .33,                  .33,                  .33,                   .33,                   .33,                   .33}                   //Dimensionless momoent of inertia of the planets. 2/5 for an homogeneous body.
#if _3D_bool
#define body_Omegx   {0.,                    0.                    0.                    0.                     0.                     0.                     0.}                    //x-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegy   {0.,                    0.                    0.                    0.                     0.                     0.                     0.}                    //y-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegz   {6.283196127164415,     3.76940333938465,     2.5139564757683734,   1.6724628412250206,    1.115898064238376,     .8389792014414005,     .5012266756963253}     //z-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#else
#define body_Omega   {} //Initial sideral rotation of the planets, in radians/unit of time
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
