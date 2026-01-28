#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/git/K2138/K2138_in_3pla_0th_degree/Numerical_integration/data/5pla_past/k2Q=1.4_0.5_0.5_0.5/"
//#define pth "/home/atipique/Documents/Aptidal_vs_intplabasic/Hamiltonian_conservation_Aptidal/"
//#define pth "/home/atipique/Documents/Aptidal_vs_intplabasic/Canonical/"



/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {16, 24, 36, 54, 81}
//#define resonance_chain {0, 0, 0, 0, 0}



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
#define how_many_planet 5              //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964  //Gravitational constant. Should be left to 4*pi^2 (the unit of time is the orbital period of a massless particle with semi-major axis the unit of length)
#define m0 1.                          //Mass of the star. Should be left to 1. (the unit of mass is the star mass)
#if GR_bool
#define c_light 1.e9                   //Speed of light in units of length per units of time. Only when GR_bool = 1
#endif
#define body_masses {1.337103473310441e-5,1.4290562159203059e-5,2.1684355013014528e-5,3.122664547488287e-5,1.7893159976480633e-5}  //Masses of the planets in units of the star mass.
#define body_sma    {1.,                  1.31052399590020,     1.71757711268338,     2.25126112381953,    2.95113657525905}       //Initial semi-major axes             of the planets.
#define body_ecc    {0.04601940592723566, 0.0977943776369584,   0.0774781222725785,   0.0464169646657162,  0.029766096365539559}   //Initial eccentricities              of the planets.
#define body_lambda {0.73461451389486,    3.05434389145679,     1.45923749525822,     5.63182099520643,    1.08316047679941}       //Initial mean longitudes             of the planets in radians.
#define body_varpi  {-4.87256797530106,  -1.73097531149534,    -4.87256793692540,    -1.73097525946379,   -4.87256789885291}       //Initial longitudes of the periapses of the planets in radians.

#if _3D_bool
#define body_inc    {0.1301637290028837,  0.068036468007366,    0.160033715530211,    0.003781880027612,   0.11022781036359}       //Initial inclinations of the planets in radians.
#define body_Omeg   {2.7301739020081763, -1.637829117636709,    3.091766271890028,    0.918782901872667,  -2.36789920376158}       //Initial longitudes of the ascending node of the planets in radians.
#endif

#if tides_bool
#define body_radii  {0.0020838706470602,  0.002810401191461818, 0.003180302820616584, 0.004175347043679761,0.0036929697668095203}  //Radii of the planets, in units of length
//#define body_k2     {1.,                  1.,                   1.,                   1.,                  1.}                     //Second Love number of the planets
#define body_k2     {14.,                 5.,                   5.,                   5.,                  0.}                     //Second Love number of the planets
//#define body_Dt     {0.001591549430918,   0.00238774309418575,  0.00358256845988597,  0.00537599602912457, 0.008068709775936776}   //Tidal timelag of the bodies, in units of time. Q=100
#define body_Dt     {0.01591549430918,    0.0238774309418575,   0.0358256845988597,   0.0537599602912457,  0.08068709775936776}    //Tidal timelag of the bodies, in units of time. Q=10
#define body_alpha  {0.33,                0.33,                 0.33,                 0.33,                0.33}                   //Dimensionless structure constant of the bodies. 2/5 for an homogeneous body.
#if _3D_bool
#define body_Omegx  {0.,                  0.,                   0.,                   0.,                  0.}                     //x-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegy  {0.,                  0.,                   0.,                   0.,                  0.}                     //y-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegz  {6.283185307179586477,4.1880552494740315416,2.7912934845404989797,1.860120421559985287,1.2393555199894398798}  //z-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#else
#define body_Omega  {6.283185307179586477,4.1880552494740315416,2.7912934845404989797,1.860120421559985287,1.2393555199894398798}  //Initial sideral rotation of the bodies, in radians/unit of time
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
