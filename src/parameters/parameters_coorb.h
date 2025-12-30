#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/Aptidal_vs_intplabasic/"



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
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the Hamiltonian of the model
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.
                             //For example, if max_deg is 3 and one_more_deg_bool is 1, the non-resonant terms will be expanded to 4^th degree. Unimportant if non_resonant_bool is 0.
#define second_mass_bool 0   //Determines if the model is expanded to second order in mass. Terms at second order in mass are truncated to degree 0 in eccentricity.
#define tides_bool 0         //Determines if there are tides raised by the star on the planets in the system. The tidal model of Couturier et al. 2021 is used
//#define RG_bool 1            //Determines if the first PPN order of General Relativity is included in the Hamiltonian. Eq. (3.55) of https://jeremycouturier.com/img/PhD_manuscript.pdf
#define _3D_bool 0           //For function UnaveragedSABAn only. Determines if the problem is 3D or coplanar
#define canonical_bool 1     //Determines if the heliocentric coordinates below are canonical (heliocentric position, barycentric speed) or non-canonical (heliocentric position, heliocentric speed)
                             //If canonical_bool is 0, the coordinates are converted to canonical for computations, and then converted back to non-canonical before output. DOES NOT WORK. LEAVE TO 1
#define toInvar_bool 1       //Determines if the system is rotated so that the total angular momentum points in the z-direction. Irrelevant if _3D_bool is 0
#define close_enc_bool 0     //Determines if close encounters are checked when integrating the complete system.


                            
/**********************************************************************************************************************************/
/****************************************  Defining some physical constants *******************************************************/
/**********************************************************************************************************************************/
/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** planet and the orbital period of a massless particle with semi-major axis 1 (hence G=4*pi^2). The canonical      ********/
/******** heliocentric coordinates are used (heliocentric positions and barycentric velocities). See Laskar & Robutel 1995 ********/
/**********************************************************************************************************************************/
#define how_many_planet 2                                             //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964                                 //Gravitational constant is set to 4*pi^2, so the orbital period of the innermost planet is 1
#define body_masses {1., 1.337103473310441e-5,1.4290562159203059e-5}  //Masses of the bodies of the system, beginning with the star.
#define body_sma        {1.,                  1.}                     //Initial semi-major axes             of the planets.
#define body_ecc        {0.0046019405927235,  0.009779437763695}      //Initial eccentricities              of the planets.
#define body_lambda     {1.05,                0.}                     //Initial mean longitudes             of the planets in radians.
#define body_varpi      {-4.87256797530106,  -1.73097531149534}       //Initial longitudes of the periapses of the planets in radians.

#if _3D_bool
#define body_inc        {0.1301637290028837,  0.068036468007366,    0.160033715530211,    0.003781880027612,   0.11022781036359}       //Initial inclinations of the planets in radians.
#define body_Omeg       {2.7301739020081763, -1.637829117636709,    3.091766271890028,    0.918782901872667,  -2.36789920376158}       //Initial longitudes of the ascending node of the planets in radians.
#endif

#if tides_bool
#define body_radii      {0.0027,              0.0031}                 //Radii of the planets, in units of length
#define body_k2         {4.5,                 4.5}                    //Second Love number of the planets
#define body_Dt         {0.008376575952205,   0.004973591971621729}   //Tidal timelag of the bodies, in units of time.
#define body_alpha      {0.33,                0.33}                   //Dimensionless structure constant of the bodies. 2/5 for an homogeneous body.
#if _3D_bool
#define body_Omegx      {0.,                  0.,                   0.,                   0.,                  0.}                     //x-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegy      {0.,                  0.,                   0.,                   0.,                  0.}                     //y-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegz      {6.283185307179586477,4.1880552494740315416,2.7912934845404989797,1.860120421559985287,1.2393555199894398798}  //z-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#else
#define body_Omega      {6.283185307179586477,6.283185307179586477}   //Initial sideral rotation of the bodies, in radians/unit of time
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
