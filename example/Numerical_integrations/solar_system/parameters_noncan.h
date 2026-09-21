#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/home/atipique/Documents/git/Aptidal/example/Numerical_integrations/solar_system/"



/*******************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers  ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios ********/
/******** Two or more consecutive p_i can be equal (coorbital planets) and p_i = 0  ********/
/******** indicates that planet n° i is not in resonance with other planets         ********/
/*******************************************************************************************/
#define resonance_chain {0, 0, 0, 0, 0, 0, 0, 0}



/****************************************/
/******** Defining some booleans ********/
/****************************************/

/******** Boolean relative to numerical simulations ********/
#define tides_bool 0         //Determines if there are tides raised by the star on the planets in the system. See Couturier et al. (2026), Appendix B.2.3.
#define GR_bool 0            //Determines if the first PPN order of General Relativity is included. See Saha&Tremaine (1994), Sect. 5. Implemented but not yet verified. Could contain errors.
#define _3D_bool 1           //Determines if the problem is 3D or coplanar
#define toInvar_bool 1       //Determines if the system is rotated so that the total angular momentum points in the z-direction. Irrelevant if _3D_bool is 0
#define close_enc_bool 0     //Determines if close encounters are checked when integrating the complete system. Does not handle large timesteps.

/******** Boolean relative to coordinate systems ********/
#define canon_input_bool 0   //Determines if the heliocentric coordinates in input  are canonical (heliocentric position, barycentric speed) or not (heliocentric position, heliocentric speed)
#define canon_output_bool 0  //Determines if the heliocentric coordinates in output are canonical (heliocentric position, barycentric speed) or not (heliocentric position, heliocentric speed)
#define ellip_input_bool 0   //Determines if the input  coordinates are elliptic or cartesian
#define ellip_output_bool 1  //Determines if the output coordinates are elliptic or cartesian

/******** Boolean relative to the averaged model ********/
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the averaged Hamiltonian
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.

                            
/**********************************************************************************************************************************/
/****************************************  Defining some physical constants *******************************************************/
/**********************************************************************************************************************************/
/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** planet and the orbital period of a massless particle with semi-major axis the unit of length.                    ********/
/**********************************************************************************************************************************/
#define how_many_planet 8              //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G .0002959122                  //Gravitational constant. Should be left to 4*pi^2 (the unit of time is the orbital period of a massless particle with semi-major axis the unit of length)
#define m0 1.                          //Mass of the star. Should be left to 1. (the unit of mass is the star mass)
#if GR_bool
#define c_light 173.144632674          //Speed of light in units of length per units of time. Only when GR_bool = 1. Here in AU/days
#endif
#define body_masses {1.660120825459e-7, 2.447838287797e-6, 3.040432648963e-6, 3.227156082932e-7, 9.547919099414e-4, 2.858856700246e-4, 4.366249613222e-5, 5.151383772654e-5} //Masses in Sun masses
#define body_sma    { .02136639564901786,   .00079811748145842455,-.017203109059379617, .00067149954213035459,-.0045683137846838275,-.0042923515984516171,.0026781050868787204,.0025792748075360731} //vx
#define body_ecc    {-.0049262994196977248,-.018491837535310501,  -.0029028419967796386,.013814037515929691,   .0058814620366561881, .0035283450596930352,.0024620045387926306,.0016684247014475172} //vy
#define body_lambda {-.13009360552399035,  -.71830229577683646,   -.17715878283537426,  1.3907159218099714,    4.0011770363934236,  6.4064093505456867,   14.431859127944875, 16.812053232861203}    //x
#define body_varpi  {-.40059371556854023,  -.046274247988195134,   .88740685950735376,  .0014012156287008626,  2.7365789474791216,  6.1746575459511632,  -12.506266486666524,-22.980107302845063}    //y

#if _3D_bool
#define body_inc    {-.0048474335356574125,-.0083697352411685674, -.0012585096152846315,.0063179004300907187,  .0026323026658345568,.0016419307685794325, .0010404080796403874,.00061881380327439954}//vz
#define body_Omeg   {-.200489313264401,     .024640644386509258,   .3847367177944262,  -.036960167462523869,   1.075512278219606,   2.2747711011338763,  -5.6816893108481299, -9.8244261231785899}   //z
#endif

#if tides_bool
#define body_radii   {}  //Radii of the planets, in units of length
#define body_k2      {}  //Second Love number of the planets
#define body_Dt      {}  //Tidal timelag of the planets, in units of time. k21/Q1 = 4.6e-2, k22/Q2 = 7.2e-2, k23/Q3 = 3.9e-2
#define body_alpha   {}  //Dimensionless momoent of inertia of the planets. 2/5 for an homogeneous body.
#if _3D_bool
#define body_Omegx   {}  //x-coordinate of the initial sideral rotation of the planets, in radians/unit of time
#define body_Omegy   {}  //y-coordinate of the initial sideral rotation of the planets, in radians/unit of time
#define body_Omegz   {}  //z-coordinate of the initial sideral rotation of the planets, in radians/unit of time
#else
#define body_Omega   {}  //Initial sideral rotation of the planets, in radians/unit of time
#endif
#endif



/******************************************/
/******** Defining some thresholds ********/
/******************************************/
#define max_deg 3            //Maximum degree in eccentricity for the averaged model. Aptidal allows up to 3.
#define max_res 9            //Maximum value of q for a resonance p : q (with gcd(p,q) = 1 and q >= p) for the averaged model. Aptidal currently allows up to 9


#define typ double           //Renaming double as typ. To color gedit, update the field <context id="types" style-ref="type"> of the file c.lang used by gedit
                             //The idea is to allow quadruple precision in a future update.
                             //Typical paths where this file is located are /usr/share/gtksourceview-5/language-specs or /usr/share/libgedit-gtksourceview-300/language-specs


#if _3D_bool
#define Nd 6
#else
#define Nd 4
#endif

#endif
