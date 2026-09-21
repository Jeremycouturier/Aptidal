/*
               AAA                PPPPPPPPPPPPPPPPP    TTTTTTTTTTTTTTTTTTTTTTT IIIIIIIIII DDDDDDDDDDDDD                   AAA                LLLLLLLLLLL             
              A:::A               P::::::::::::::::P   T:::::::::::::::::::::T I::::::::I D::::::::::::DDD               A:::A               L:::::::::L             
             A:::::A              P::::::PPPPPP:::::P  T:::::::::::::::::::::T I::::::::I D:::::::::::::::DD            A:::::A              L:::::::::L             
            A:::::::A             PP:::::P     P:::::P T:::::TT:::::::TT:::::T II::::::II DDD:::::DDDDD:::::D          A:::::::A             LL:::::::LL             
           A:::::::::A              P::::P     P:::::P TTTTTT  T:::::T  TTTTTT   I::::I     D:::::D    D:::::D        A:::::::::A              L:::::L               
          A:::::A:::::A             P::::P     P:::::P         T:::::T           I::::I     D:::::D     D:::::D      A:::::A:::::A             L:::::L               
         A:::::A A:::::A            P::::PPPPPP:::::P          T:::::T           I::::I     D:::::D     D:::::D     A:::::A A:::::A            L:::::L               
        A:::::A   A:::::A           P:::::::::::::PP           T:::::T           I::::I     D:::::D     D:::::D    A:::::A   A:::::A           L:::::L               
       A:::::A     A:::::A          P::::PPPPPPPPP             T:::::T           I::::I     D:::::D     D:::::D   A:::::A     A:::::A          L:::::L               
      A:::::AAAAAAAAA:::::A         P::::P                     T:::::T           I::::I     D:::::D     D:::::D  A:::::AAAAAAAAA:::::A         L:::::L               
     A:::::::::::::::::::::A        P::::P                     T:::::T           I::::I     D:::::D     D:::::D A:::::::::::::::::::::A        L:::::L               
    A:::::AAAAAAAAAAAAA:::::A       P::::P                     T:::::T           I::::I     D:::::D    D:::::D A:::::AAAAAAAAAAAAA:::::A       L:::::L         LLLLLL
   A:::::A             A:::::A    PP::::::PP                 TT:::::::TT       II::::::II DDD:::::DDDDD:::::D A:::::A             A:::::A    LL:::::::LLLLLLLLL:::::L
  A:::::A               A:::::A   P::::::::P                 T:::::::::T       I::::::::I D:::::::::::::::DD A:::::A               A:::::A   L::::::::::::::::::::::L
 A:::::A                 A:::::A  P::::::::P                 T:::::::::T       I::::::::I D::::::::::::DDD  A:::::A                 A:::::A  L::::::::::::::::::::::L
AAAAAAA                   AAAAAAA PPPPPPPPPP                 TTTTTTTTTTT       IIIIIIIIII DDDDDDDDDDDDD    AAAAAAA                   AAAAAAA LLLLLLLLLLLLLLLLLLLLLLLL
*/

/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    parameters.h                                                ********/
/******** @brief   Main header file. The user who does not modify Aptidal      ********/
/********          shall only modify this file                                 ********/
/******** @author  Jérémy COUTURIER <jeremycouturier.com>                      ********/
/********                                                                      ********/
/******** @section LICENSE                                                     ********/
/******** Copyright (c) 2026 Jérémy COUTURIER                                  ********/
/********                                                                      ********/
/******** Aptidal is free software. You can redistribute it and/or modify      ********/
/******** it under the terms of the GNU General Public License as published by ********/
/******** the Free Software Foundation, either version 3 of the License, or    ********/
/******** (at your option) any later version.                                  ********/
/********                                                                      ********/
/******** Aptidal is distributed in the hope that it will be useful,           ********/
/******** but WITHOUT ANY WARRANTY; without even the implied warranty of       ********/
/******** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ********/
/******** GNU General Public License for more details.                         ********/
/********                                                                      ********/
/******** You should have received a copy of the GNU General Public License    ********/
/******** along with Aptidal. If not, see <http://www.gnu.org/licenses/>.      ********/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/*********************************************************************/
/******** Defining the output path. Must end with / and exist ********/
/*********************************************************************/
#define pth "/path/towards/output/location/with/slash/at/the/end/"



/********************************************************************************************/
/******** Defining the resonance chain p_1 : p_2 : ... : p_n. The p_i are integers   ********/
/******** given in increasing order whose ratio indicate the orbital periods ratios  ********/
/******** Two consecutive p_i can be equal (coorbital planets) and p_i = 0 indicates ********/
/******** that planet n° i is not in resonance with other planets. When working with ********/
/******** the complete Hamiltonian, this only affects the output coordinates.        ********/
/********************************************************************************************/
#define resonance_chain {8, 12, 18, 27}



/****************************************/
/******** Defining some booleans ********/
/****************************************/

/******** Boolean relative to numerical simulations ********/
#define tides_bool 1         //Determines if there are tides raised by the star on the planets in the system. See Couturier et al. (2026), Appendix B.2.3.
#define GR_bool 0            //Determines if the first PPN order of General Relativity is included. See Saha&Tremaine (1994), Sect. 5. Implemented but not yet verified. Could contain errors.
#define _3D_bool 0           //Determines if the problem is 3D or coplanar. 3D is only possible with the complete Hamiltonian
#define toInvar_bool 0       //Determines if the system is rotated so that the total angular momentum points in the z-direction. Irrelevant if _3D_bool is 0. Planet's rotations are not rotated.
#define close_enc_bool 0     //Determines if close encounters are checked when integrating the complete system. Does not handle large timesteps.

/******** Boolean relative to coordinate systems ********/
#define canon_input_bool 1   //Determines if the heliocentric coordinates in input  are canonical (heliocentric position, barycentric speed) or not (heliocentric position, heliocentric speed)
#define canon_output_bool 1  //Determines if the heliocentric coordinates in output are canonical (heliocentric position, barycentric speed) or not (heliocentric position, heliocentric speed)
#define ellip_input_bool 1   //Determines if the input  coordinates are elliptic or cartesian. If cartesian (0), the center of mass must have zero speed.
#define ellip_output_bool 1  //Determines if the output coordinates are elliptic or cartesian.

/******** Boolean relative to the averaged model ********/
#define non_resonant_bool 1  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are taken into account in the averaged Hamiltonian
#define one_more_deg_bool 0  //Determines if non-resonant terms (those associated with the inequality (ki,kj)=(0,0)) are pushed one degree further in eccentricity than the resonant terms.

                            
/**********************************************************************************************************************************/
/****************************************  Defining some physical constants *******************************************************/
/**********************************************************************************************************************************/
/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** planet and the orbital period of a massless particle with semi-major axis the unit of length.                    ********/
/**********************************************************************************************************************************/
#define how_many_planet 4              //Number of planets in the system. Must match the length of the resonance chain. Minimum is 2
#define G 39.478417604357434475337964  //Gravitational constant. Should be left to 4*pi^2 (the unit of time is the orbital period of a massless particle with semi-major axis the unit of length)
#define m0 1.                          //Mass of the star. Should be left to 1. (the unit of mass is the star mass)
#if GR_bool
#define c_light 173.14463267424        //Speed of light in units of length per units of time. Only when GR_bool = 1. Here in AU/days
#endif
#define body_masses {1.103e-5,            1.195e-5,             2.422e-5,             1.737e-5}             //Masses of the planets in stellar masses
#define body_sma    {1.,                  1.31063388956530,     1.71785897200592,     2.25140923008155}     //Initial and nominal semi-major axes. vx coordinate if ellip_input_bool is 0
#define body_ecc    {.025969963962830640, .059382364145342155,  .037409363918984279,  .029034289643334182}  //Initial eccentricities. vy coordinate if ellip_input_bool is 0
#define body_lambda {2.65725525333552,    3.78233356367910,     3.48546755071036,     .14584719428602}      //Initial mean longitudes. x coordinate if ellip_input_bool is 0
#define body_varpi {-0.25011308400757,   -3.39240540126691,    -.25098431811849,     -3.39237472136929}     //Initial longitudes of the periapses. y coordinate if ellip_input_bool is 0

#if _3D_bool
#define body_inc    {}  //Initial inclinations of the planets in radians. vy coordinate if ellip_input_bool is 0
#define body_Omeg   {}  //Initial longitudes of the ascending node of the planets in radians. y coordinate if ellip_input_bool is 0
#endif

#if tides_bool
#define body_radii  {.0016682,            .0018038,             .0036077,             .0026040}             //Radii of the planets, in units of length
#define body_k2     {1.5,                 1.5,                  1.5,                  1.5}                  //Second Love number of the planets
#define body_Dt     {.06578400722471152,  .06845556018452882,   .08839218851765285,   .17204735546952463}   //Tidal timelag of the bodies, in units of time. k2/Q = {6.2e-1, 4.3e-1, 3.7e-1, 4.8e-1}
#define body_alpha  {.33,                 .33,                  .33,                  .33}                  //Dimensionless structure constant of the bodies. 2/5 for an homogeneous body.
#if _3D_bool
#define body_Omegx  {}  //x-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegy  {}  //y-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#define body_Omegz  {}  //z-coordinate of the initial sideral rotation of the bodies, in radians/unit of time
#else
#define body_Omega  {6.283219958851004,   4.187553542589783,    2.790640330614937,    1.8599530293662827}   //Initial sideral rotation of the bodies, in radians/unit of time
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
