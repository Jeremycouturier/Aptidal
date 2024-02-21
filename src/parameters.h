#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_


#define typ double      // Renaming double as typ so it can be easily changed to long double, for example.
                        // Open the file /usr/share/gtksourceview/language-specs/c.lang and then find the field
                        // <context id="types" style-ref="type">. Add the line <keyword>typ</keyword> to it
                        // and color gedit with C. The keyword typ should now be colored after a reboot
#define absolute fabs   // The function returning the absolute value of a typ. Choose fabs (resp. fabsf or fabsl) if typ is double (resp. float or long double)
#define integral floor  // The function returning the floor of a typ. Choose floor (resp. floorf or floorl) if typ is double (resp. float or long double)



/**************************************************/
/******** Defining some physical constants ********/
/**************************************************/

/******** By default, the units of mass, length and time are the star's mass, the initial semi-major axis of the innermost ********/
/******** orbiting body and the initial orbital period of the innermost orbiting body (hence G=4*pi^2). The canonical      ********/
/******** heliocentric coordinates are used (heliocentric positions and barycentric velocities). See Laskar & Robutel 1995 ********/
#define how_many_planet 5                                                           //Number of orbiting bodies in the system
#define G 39.47841760435743                                                         //Gravitational constant is set to 4*pi^2, so the orbital period of the closest orbiting body is 1
#define body_masses {1.0 , 0.00238,  0.000488, 0.000793, 0.00132,  0.00056}         //Masses of the bodies of the system, beginning with the star.
#define body_sma          {1.0,      1.5874,   2.071,    2.08,     2.78}            //Initial semi-major axes             of the orbiting bodies.
#define body_ecc          {0.0,      0.0,      0.0,      0.0,      0.0}             //Initial eccentricities              of the orbiting bodies.
#define body_lambda       {0.0,      0.0,      0.0,      0.0,      0.0}             //Initial mean longitudes             of the orbiting bodies.
#define body_varpi        {0.0,      0.0,      0.0,      0.0,      0.0}             //Initial longitude of the pericenter of the orbiting bodies.
#define body_radii  {0.12, 0.006744, 0.010036, 0.008982, 0.00733,  0.00438}         //Radii               of the bodies of the system, beginning with the star. Unimportant if no tides
#define body_k2     {0.0,  0.8,      1.32,     0.87,     2.77,     2.33}            //Second Love numbers of the bodies of the system, beginning with the star. Unimportant if no tides
#define body_Q  {10000.0,  234.8,    78.9,     1633.0,   687.9,    1239.8}          //Quality factors     of the bodies of the system, beginning with the star. Unimportant if no tides
//#define resonance_chain   {2,        5,        5,        6,        0}               //The resonance chain p_1:p_2:p_3: ... :p_n of the system.
#define resonance_chain   {1, 2, 3, 3, 5}



/******************************************/
/******** Defining some thresholds ********/
/******************************************/
#define max_deg 3                                                                   //Maximum degree in eccentricity in the development of the Hamiltonien. Aptidal allows up to 3
#define max_res 9                                                                   //Maximum value of q for a resonance p:q (with gcd(p,q) = 1). Aptidal currently allows up to 9


/******** Some strings ********/
#define string1 "sq2DsL * cos(p*ll - (p+1)*ll' + vp)"
#define string2 "sq2DsL' * cos(p*ll - (p+1)*ll' + vp')"
#define string3 "sq2DsL**2 * cos(2*p*ll - 2*(p+1)*ll' + 2*vp)"
#define string4 "sq2DsL*sq2DsL' * cos(2*p*ll - 2*(p+1)*ll' + vp + vp')"
#define string5 "sq2DsL'**2 * cos(2*p*ll - 2*(p+1)*ll' + 2*vp')"
#define string6 "sq2DsL**3 * cos(p*ll - (p+1)*ll' + vp)"
#define string7 "sq2DsL**2*sq2DsL' * cos(p*ll - (p+1)*ll' + vp')"
#define string8 "sq2DsL**2*sq2DsL' * cos(p*ll - (p+1)*ll' + 2*vp - vp')"
#define string9 "sq2DsL*sq2DsL'**2 * cos(p*ll - (p+1)*ll' + vp)"
#define string10 "sq2DsL*sq2DsL'**2 * cos(p*ll - (p+1)*ll' + 2*vp' - vp)"
#define string11 "sq2DsL'**3 * cos(p*ll - (p+1)*ll' + vp')"
#define string12 "sq2DsL**3 * cos(3*p*ll - 3*(p+1)*ll' + 3*vp)"
#define string13 "sq2DsL**2*sq2DsL' * cos(3*p*ll - 3*(p+1)*ll' + 2*vp + vp')"
#define string14 "sq2DsL*sq2DsL'**2 * cos(3*p*ll - 3*(p+1)*ll' + 2*vp' + vp)"
#define string15 "sq2DsL'**3 * cos(3*p*ll - 3*(p+1)*ll' + 3*vp')"




















#endif
