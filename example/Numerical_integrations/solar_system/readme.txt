In this example, we integrate the solar system composed of the Sun and eight planets.
We use a SABAH1064 integrator and initial conditions provided in Appendix of Laskar&Gastineau (2009)

You can use either of the parameters files (you have to rename them to parameters.h).
One of them uses canonical heliocentric coordinates whereas the other uses non-canonical heliocentric coordinates.
The booleans are already properly specified in each of them, and the output path is the only parameter that you need to change

The Earth is replaced by the Earth-Moon barycenter. Likewise, planets with satellites are replaced by the barycenter of the system planet+satellites.
We use a system of units where the unit of time is the day, the unit of length is the AU and the unit of mass is the Solar mass.
