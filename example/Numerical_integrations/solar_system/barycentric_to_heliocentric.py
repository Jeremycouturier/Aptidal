#Converting the barycentric initial condition S_0 of Laskar&Gastineau (2009) to both heliocentric canonical and heliocentric non canonical coordinates
#Pluto is excluded

import numpy as np


#Barycentric initial condition S_0 of Laskar&Gastineau (2009) in the ICRF at J2000 in AU and AU/days
masses = np.array([1.660120825459e-7, 2.447838287797e-6, 3.040432648963e-6, 3.227156082932e-7, 9.547919099414e-4, 2.858856700246e-4, 4.366249613222e-5, 5.151383772654e-5, 1.])
x  = np.array([-1.3723006138719493E-01,-7.2543875164004101E-01,-1.8429523869857883E-01,+1.3835794659467666E+00,+3.9940405805302190E+00,+6.3992728946824826E+00,+1.4424722672081669E+01,+1.6804916776997999E+01,-7.1364558632045867E-03])
y  = np.array([-4.0324074993791997E-01,-4.8921282357574918E-02,+8.8475982513797391E-01,-1.2458187406789204E-03,+2.7339319131097417E+00,+6.1720105115817834E+00,-1.2508913521035904E+01,-2.2982754337214441E+01,-2.6470343693797828E-03])
z  = np.array([-2.0141230251495065E-01,+2.3717655135959591E-02,+3.8381372854387652E-01,-3.7883156713073539E-02,+1.0745892889690563E+00,+2.2738481118833267E+00,-5.6826123000986799E+00,-9.8253491124291390E+00,-9.2298925054966601E-04])
vx = np.array([+2.1371774109527775E-02,+8.0349594196834088E-04,-1.7197730598869702E-02,+6.7687800264027093E-04,-4.5629353241739114E-03,-4.2869731379417009E-03,+2.6834835473886366E-03,+2.5846532680459892E-03,+5.3784605099164219E-06])
vy = np.array([-4.9330576077080499E-03,-1.8498595723320826E-02,-2.9096001847899629E-03,+1.3807279327919367E-02,+5.8747038486458630E-03,+3.5215868716827110E-03,+2.4552463507823063E-03,+1.6616665134371927E-03,-6.7581880103245664E-06])
vz = np.array([-4.8504663887848996E-03,-8.3727680942960545E-03,-1.2615424684121190E-03,+6.3148675769632316E-03,+2.6292698127070692E-03,+1.6388979154519450E-03,+1.0373752265128999E-03,+6.1578095014691201E-04,-3.0328531274874777E-06])


#The center of mass is not at zero because I removed Pluto. There is no point in integrating Pluto since dozens of other dwarf planets in the solar system have similar masses.
com_x  = 0.
com_y  = 0.
com_z  = 0.
com_vx = 0.
com_vy = 0.
com_vz = 0.
mass   = 0.
for i in range(9):
      com_x  = com_x  + masses[i]*x[i]
      com_y  = com_y  + masses[i]*y[i]
      com_z  = com_z  + masses[i]*z[i]
      com_vx = com_vx + masses[i]*vx[i]
      com_vy = com_vy + masses[i]*vy[i]
      com_vz = com_vz + masses[i]*vz[i]
      mass   = mass   + masses[i]

#Cancelling the center of mass
x  = x  - com_x/mass
y  = y  - com_y/mass
z  = z  - com_z/mass
vx = vx - com_vx/mass
vy = vy - com_vy/mass
vz = vz - com_vz/mass

#Heliocentric canonical coordinates. The positions are heliocentric but the speeds are barycentric.
x = x - x[-1]
y = y - y[-1]
z = z - z[-1]
print(f"x   = ({x[0]:.17}, {x[1]:.17}, {x[2]:.17}, {x[3]:.17}, {x[4]:.17}, {x[5]:.17}, {x[6]:.17}, {x[7]:.17})")
print(f"y   = ({y[0]:.17}, {y[1]:.17}, {y[2]:.17}, {y[3]:.17}, {y[4]:.17}, {y[5]:.17}, {y[6]:.17}, {y[7]:.17})")
print(f"z   = ({z[0]:.17}, {z[1]:.17}, {z[2]:.17}, {z[3]:.17}, {z[4]:.17}, {z[5]:.17}, {z[6]:.17}, {z[7]:.17})")
print(f"vx  = ({vx[0]:.17}, {vx[1]:.17}, {vx[2]:.17}, {vx[3]:.17}, {vx[4]:.17}, {vx[5]:.17}, {vx[6]:.17}, {vx[7]:.17})")
print(f"vy  = ({vy[0]:.17}, {vy[1]:.17}, {vy[2]:.17}, {vy[3]:.17}, {vy[4]:.17}, {vy[5]:.17}, {vy[6]:.17}, {vy[7]:.17})")
print(f"vz  = ({vz[0]:.17}, {vz[1]:.17}, {vz[2]:.17}, {vz[3]:.17}, {vz[4]:.17}, {vz[5]:.17}, {vz[6]:.17}, {vz[7]:.17})")
print()

#Heliocentric non-canonical coordinates. Both positions and speeds are heliocentric.
vx = vx - vx[-1]
vy = vy - vy[-1]
vz = vz - vz[-1]
print(f"x   = ({x[0]:.17}, {x[1]:.17}, {x[2]:.17}, {x[3]:.17}, {x[4]:.17}, {x[5]:.17}, {x[6]:.17}, {x[7]:.17})")
print(f"y   = ({y[0]:.17}, {y[1]:.17}, {y[2]:.17}, {y[3]:.17}, {y[4]:.17}, {y[5]:.17}, {y[6]:.17}, {y[7]:.17})")
print(f"z   = ({z[0]:.17}, {z[1]:.17}, {z[2]:.17}, {z[3]:.17}, {z[4]:.17}, {z[5]:.17}, {z[6]:.17}, {z[7]:.17})")
print(f"vx  = ({vx[0]:.17}, {vx[1]:.17}, {vx[2]:.17}, {vx[3]:.17}, {vx[4]:.17}, {vx[5]:.17}, {vx[6]:.17}, {vx[7]:.17})")
print(f"vy  = ({vy[0]:.17}, {vy[1]:.17}, {vy[2]:.17}, {vy[3]:.17}, {vy[4]:.17}, {vy[5]:.17}, {vy[6]:.17}, {vy[7]:.17})")
print(f"vz  = ({vz[0]:.17}, {vz[1]:.17}, {vz[2]:.17}, {vz[3]:.17}, {vz[4]:.17}, {vz[5]:.17}, {vz[6]:.17}, {vz[7]:.17})")
