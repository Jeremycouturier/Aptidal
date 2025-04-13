import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random

in_path  = "./continuedSeparatrix.txt"
out_path = "./continuedSeparatrix5col.txt"
out_file = open(out_path, "w")

delt, Xmin, Xmax = np.loadtxt(in_path, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2]))


def topologie(delta):
      #If delta <= 1, returns the X-coordinate of the unique (elliptic) fixed point.
      #Else, returns [X1, X2, X3, X4, X5] where X1 is the resonant elliptic fixed point,
      #X2 is the external circulation elliptic fixed point, and X3 to X5 are
      #the three X-coordinates of the separatrix by increasing order
      if (delta <= 1.):
            u = (1. + m.sqrt(1. - delta**3))**(1./3.)
            v =  1. - m.sqrt(1. - delta**3)
            if (v < 0.):
                  v = -((-v)**(1./3.))
            else:
                  v = v**(1./3.)
            return [u + v]
      else:
            u3     = 1. + cm.sqrt(1. - delta**3)
            v3     = 1. - cm.sqrt(1. - delta**3)
            j      = cm.exp( 2.0*1j*cm.pi/3.0)
            jb     = cm.exp(-2.0*1j*cm.pi/3.0)
            mod_u3 = cm.polar(u3)[0]
            arg_u3 = cm.polar(u3)[1]
            mod_v3 = cm.polar(v3)[0]
            arg_v3 = cm.polar(v3)[1]
            u      = mod_u3**(1.0/3.0)*cm.exp(1j*arg_u3/3.0)
            v      = mod_v3**(1.0/3.0)*cm.exp(1j*arg_v3/3.0)
            X1 = np.real(u + v)
            X2 = np.real(j*u  + jb*v)
            X3 = np.real(jb*u + j*v)
            if (X2**2 > delta and X2**2 < 3.*delta):
                  X_hyp = X2
                  X_ell = X3
            else:
                  X_hyp = X3
                  X_ell = X2
            if (delta < 1.0001):
                  return [X1, X_ell, X_hyp, X_ell + 0.5*(X_ell - X_hyp), 3.]
            H   = 1.5*delta*X_hyp**2 - 0.25*X_hyp**4 + 2.*X_hyp
            sol = 0
            Sol = []
            while (sol < 3):
                  X0 = random.random()
                  X0 = 4.*X_hyp*X0 - 2.*X_hyp
                  h  = 1.
                  it = 0
                  while(abs(h) > 1.e-11 and it < 1000):
                        num   = X0**4 - 6.*delta*X0**2 - 8.*X0 + 4.*H
                        denom = 4.*X0**3 - 12.*delta*X0 - 8.
                        h     = -num/denom
                        X0    = X0 + h
                        it    = it + 1
                  new_sol = 1
                  if (it == 1000):
                        new_sol = 0
                  for i in range(sol):
                        if (abs(X0 - Sol[i]) < 1.e-6):
                              new_sol = 0
                  sol = sol + new_sol
                  if (new_sol):
                        Sol.append(X0)
            Sol.sort()
            return [X1, X_ell, Sol[0], Sol[1], Sol[2]]


for i in range(len(delt)):
      delta = delt[i]
      xmin  = Xmin[i]
      xmax  = Xmax[i]
      if (delta > 1.):
            [X1, X2, X3, X4, X5] = topologie(delta)
            out_file.write(str(round(delta,3)) + " " + str(round(X4,5)) + " " + str(round(X5,5)) + " " + str(round(X1,5)) + " " + str(round(X2,5)) + " " + str(round(X3,5)))
            out_file.write('\n')
      elif (delta < 1.):
            [X1] = topologie(delta)
            out_file.write(str(round(delta,3)) + " " + str(round(xmin,5)) + " " + str(round(xmax,5)) + " " + str(round(X1,5)) + " " + str(0.0) + " " + str(0.0))
            out_file.write('\n')
      else:
            out_file.write(str(round(delta,3)) + " " + str(round(xmin,5)) + " " + str(round(xmax,5)) + " " + str(2.0) + " " + str(-1.0) + " " + str(-1.0))
            out_file.write('\n')


out_file.close()










