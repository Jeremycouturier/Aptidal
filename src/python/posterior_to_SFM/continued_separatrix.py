import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random

out_path = "/home/jeremy/Documents/CleanPosteriors/continuedSeparatrix.txt"
out_file = open(out_path, "w")

def area(delta, Xmax):

      H  = 1.5*delta*Xmax**2 - 0.25*Xmax**4 + 2.*Xmax
      A  = 0.
      S  = 1
      X  = Xmax - 1.e-14
      dX = 4.026416e-5
      Xmin_was_set = 0
      while(S):
            D = 8.*X + 9.*delta**2 - 4.*H
            x = X - dX
            d = 8.*x + 9.*delta**2 - 4.*H
            if (D < 0.):
                  if (Xmin_was_set):
                        return [A, Xmin]
                  else:
                        return [A, X]
            if (d < 0.):
                  d = 0.
            Y2lgr = 3.*delta - X**2 + m.sqrt(D)
            Y2sml = 3.*delta - X**2 - m.sqrt(D)
            y2lgr = 3.*delta - x**2 + m.sqrt(d)
            y2sml = 3.*delta - x**2 - m.sqrt(d)
            if (Y2lgr < 0.):
                  S = 0
            else:
                  Ylgr = m.sqrt(Y2lgr)
                  if (Y2sml < 0. and y2lgr >= 0. and y2sml < 0.):
                        Ysml = 0.
                        ylgr = m.sqrt(y2lgr)
                        ysml = 0.
                  elif (Y2sml < 0. and y2lgr >= 0. and y2sml >= 0.):
                        Ysml = 0.
                        ylgr = m.sqrt(y2lgr)
                        ysml = m.sqrt(y2sml)
                        Xmin = X
                        Xmin_was_set = 1
                  elif (Y2sml >= 0. and y2lgr >= 0. and y2sml >= 0.):
                        Ysml = m.sqrt(Y2sml)
                        ylgr = m.sqrt(y2lgr)
                        ysml = m.sqrt(y2sml)
                  elif (Y2sml >= 0. and y2lgr < 0. and y2sml < 0.):
                        Ysml = m.sqrt(Y2sml)
                        ylgr = 0.
                        ysml = 0.
                  elif (Y2sml < 0. and y2lgr < 0. and y2sml < 0.):
                        Ysml = 0.
                        ylgr = 0.
                        ysml = 0.
                        Xmin = X
                        Xmin_was_set = 1
                  else:
                        print("Error : Not in one of expected case\n")
                        print("X      = ", X)
                        print("X - dX = ", x)
                        print("Y^2(X)      = ", Y2lgr, Y2sml)
                        print("Y^2(X - dX) = ", y2lgr, y2sml)
                  A = A + dX*(Ylgr + ylgr - Ysml - ysml)
            X = x
      if (Xmin_was_set):
            return [A, Xmin]
      else:
            return [A, X]


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


Obj = 6.*m.pi #Area of the cardioid

d_delta   = 2.e-3
delta     = 50.
while (delta > 1.00001):
      [X1, X2, X3, X4, X5] = topologie(delta)
      out_file.write(str(round(delta,3)) + " " + str(round(X4,5)) + " " + str(round(X5,5)))
      out_file.write('\n')
      delta = delta - d_delta

previousX = 2.488670410156049#3.
out_file.write(str(1.0) + " " + str(-1.0) + " " + str(3.0))
out_file.write('\n')
delta = -15#0.998

while (delta >= -20.0001):
      LargeTrial = previousX
      SmallTrial = previousX - 0.005
      print("delta = ", delta)
      LargeArea  = area(delta, LargeTrial)[0]
      SmallArea  = area(delta, SmallTrial)[0]
      if (LargeArea < Obj or SmallArea > Obj):
            raise Exception("Error in dichotomy\n")
      error = 1.
      
      while (error > 9.e-6):
            MiddleTrial = 0.5*(SmallTrial + LargeTrial)
            MiddleArea  = area(delta, MiddleTrial)[0]
            if (MiddleArea > Obj):
                  LargeTrial = MiddleTrial
            else:
                  SmallTrial = MiddleTrial
            error = LargeTrial - SmallTrial

      Xmax = 0.5*(SmallTrial + LargeTrial)
      [Area, Xmin] = area(delta, Xmax)
      #Area = area(delta, Xmax)
      if (abs(Area - Obj)/Obj > 1.e-3):
            print("Area = ", Area)
            raise Exception("Area error is too large\n")
      print("Xmax, Xmin = ", Xmax, Xmin)
      out_file.write(str(round(delta,3)) + " " + str(round(Xmin,5)) + " " + str(round(Xmax,5)))
      out_file.write('\n')
      previousX = Xmax
      delta     = delta - d_delta

out_file.close()










