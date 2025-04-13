import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random

in_path  = "../continuedSeparatrix.txt"

delt, Xmin, Xmax, Xint, Xext, Xhyp = np.loadtxt(in_path, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5]))
delt = np.flip(delt)
Xmin = np.flip(Xmin)
Xmax = np.flip(Xmax)
Xint = np.flip(Xint)
Xext = np.flip(Xext)
Xhyp = np.flip(Xhyp)

def rk4(X, delta, estT):
      #Starting from (X,0), integrates with a Runge-Kutta 4 until Y changes sign (half a period)
      #Returns the corresponding time
      #estT is a rough order of magnitude of the period

      X0    = X
      Y0    = 0.
      count = 0
      dt    = estT/192.
      go_on = 1 
      
      while(go_on):
            Zx    = X0
            Zy    = Y0
            k1_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k1_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            Zx    = X0 + 0.5*dt*k1_x
            Zy    = Y0 + 0.5*dt*k1_y
            k2_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k2_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            Zx    = X0 + 0.5*dt*k2_x
            Zy    = Y0 + 0.5*dt*k2_y
            k3_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k3_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            Zx    = X0 + dt*k3_x
            Zy    = Y0 + dt*k3_y
            k4_x  = -3.*delta*Zy + Zy*(Zx**2 + Zy**2)
            k4_y  =  3.*delta*Zx - Zx*(Zx**2 + Zy**2) + 2.
            oldY  = Y0
            X0    = X0 + dt/6.*(k1_x + 2.*k2_x + 2.*k3_x + k4_x)
            Y0    = Y0 + dt/6.*(k1_y + 2.*k2_y + 2.*k3_y + k4_y)
            count = count + 1
            if (Y0 == 0. or oldY*Y0 < 0.):
                  go_on = 0
                  time  = count*dt - dt*abs(Y0/(Y0 - oldY))
      return time

def freq(X, delta, estT):
      #Returns the frequency of the orbit with initial condition (X, 0)
      
      # Geting the topology
      N       = len(delt[delt < delta])
      xintmin = Xint[N - 1]
      xintmax = Xint[N]
      xextmin = Xext[N - 1]
      xextmax = Xext[N]
      xminmin = Xmin[N - 1]
      xminmax = Xmin[N]
      xmaxmin = Xmax[N - 1]
      xmaxmax = Xmax[N]
      xhypmin = Xhyp[N - 1]
      xhypmax = Xhyp[N]
      Dmin    = delt[N - 1]
      Dmax    = delt[N]
      if (delta > Dmax or delta < Dmin):
            print("N = ", N)
            print("delta = ", delta)
            print("Dmin = ", Dmin)
            print("Dmax = ", Dmax)
            raise Exception("Problem with Dmin or Dmax\n")
      t    = (Dmax - delta)/(Dmax - Dmin)
      xint = t*xintmin + (1. - t)*xintmax
      xext = t*xextmin + (1. - t)*xextmax
      xmin = t*xminmin + (1. - t)*xminmax
      xmax = t*xmaxmin + (1. - t)*xmaxmax
      xhyp = t*xhypmin + (1. - t)*xhypmax
      
      #Case disjonction
      if  (abs(X - xint) < 0.0005 or abs(X - xext) < 0.0005): #Quadratic approximation. Frequency obtained from eigenvalues
            Xbar = xint
            return np.sqrt((Xbar**2 - 3.*delta)*(3.*Xbar**2 - 3.*delta))
      elif(abs(X - xext) < 0.0005 and delta >= 1.):
            Xbar = xext
            return np.sqrt((Xbar**2 - 3.*delta)*(3.*Xbar**2 - 3.*delta))
      elif(delta > 0.9999 and (abs(X - xmin) < 0.0001 or abs(X - xmax) < 0.0001 or abs(X - xhyp) < 0.0001)): #Too close from the separatrix for confort
            return 0.
            
      else: #Far from separatrix and fixed points. Computing frequency from integral
            P    = 2.*rk4(X, delta, estT)
            return 2.*np.pi/P


ds = np.linspace(48., 50., 1001, endpoint = True)


out_path = '../frequency_48_50.txt'
out_file = open(out_path, "w")

for d in ds:

      # Retrieving X at the internal resonance
      N       = len(delt[delt < d])
      xintmin = Xint[N - 1]
      xintmax = Xint[N]
      Dmin    = delt[N - 1]
      Dmax    = delt[N]
      if (d > Dmax or d < Dmin):
            print("N = ", N)
            print("delta = ", d)
            print("Dmin = ", Dmin)
            print("Dmax = ", Dmax)
            raise Exception("Problem with Dmin or Dmax\n")
      t    = (Dmax - d)/(Dmax - Dmin)
      xint = t*xintmin + (1. - t)*xintmax
      
      #Computing frequency
      nu = freq(xint, d, 1.)
      P  = 2.*np.pi/nu
      up_freq   = []
      down_freq = []
      xint  = xint - np.mod(xint, 0.01)
      xup   = xint
      xdown = xint - 0.01
      Pup   = P
      Pdown = P
      while(xup < 14.0000001):
            nu  = freq(xup, d, Pup)
            if (nu != 0.):
                  Pup = 2.*np.pi/nu
            else:
                  Pup = 0.1
            up_freq.append(nu)
            xup = xup + 0.01
      while(xdown > -14.0000001):
            nu    = freq(xdown, d, Pdown)
            if (nu != 0.):
                  Pdown = 2.*np.pi/nu
            else:
                  Pdown = 0.1
            down_freq.append(nu)
            xdown = xdown - 0.01
      
      up_freq   = np.array(up_freq)
      down_freq = np.array(down_freq)
      freqs     = np.concatenate((np.flip(down_freq), up_freq))
      N         = len(freqs)
      for i in range(N - 1):
            out_file.write(str(round(freqs[i],5)) + " ")
      out_file.write(str(round(freqs[-1],5)) + "\n")
      print('delta = ', d)

out_file.close()










