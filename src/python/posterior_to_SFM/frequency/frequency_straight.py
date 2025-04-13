import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random

delt, Xmin, Xmax, Xint, Xext, Xhyp = np.loadtxt('../continuedSeparatrix.txt', dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5]))
delt = np.flip(delt)
Xmin = np.flip(Xmin)
Xmax = np.flip(Xmax)
Xint = np.flip(Xint)
Xext = np.flip(Xext)
Xhyp = np.flip(Xhyp)

freqs = np.loadtxt('/home/jeremy/Documents/GRSW/frequency_fixed.csv', dtype = np.float64, delimiter=',', unpack=True)

      
def levelLine2X(X, Y, delta):
      #Returns X1 and X2 such that (X1, 0) and (X2, 0) are on the same level line as (X, Y)
      
      H   = 1.5*delta*(X**2 + Y**2) - 0.25*(X**2 + Y**2)**2 + 2.*X
      D   = m.sqrt(X**2 + Y**2)
      sol = 0
      Sol = []
      #while (sol < 2):
      for mm in range(200):
            X0 = random.random()
            X0 = 4.*D*X0 - 2.*D
            h  = 1.
            it = 0
            cv = 1
            while(abs(h) > 1.e-11):
                  num   = X0**4 - 6.*delta*X0**2 - 8.*X0 + 4.*H
                  denom = 4.*X0**3 - 12.*delta*X0 - 8.
                  h     = -num/denom
                  X0    = X0 + h
                  it    = it + 1
                  if (it > 30):
                        h  = 1.e-200
                        cv = 0
            new_sol = cv
            for i in range(sol):
                  if (abs(X0 - Sol[i]) < 1.e-9):
                        new_sol = 0
            sol = sol + new_sol
            if (new_sol):
                  Sol.append(X0)
      if (sol == 0):
            raise Exception("Could not find at least one solution in function levelLine2X")
      elif (sol == 1):
            print("Warning : Could not find at least two solutions in function levelLine2X")
            return [Sol[0], Sol[0]]
      elif (sol == 4):
            D1 = abs(D - abs(Sol[0]))
            D2 = abs(D - abs(Sol[1]))
            D3 = abs(D - abs(Sol[2]))
            D4 = abs(D - abs(Sol[3]))
            DD = np.sort(np.array([D1, D2, D3, D4]))
            if (DD[0] == D1):
                  S1 = Sol[0]
            elif (DD[0] == D2):
                  S1 = Sol[1]
            elif (DD[0] == D3):
                  S1 = Sol[2]
            else:
                  S1 = Sol[3]
            if (DD[1] == D1):
                  S2 = Sol[0]
            elif (DD[1] == D2):
                  S2 = Sol[1]
            elif (DD[1] == D3):
                  S2 = Sol[2]
            else:
                  S2 = Sol[3]
            return [S1, S2]
      else:
            return [Sol[0], Sol[1]]

def topologie(delta):
      #Returns [Xmin, Xmax, Xint, Xext, Xhyp] as a function of delta
      #Instead of a direct calculation, extrapolates from file './continuedSeparatrix.txt'
            
      N       = len(delt[delt < delta])
      xminmin = Xmin[N - 1]
      xminmax = Xmin[N]
      xmaxmin = Xmax[N - 1]
      xmaxmax = Xmax[N]
      xintmin = Xint[N - 1]
      xintmax = Xint[N]
      xextmin = Xext[N - 1]
      xextmax = Xext[N]
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
      xmin = t*xminmin + (1. - t)*xminmax
      xmax = t*xmaxmin + (1. - t)*xmaxmax
      xint = t*xintmin + (1. - t)*xintmax
      xext = t*xextmin + (1. - t)*xextmax
      xhyp = t*xhypmin + (1. - t)*xhypmax
      return [xmin, xmax, xint, xext, xhyp]

def freq(X, delta):
      #Returns the frequency of the trajectory starting at (X, 0)
      #Instead of a direct calculation, extrapolates from file './frequency.txt'
      
      if (delta < -20. or delta >= 50.):
            raise Exception("Out of bounds.")
      if (X < -14. or X >= 14.):
            raise Exception("Out of bounds.")
            
      row = 35000.*(delta + 20.)/70.
      row = round(row, 9)
      row = int(row)
      col = 2800.*(X + 14.)/28.
      col = round(col, 9)
      col = int(col)
      
      Dmin = delta - np.mod(delta + 1.e-13, 0.002)
      Dmax = Dmin + 0.002
      xmin = X - np.mod(X + 1.e-13, 0.01)
      xmax = xmin + 0.01
      
      nu_min_min = freqs[col, row]
      nu_min_max = freqs[col + 1, row]
      nu_max_min = freqs[col, row + 1]
      nu_max_max = freqs[col + 1, row + 1]
      tD = (Dmax - delta)/0.002
      tX = (xmax - X)/0.01
      if (tD < -1.e-14 or tD > 1.+1.e-14 or tX < -1.e-14 or tX > 1.+1.e-14):
            print('Dmax =', Dmax)
            print('delta =', delta)
            print('xmax =', xmax)
            print('X =', X)
            raise Exception("Problem with tD or tX")
      
      nu_min = tD*nu_min_min + (1. - tD)*nu_max_min
      nu_max = tD*nu_min_max + (1. - tD)*nu_max_max
      nu1    = tX*nu_min + (1. - tX)*nu_max
      nu_min = tX*nu_min_min + (1. - tX)*nu_min_max
      nu_max = tX*nu_max_min + (1. - tX)*nu_max_max
      nu2    = tD*nu_min + (1. - tD)*nu_max #Should be equal to nu1
      if (abs(nu1 - nu2) > 1.e-9):
            raise Exception("Some weird sh*t going on")
      return 0.5*(nu1 + nu2)
      
def target(X, x, delta):

      #Computes the target, that is, the point (X_target, delta_target) on the corresponding family of elliptic fixed point with the same frequency
      #Returns [X_target, delta_target, distance] where distance is the distance between the target and the current point
      #(X, 0) and (x, 0) are on the same level line
      
      nu1 = freq(X, delta)
      nu2 = freq(x, delta)
      if (abs((nu1 - nu2)/nu1) > 1.e-2):
            print('Problem with frequencies: nu1 = ', nu1, 'nu2 = ', nu2)
      nu = 0.5*(nu1 + nu2)
      [xmin, xmax, xint, xext, xhyp] = topologie(delta)
      if (nu <= 3.**(0.5)*2.**(2./3.) and X >= xmin):
            distance = np.sqrt((X - 2.**(1./3.))**2 + delta**2)
            return [2.**(1./3.), 0., distance]
      if (X >= xmin):
            if ((X <= xmax and delta >= 1.) or (X <= xmax and delta < 1. and delta > 0. and x >= 0.)): #Corresponding branch is the right part of the internal branch
                  [X_target, delta_target] = dicho('internal', 0., 49.999, nu)
            else: #Corresponding branch is the left part of the internal branch
                  [X_target, delta_target] = dicho('internal', -20., 0., nu)
      else: #Corresponding branch is the external circulation
            [X_target, delta_target] = [0., 0.]
      distance = np.sqrt((X - X_target)**2 + (delta - delta_target)**2)
      return [X_target, delta_target, distance]

def dicho(branch, min_delta, max_delta, frequency):
      #Returns the [X, delta] value on the corresponding branch having frequency equal to frequency
      #branch is either 'external', 'internal' or 'max'
      
      eps      = 2.5e-5
      
      if (branch == 'internal'):
            X_left  = topologie(min_delta)[2]
            X_right = topologie(max_delta)[2]
      elif (branch == 'external'):
            X_left  = topologie(min_delta)[3]
            X_right = topologie(max_delta)[3]
      else:
            X_left  = topologie(min_delta)[1]
            X_right = topologie(max_delta)[1]
      nu_left  = freq(X_left,  min_delta)
      nu_right = freq(X_right, max_delta)
      if (frequency > max(nu_left, nu_right) or frequency < min(nu_left, nu_right)):
            print('nu_left =', nu_left, 'nu_right =', nu_right, 'frequency = ', frequency)
            raise Exception("Problem with dichotomy")
      
      while(abs(max_delta - min_delta) > eps):
            middle_delta = 0.5*(min_delta + max_delta)
            if (branch == 'internal'):
                  middle_X = topologie(middle_delta)[2]
            elif (branch == 'external'):
                  middle_X = topologie(middle_delta)[3]
            else:
                  middle_X = topologie(middle_delta)[1]
            middle_freq = freq(middle_X, middle_delta)
            if ((frequency < nu_right and frequency > middle_freq) or (frequency > nu_right and frequency < middle_freq)):
                  X_left    = middle_X
                  nu_left   = middle_freq
                  min_delta = middle_delta
            else:
                  X_right   = middle_X
                  nu_right  = middle_freq
                  max_delta = middle_delta
      return [0.5*(X_left + X_right), 0.5*(min_delta + max_delta)]


XX = np.linspace(-10., 30., 4001)
YY = np.linspace(0., 10., 1001)
ZZ = np.zeros((1001, 4001))
hits = np.zeros((1001, 4001))

dlt = np.arange(-5., 15., 0.01)
Xs  = np.arange(0., 6., 0.01)

print(freq(-6.5, -5.))
print(freq(-6.5, 20.))
print(freq(6.5, -5.))
print(freq(6.5, 20.))

[x1, x2] = levelLine2X(0., 0., -5.)
print(x1, x2)
[x1, x2] = levelLine2X(0.1, 0., -5.)
print(x1, x2)
'''
for i in range(600):
      for j in range(2000):
            delta    = dlt[j]
            X        = Xs[i]
            [x1, x2] = levelLine2X(X, 0., delta)
            X        = max(x1, x2)
            x        = min(x1, x2)
            if (delta <= 50. and delta >= -20. and X <= 14. and x > -14.):
                  nu1      = freq(X, delta)
                  nu2      = freq(x, delta)
                  if (abs((nu1 - nu2)/nu1) > 1.e-2):
                        print('Problem with frequencies: nu1 = ', nu1, 'nu2 = ', nu2)
                  nu = 0.5*(nu1 + nu2)
                  [xmin, xmax, xint, xext, xhyp] = topologie(delta)
                  [X_target, delta_target, distance] = target(X, x, delta)
                  Normalized_delta = delta_target
                  if (nu <= 3.**(0.5)*2.**(2./3.) and X >= xmin):
                        S = xmax
                        D = np.sqrt((S - 2.**(1./3.))**2 + delta**2)
                        NormalizedX = distance/D
                  elif (X >= xmin):
                        if ((X <= xmax and delta >= 1.) or (X <= xmax and delta < 1. and delta > 0. and x >= 0.)):
                              S = xmax
                              D = np.sqrt((S - 2.**(1./3.))**2 + delta**2)
                              NormalizedX = distance/D
                        else:
                              [S, dl] = dicho('max', -20., 1., nu)
                              D = np.sqrt((S - X_target)**2 + (dl - delta_target)**2)
                              NormalizedX = distance/D
                  else:
                        NormalizedX = -1.
                  row = 4000.*(Normalized_delta + 10.)/40.
                  row = round(row, 9)
                  row = int(row)
                  col = 1000.*(NormalizedX + 0.)/10.
                  col = round(col, 9)
                  col = int(col)
                  if (row >= 0 and col >= 0 and row <= 4000 and col <= 1000):
                        hit = hits[col, row]
                        previous_value = ZZ[col, row]
                        new_value      = nu
                        ZZ[col, row]   = (hit * previous_value + new_value)/(hit + 1.)
                        hits[col, row] = hits[col, row] + 1.
      print(i+1,'/',600)

out_path = './frequency_straight.csv'
out_file = open(out_path, 'w')

N = ZZ.shape[1]

for i in range(N):
      row = ZZ[:,i]
      n   = len(row)
      for j in range(n - 1):
            out_file.write(str(row[j]) + ",")
      out_file.write(str(row[-1]) + "\n")
      print(i+1,'/',N)
      
out_file.close()
'''
      






