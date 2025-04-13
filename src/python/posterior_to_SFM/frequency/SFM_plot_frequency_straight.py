import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random

error_bar     = 0
radLeleu      = "/home/jeremy/Documents/CleanPosteriors/Leleu/"
radHadden_hm  = "/home/jeremy/Documents/CleanPosteriors/Hadden/high_mass/"
radHadden_df  = "/home/jeremy/Documents/CleanPosteriors/Hadden/default/"

N_pairs       = 5
pair_in_chain = [[0, 1], [0, 1], [0, 1], [0, 1], [0, 1]]
dummy_fst_col = [1, 1, 1, 0, 0] #Set to 1 if the first column of the posterior file is random
paths         = ["Kepler49bc", "Kepler85bc", "Kepler54bc", "Kepler1705bc", "Kepler1972bc"]
rad           = [radLeleu, radLeleu, radLeleu, radLeleu, radLeleu]
resonances    = [2, 2, 2, 4, 2] #Value of p such that the resonance is p : p + 1
plot_color    = ['purple', 'gold', 'blue', 'red', 'green']
lbl           = ['Kepler-49 (2:3)', 'Kepler-85 (2:3)', 'Kepler-54 (2:3)', 'Kepler-1705 (4:5)', 'Kepler-1972 (2:3)']
marker        = ['o', 'o', 'o', 'o', 'o']

delt, Xmin, Xmax, Xint, Xext, Xhyp = np.loadtxt('../continuedSeparatrix.txt', dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5]))
delt = np.flip(delt)
Xmin = np.flip(Xmin)
Xmax = np.flip(Xmax)
Xint = np.flip(Xint)
Xext = np.flip(Xext)
Xhyp = np.flip(Xhyp)

freqs = np.loadtxt('/home/jeremy/Documents/GRSW/frequency_fixed.csv', dtype = np.float64, delimiter=',', unpack=True)

def ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2):
      # Assumption : e1, e2, vp1, vp2, m1, m2, T1, T2 are np.arrays obtained from MCMC posteriors
      # The masses m1 and m2 are relative to the star. e1 and e2 are the eccentricities.
      # vp1 and vp2 are the longitudes of the periapsis and T1 and T2 are the periods
      # The function returns the triplet (X, Y, \delta). The resonance is p:p+1

      # Period of inner planet is normalized to 1
      T2    = T2/T1
      T1    = T1/T1

      # Getting semi-major axes and Lambda
      G     = 4.*m.pi**2
      beta1 = m1/(1. + m1)
      beta2 = m2/(1. + m2)
      mu1   = G*(1. + m1)
      mu2   = G*(1. + m2) 
      n1    = 2.*m.pi/T1
      n2    = 2.*m.pi/T2
      n10   = 2.*m.pi
      n20   = p*n10/(p + 1)
      a1    = (mu1/n1**2)**(1./3.)
      a2    = (mu2/n2**2)**(1./3.)
      Lbd1  = beta1*np.sqrt(mu1*a1)
      Lbd2  = beta2*np.sqrt(mu2*a2)
      D1    = Lbd1*(1. - np.sqrt(1. - e1**2))
      D2    = Lbd2*(1. - np.sqrt(1. - e2**2))

      #Defining the exact resonance
      a10   = (mu1/n10**2)**(1./3.)
      a20   = (mu2/n20**2)**(1./3.)
      Lbd10 = beta1*np.sqrt(mu1*a10)
      Lbd20 = beta2*np.sqrt(mu2*a20)

      # Getting G and Gamma and normalizing
      G     = Lbd1 + Lbd2 - D1 - D2
      Gamma = (p + 1)*Lbd1 + p*Lbd2
      g     = G/Gamma
      d1    = D1/Gamma
      d2    = D2/Gamma
      C1    = Gamma/Lbd10
      C2    = Gamma/Lbd20

      #Getting alpha, beta, gamma, delta, R and S
      f1    = f1s[p - 1]
      f2    = f2s[p - 1]
      R     = (f1**2*C1*d1 + f2**2*C2*d2 + 2.*f1*f2*np.sqrt(C1*d1*C2*d2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)
      S     = (f1**2*C1*d2 + f2**2*C2*d1 - 2.*f1*f2*np.sqrt(C1*d1*C2*d2)*np.cos(vp1 - vp2))/(f1**2*C1 + f2**2*C2)
      alpha = -3.*n10*p*((g + S)*(p*C1 + (p + 1)*C2) - C1 - C2)
      beta  = 1.5*n10*p*(p*C1 + (p + 1)*C2)
      gamma = m1*n20/C2*np.sqrt(f1**2*C1 + f2**2*C2)
      delta = alpha*(4./(27.*beta*gamma**2))**(1./3.)

      #Getting X and Y
      K     = (2.*beta/gamma)**(-2./3.)
      omega = beta*(2.*beta/gamma)**(-4./3.)
      Sigma = R/K
      xi    = -p*lbd1 + (p + 1)*lbd2
      sig1  = xi - vp1
      sig2  = xi - vp2
      u1    = np.sqrt(2.*d1)*np.cos(sig1)
      u2    = np.sqrt(2.*d2)*np.cos(sig2)
      v1    = np.sqrt(2.*d1)*np.sin(sig1)
      v2    = np.sqrt(2.*d2)*np.sin(sig2)
      z     = f2*np.sqrt(C2)/(f1*np.sqrt(C1))
      cophi = 1./np.sqrt(1. + z**2)
      siphi = z /np.sqrt(1. + z**2)
      x1    = cophi*u1 + siphi*u2
      y1    = cophi*v1 + siphi*v2
      cossig= x1/np.sqrt(2.*R)
      sinsig= y1/np.sqrt(2.*R)
      X     = np.sqrt(2.*Sigma)*cossig
      Y     = np.sqrt(2.*Sigma)*sinsig

      return [X, Y, delta]
      
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
      if (tD < 0. or tD > 1. or tX < 0. or tX > 1.):
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
      if (abs(nu1 - nu2) > 1.e-4):
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
            [X_target, delta_target] = dicho('external', 1., 49.999, nu)
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


delta_min = -8#-3.3 #To be chosen by trial and error
delta_max = 26#17.4
normX_min = 0#-3.2
normX_max = 2.5#7.4

f1s = [ 1.1904936978,  2.0252226899,  2.8404318567,  3.6496182441,  4.4561427851]
f2s = [-0.4283898341, -2.4840051833, -3.2832567218, -4.0837053718, -4.8847062975]

#Plotting
fig, ((ax1)) = py.subplots(1, 1, sharex=True, sharey=True, gridspec_kw={'width_ratios': [1]}, constrained_layout=False)
py.subplots_adjust(left=0.26, right=0.71, bottom=0.1, top=0.95)
ax1.axhline(y = 1, xmin = delta_min, xmax = delta_max, color = 'lightpink', linewidth = 6, linestyle = '-', alpha = 1, label = 'Continued separatrix')

for pair in range(N_pairs):

      print ('Pair n° ', pair)
      dummy = dummy_fst_col[pair]
      cols  = [5*pair_in_chain[pair][0] + 0 + dummy, 5*pair_in_chain[pair][0] + 1 + dummy, 5*pair_in_chain[pair][0] + 2 + dummy, 5*pair_in_chain[pair][0] + 3 + dummy, 5*pair_in_chain[pair][0] + 4 + dummy, 5*pair_in_chain[pair][1] + 0 + dummy, 5*pair_in_chain[pair][1] + 1 + dummy, 5*pair_in_chain[pair][1] + 2 + dummy, 5*pair_in_chain[pair][1] + 3 + dummy, 5*pair_in_chain[pair][1] + 4 + dummy]
      lbd1, T1, k1, h1, m1, lbd2, T2, k2, h2, m2 = np.loadtxt(rad[pair] + paths[pair], dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array(cols))
      lbd1 = lbd1*m.pi/180.
      lbd2 = lbd2*m.pi/180.
      e1   = np.sqrt(k1**2 + h1**2)
      vp1  = np.arctan2(h1, k1)
      e2   = np.sqrt(k2**2 + h2**2)
      vp2  = np.arctan2(h2, k2)
      p    = resonances[pair]
      n10  = 2.*m.pi
      n20  = p*n10/(p + 1)
      m1   = 10.**(m1)
      m2   = 10.**(m2)


      #Getting the posterior points
      x1s = []
      x2s = []
      N   = len(m1)
      [X, Y, Ds] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2)
      for i in range(N):
            [x1, x2] = levelLine2X(X[i], Y[i], Ds[i])
            x1s.append(x1)
            x2s.append(x2)
            if ((i + 1)%20 == 0):
                  print(i + 1, '/', N)
      x1s = np.array(x1s)
      x2s = np.array(x2s)
      Ds  = np.array(Ds)
      NormalizedX = []
      NormalizedD = []
      for i in range(len(Ds)):
            X = max(x1s[i], x2s[i])
            x = min(x1s[i], x2s[i])
            d = Ds[i]
            nu1 = freq(X, d)
            nu2 = freq(x, d)
            if (abs(nu1 - nu2) > 1.e-4):
                  print('Problem with frequencies: nu1 = ', nu1, 'nu2 = ', nu2)
            nu = 0.5*(nu1 + nu2)
            if (d <= 50. and d >= -20. and X <= 14. and x > -14.):
                  [xmin, xmax, xint, xext, xhyp] = topologie(d)
                  [X_target, delta_target, distance] = target(X, x, d)
                  NormalizedD.append(delta_target)
                  if (nu <= 3.**(0.5)*2.**(2./3.) and X >= xmin):
                        S = xmax
                        D = np.sqrt((S - 2.**(1./3.))**2 + d**2)
                        NormalizedX.append(distance/D)
                  elif (X >= xmin):
                        if ((X <= xmax and d >= 1.) or (X <= xmax and d < 1. and d > 0. and x >= 0.)):
                              S = xmax
                              D = np.sqrt((S - 2.**(1./3.))**2 + d**2)
                              NormalizedX.append(distance/D)
                        else:
                              [S, dlt] = dicho('max', -20., 1., nu)
                              D = np.sqrt((S - X_target)**2 + (dlt - delta_target)**2)
                              NormalizedX.append(distance/D)
                  else:
                        S = xmin
                        D = np.sqrt((S + 1.)**2 + (d - 1.)**2)
                        NormalizedX.append(-distance/D)
            else:
                  print("Warning : Point n°", i, " is out of bound\n")
      NormalizedD = np.array(NormalizedD)
      NormalizedX = np.array(NormalizedX)
      if (error_bar):
            averageD = np.average(NormalizedD)
            averageX = np.average(NormalizedX)
            stdD     = np.std    (NormalizedD)
            stdX     = np.std    (NormalizedX)
            ax1.errorbar(averageD, averageX, yerr=stdX, xerr=stdD, marker = marker[pair], markersize = 14., fmt='', color = plot_color[pair], ecolor=plot_color[pair], elinewidth=2., capsize=None, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None, data=None, label = lbl[pair])
      else:
            ax1.scatter(NormalizedD, NormalizedX, c = plot_color[pair], marker = marker[pair],  s = 80, alpha = 0.4, label = lbl[pair])

ax1.set_xlim(xmin = delta_min, xmax = delta_max)
ax1.set_ylim(ymin = normX_min,     ymax = normX_max)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.set_xlabel(xlabel="$\delta$", fontsize=25, labelpad=3)
ax1.set_ylabel(ylabel="Normalized distance to resonance", fontsize=25, labelpad=4, rotation=90)
ax1.grid(linewidth=0.3, alpha = 0.45)
py.legend(fontsize = 20)
py.show()


'''ax1.set_xlim(xmin = delta_min, xmax = delta_max)
ax1.set_ylim(ymin = X_min,     ymax = X_max)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.set_xlabel(xlabel="$\delta$", fontsize=25, labelpad=3)
ax1.set_ylabel(ylabel="$X$",      fontsize=25, labelpad=4, rotation=0)
#ax1.grid(linewidth=0.3, alpha = 0.45)
ax1.plot(delta_min - 1., 0., color = 'grey', linewidth = 2., linestyle = '-', alpha = 0.5, label = 'Frequency')
ldn = [2.2, 2.6, 3.**(0.5)*2.**(2./3.), 3., 3.5, 4., 5., 6., 8., 12., 18., 26., 36., 48., 62.]
cs = ax1.contour(X, Y, Z, ldn, linewidths = 2., colors = ['grey'], linestyles = 'solid', alpha = 0.5)
#cs = ax1.contour(X, Y, Z, linewidths = 2., colors = ['grey'], linestyles = 'solid', alpha = 0.5)
ax1.clabel(cs, inline=1, inline_spacing=10, fmt='%.1f', manual=True, fontsize=20)
ax1.plot(delt, Xint, color = 'black',     linewidth = 4, linestyle = '-',  alpha = 1, label = 'Elliptic')
ax1.plot(delt[delt >= 1.], Xext[delt >= 1.], color = 'black',     linewidth = 4, linestyle = '-',  alpha = 1)
ax1.plot(delt[delt >= 1.], Xhyp[delt >= 1.], color = 'lightpink', linewidth = 4, linestyle = '--', alpha = 1, label = 'Hyperbolic')
ax1.plot(delt, Xmin, color = 'lightpink', linewidth = 4, linestyle = '-',  alpha = 1, label = 'Separatrix')
ax1.plot(delt, Xmax, color = 'lightpink', linewidth = 4, linestyle = '-',  alpha = 1)
py.legend(fontsize = 25)
py.show()'''




