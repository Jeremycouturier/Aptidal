import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np
import random

straight    = 1
error_bar   = 1
names_print = 1

radical     = '/home/jeremy/Documents/CleanPosteriors/ForClassification/'
sample_path = ['Agol_TRAPPIST1','Hadden_Kepler-176','Hadden_Kepler-307','Hadden_Kepler-53','Almenara_Kepler-117','Hadden_Kepler-177','Hadden_Kepler-310','Hadden_Kepler-54','Almenara_Kepler-419','Hadden_Kepler-18','Hadden_Kepler-31','Hadden_Kepler-549','Almenara_WASP-148','Hadden_Kepler-223','Hadden_Kepler-32','Hadden_Kepler-55','Almenara_WASP-47','Hadden_Kepler-23','Hadden_Kepler-324','Hadden_Kepler-56','Borsato_TOI-1130','Hadden_Kepler-238','Hadden_Kepler-33','Hadden_Kepler-57','Hadden_Kepler-24','Hadden_Kepler-345','Hadden_Kepler-58','Dai_TOI1136','Hadden_Kepler-25','Hadden_Kepler-359','Hadden_Kepler-60','Hadden_Kepler-105','Hadden_Kepler-26','Hadden_Kepler-36','Hadden_Kepler-79','Hadden_Kepler-11','Hadden_Kepler-27','Hadden_Kepler-396','Hadden_Kepler-80','Hadden_Kepler-1126','Hadden_Kepler-277','Hadden_Kepler-444','Hadden_Kepler-81','Hadden_Kepler-114','Hadden_Kepler-279','Hadden_Kepler-48','Hadden_Kepler-84','Hadden_Kepler-122','Hadden_Kepler-28','Hadden_Kepler-49','Hadden_Kepler-85','Hadden_Kepler-127','Hadden_Kepler-29','Hadden_Kepler-51','Hadden_Kepler-9','Hadden_Kepler-128','Hadden_Kepler-30','Hadden_Kepler-52','Lindor_Sol','Hadden_Kepler-138','Hadden_Kepler-305','Hadden_Kepler-526','Wang_TOI1338']
names       = ['TRAPPIST1','K176','K307','K53','K117','K177','K310','K54','K419','K18','K31','K549','WASP148','K223','K32','K55','WASP47','K23','K324','K56','TOI1130','K238','K33','K57','K24','K345','K58','TOI1136','K25','K359','K60','K105','K26','K36','K79','K11','K27','K396','K80','K1126','K277','K444','K81','K114','K279','K48','K84','K122','K28','K49','K85','K127','K29','K51','K9','K128','K30','K52','VenusEarth','K138','K305','K526','TOI1338']

#sample_path = ['Hadden_Kepler-176','Hadden_Kepler-307','Hadden_Kepler-53','Hadden_Kepler-177','Hadden_Kepler-54', 'Hadden_Kepler-18']
#names       = ['K176','K307','K53','K177','K54', 'K18']

delta_min   = -30. #To be chosen by trial and error
delta_max   = 25.
normX_min   = 0.
normX_max   = 4.5
color_min   = -6.
color_max   = -2.

'''delta_min   = -25 #To be chosen by trial and error
delta_max   = 40
normX_min   = -12
normX_max   = 12'''

#Hamiltonian coefficients
f1s = [ 1.1904936978,  2.0252226899,  2.8404318567,  3.6496182441,  4.4561427851]
f2s = [-0.4283898341, -2.4840051833, -3.2832567218, -4.0837053718, -4.8847062975]


delt, Xmin, Xmax, Xint, Xext, Xhyp = np.loadtxt('../continuedSeparatrix.txt', dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5]))
delt = np.flip(delt)
Xmin = np.flip(Xmin)
Xmax = np.flip(Xmax)
Xint = np.flip(Xint)
Xext = np.flip(Xext)
Xhyp = np.flip(Xhyp)

freqs = np.loadtxt('/home/jeremy/Documents/GRSW/frequency_fixed.csv', dtype = np.float64, delimiter=',', unpack=True)

def f(e):
      return -1.0/9.0*(16.0*e**2-32.0*e+7.0)

def g(e):
      return 4.0*e*(1.0-e)
      
def h(e):
      return 64.0*e*(0.25-e)

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
            
def rk4(X, Y, delta, dt):
      #Starting from (X,Y), integrates with a Runge-Kutta 4 until Y has changed sign twice
      #Returns the two values of X corresponding to Y=0 and the frequency of the orbit

      X0    = X
      Y0    = 0.
      count = 0
      go_on = 1
      first_time = 1
      
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
            oldX  = X0
            X0    = X0 + dt/6.*(k1_x + 2.*k2_x + 2.*k3_x + k4_x)
            Y0    = Y0 + dt/6.*(k1_y + 2.*k2_y + 2.*k3_y + k4_y)
            count = count + 1
            if (count >= 10000):
                  return [0., 0., 1000.]
            if (Y0 == 0. or oldY*Y0 < 0.):
                  t    = abs(Y0/(Y0 - oldY))
                  tOld = (count - 1)*dt
                  tNew = count*dt
                  if (first_time):
                        first_time = 0
                        X1 = t*oldX + (1. - t)*X0
                        t1 = t*tOld + (1. - t)*tNew
                  else:
                        go_on = 0
                        X2 = t*oldX + (1. - t)*X0
                        t2 = t*tOld + (1. - t)*tNew
      Period    = 2.*(t2 - t1)
      frequency = 2.*np.pi/Period
      return [X1, X2, frequency]

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
      
def target(X, x, delta, nu):

      #Computes an returns [X_ell, X_sep, X_ref, delta_ell, delta_sep, delta_ref] where (X_ell, delta_ell) is the point on the branch of elliptic fixed points
      #with the same frequency (or with frequency sqrt(3)*2**(2/3)), (X_sep, delta_sep) is the point on the (continued) separatrix with either the same frequency or with the same delta
      #and (X_ref, delta_ref) is the point on the branch of elliptic fixed points used to define the distance to (X_sep, delta_sep)
      #(X, 0) and (x, 0) are on the same level line with X > x and have frequency nu
      
      [xmin, xmax, xint, xext, xhyp] = topologie(delta)
      if (nu <= 3.**(0.5)*2.**(2./3.) and X >= xmin):
            return [2.**(1./3.), xmax, 2.**(1./3.), 0., delta, 0.]
      if (X >= xmin):
            if ((X <= xmax and delta >= 1.) or (X <= xmax and delta < 1. and delta > 0. and x >= 0.)): #Corresponding branch is the right part of the internal branch
                  [X_ell, delta_ell] = dicho('internal', 0., 49.999, nu)
                  return [X_ell, xmax, 2.**(1./3.), delta_ell, delta, 0.]
            else: #Corresponding branch is the left part of the internal branch
                  [X_ell, delta_ell] = dicho('internal', -100., 0., nu)
                  [X_sep, delta_sep] = dicho('max',      -100., 1., nu)
                  return [X_ell, X_sep, X_ell, delta_ell, delta_sep, delta_ell]
      else: #Corresponding branch is the external circulation
            if (nu > 3.**(0.5)*2.**(2./3.)):
                  [X_ell, delta_ell] = dicho('internal', -100., 0., nu)
                  [X_sep, delta_sep] = dicho('max',      -100., 1., nu)
                  return [X_ell, xmin, -1., delta_ell, delta, 1.]
            else:
                  [X_sep, delta_sep] = dicho('max',      -100., 1., nu)
                  return [2.**(1./3.), X_sep, 2.**(1./3.), 0., delta_sep, 0.]

def freqOnBranch(branch, delta):
      #Computes the frequency at delta on the branch 'internal', 'external', or 'max'
      #Returns the frequency and the corresponding value of X
      if (delta >= 50.):
            raise Exception("delta >= 50 in function freqOnBranch")
      if (branch == 'internal'):
            if (delta >= -20.):
                  X  = topologie(delta)[2]
                  nu = m.sqrt((3.*delta - X**2)*(3.*delta - 3.*X**2))
            else:
                  X  = -2./(3.*delta)
                  nu = -3.*delta + 24./(27.*delta**2)
      elif (branch == 'external'):
                  X  = topologie(delta)[3]
                  nu = m.sqrt((3.*delta - X**2)*(3.*delta - 3.*X**2))
      else:
            if (delta >= -20.):
                  X  = topologie(delta)[1]
                  nu = freq(X, delta)
            else:
                  X  = m.sqrt(6.) - 2./(3.*delta)
                  nu = 6. - 3.*delta + 24./(27.*delta**2)
      return [nu, X]
      
def dicho(branch, min_delta, max_delta, frequency):
      #Returns the [X, delta] value on the corresponding branch having frequency equal to frequency
      #branch is either 'external', 'internal' or 'max'
      
      eps      = 2.5e-5
      
      [nu_left,  X_left]  = freqOnBranch(branch, min_delta)
      [nu_right, X_right] = freqOnBranch(branch, max_delta)
      if (frequency > max(nu_left, nu_right) or frequency < min(nu_left, nu_right)):
            print('frequency = ', frequency)
            print('X_left    = ', X_left)
            print('X_right   = ', X_right)
            print('min_delta = ', min_delta)
            print('max_delta = ', max_delta)
            print('nu_left   = ', nu_left)
            print('nu_right  = ', nu_right)
            raise Exception("Problem with dichotomy")
      
      while(abs(max_delta - min_delta) > eps):
            middle_delta = 0.5*(min_delta + max_delta)
            [middle_freq, middle_X] = freqOnBranch(branch, middle_delta)
            if ((frequency < nu_right and frequency > middle_freq) or (frequency > nu_right and frequency < middle_freq)):
                  X_left    = middle_X
                  nu_left   = middle_freq
                  min_delta = middle_delta
            else:
                  X_right   = middle_X
                  nu_right  = middle_freq
                  max_delta = middle_delta
      return [0.5*(X_left + X_right), 0.5*(min_delta + max_delta)]

def couleur(Min, Max, value):
      e = (value - Min)/(Max - Min)
      if (e < 0. or e > 1.):
            print("value = ", value)
            raise Exception("value is not between Min and Max in function couleur")
      if (e < 0.25):
            color = (0.0, g(e), h(e))
      elif (e <= 1.0):
            color = (f(e), g(e), 0.0)
      return color

def artisanal_colorbar(Min, Max, lbl, ticks, tick_lbls):

      bottom = normX_min + 0.25*(normX_max - normX_min)
      top    = normX_max - 0.05*(normX_max - normX_min)
      left   = delta_max - 0.09*(delta_max - delta_min)
      right  = delta_max - 0.05*(delta_max - delta_min)
      ecc = np.linspace(0.0,1.0,1000)
      where_to_plot = np.linspace(bottom,top,1000)
      for p in range(1000):
            if (p < 250):
                  color = (0.0, g(ecc[p]), h(ecc[p]))
            else:
                  color = (f(ecc[p]), g(ecc[p]), 0.0)
            py.plot([left,right], [where_to_plot[p],where_to_plot[p]], '-', color=color, linewidth=2)
      ax1.text(right + 0.2*(right - left), bottom, lbl.center(75), color="black", fontsize=20, alpha = 1., rotation=270)
      ax1.plot([left,left],[bottom,top],'-',color='black',linewidth=1)
      ax1.plot([right,right],[bottom,top],'-',color='black',linewidth=1)
      ax1.plot([left,right],[bottom,bottom],'-',color='black',linewidth=1)
      ax1.plot([left,right],[top,top],'-',color='black',linewidth=1)
      for i in range(len(ticks)):
            tick = ticks[i]
            tick_lbl = tick_lbls[i]
            if (tick >= Min and tick < Max):
                  p = 1000*(tick - Min)/(Max - Min)
                  p = round(p, 9)
                  p = int(p)
                  py.plot([left,right], [where_to_plot[p],where_to_plot[p]], '-', color="black", linewidth=0.5, alpha = 0.4)
                  ax1.text(left - 1.7*(right - left), where_to_plot[p] - 0.05, tick_lbl, fontsize=18, alpha = 1.)
            else:
                  print("Warning: Tick out of range")

#Plotting
fig, ((ax1)) = py.subplots(1, 1, sharex=True, sharey=True, gridspec_kw={'width_ratios': [1]}, constrained_layout=False)
py.subplots_adjust(left=0.26, right=0.71, bottom=0.1, top=0.95)

if (straight and error_bar):
      avD   = []
      avX   = []
      stddD = []
      stddX = []
      colr  = []
      lbl   = []

count = 0
for pth in sample_path:
      path   = radical + pth + '_0_samples.csv'
      print("Doing sample " + names[count])
      sample = np.loadtxt(path, dtype = np.float64, delimiter=',', unpack=True)
      row    = sample[:,0]
      n      = len(row)
      k      = n%8
      N      = (n - k)//8 #Number of planets
      ps     = [] # p such that the resonance is p : p+1 for the pair in question
      pairs  = [] # The corresponding pair of planets
      for j in range(N):
            for i in range(j): # Looking at pair [i,j]
                  Pi = sample[2 + 8*i,:]
                  Pj = sample[2 + 8*j,:]
                  PR = np.average(Pj/Pi)
                  if (abs((PR - 2.)/2.) < 0.03):
                        ps.append(1)
                        pairs.append([i,j])
                  elif (abs((PR - 1.5)/1.5) < 0.02):
                        ps.append(2)
                        pairs.append([i,j])
                  elif (abs((PR - 4./3.)/(4./3.)) < 0.015):
                        ps.append(3)
                        pairs.append([i,j])
                  elif (abs((PR - 1.25)/1.25) < 0.01):
                        ps.append(4)
                        pairs.append([i,j])
                  elif (abs((PR - 1.2)/1.2) < 0.01):
                        ps.append(5)
                        pairs.append([i,j])
      N_pairs = len(ps)
      for pair in range(N_pairs):
            I    = pairs[pair][0]
            J    = pairs[pair][1]
            print('pair '+str(I)+str(J))
            lbd1 = sample[1 + 8*I,:]
            lbd2 = sample[1 + 8*J,:]
            P1   = sample[2 + 8*I,:]
            P2   = sample[2 + 8*J,:]
            k1   = sample[3 + 8*I,:]
            k2   = sample[3 + 8*J,:]
            h1   = sample[4 + 8*I,:]
            h2   = sample[4 + 8*J,:]
            I1   = sample[5 + 8*I,:]
            I2   = sample[5 + 8*J,:]
            O1   = sample[6 + 8*I,:]
            O2   = sample[6 + 8*J,:]
            m1   = sample[7 + 8*I,:]
            m2   = sample[7 + 8*J,:]
            R1   = sample[8 + 8*I,:]
            R2   = sample[8 + 8*J,:]
            Mstar= sample[8*N + 1,:]
            Rstar= sample[8*N + 2,:]
            lbd1 = lbd1*np.pi/180.
            lbd2 = lbd2*np.pi/180.
            e1   = np.sqrt(k1**2 + h1**2)
            e2   = np.sqrt(k2**2 + h2**2)
            vp1  = np.arctan2(h1, k1)
            vp2  = np.arctan2(h2, k2)
            p    = ps[pair]
            n10  = 2.*m.pi
            n20  = p*n10/(p + 1)
            clr  = np.log10(np.average(m1) + np.average(m2))
            x1s  = []
            x2s  = []
            nus  = []
            n    = len(m1)
            [X, Y, Ds] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, P1, P2, lbd1, lbd2)
            for i in range(n):
                  XXX    = np.sqrt(X[i]**2 + Y[i]**2)
                  Dltt   = Ds[i]
                  if (Dltt >= 50. or Dltt <= -20. or XXX >= 14. or XXX <= -14.):
                        dt = 0.00001
                  else:
                        est_nu = freq(XXX, Dltt)
                        dt     = min(0.04, 2.*np.pi/est_nu/50.)
                  #[x1, x2] = levelLine2X(X[i], Y[i], Ds[i])
                  [x1, x2, nu] = rk4(X[i], Y[i], Ds[i], dt)
                  x1s.append(x1)
                  x2s.append(x2)
                  nus.append(nu)
            x1s = np.array(x1s)
            x2s = np.array(x2s)
            nus = np.array(nus)
            Ds  = np.array(Ds)
            if (straight):
                  NormalizedX = []
                  NormalizedD = []
                  for i in range(len(Ds)):
                        X  = max(x1s[i], x2s[i])
                        x  = min(x1s[i], x2s[i])
                        nu = nus[i]
                        d  = Ds[i]
                        if (d >= -20. and d < 50. and nu < 300.00008):
                              [X_ell, X_sep, X_ref, delta_ell, delta_sep, delta_ref] = target(X, x, d, nu)
                              distance = np.sqrt((X     - X_ell)**2 + (d         - delta_ell)**2)
                              normaD   = np.sqrt((X_sep - X_ref)**2 + (delta_sep - delta_ref)**2)
                              NormalizedD.append(delta_ell)
                              NormalizedX.append(distance/normaD)
                  NormalizedD = np.array(NormalizedD)
                  NormalizedX = np.array(NormalizedX)
                  if (len(NormalizedX)/len(m1) > 0.5): #Plot happens only if at least half of the posterior points were in bounds
                        if (error_bar):
                              averageD = np.average(NormalizedD)
                              averageX = np.average(NormalizedX)
                              stdD     = np.std    (NormalizedD)
                              stdX     = np.std    (NormalizedX)
                              colr.append(clr)
                              avD.append(averageD)
                              avX.append(averageX)
                              stddD.append(stdD)
                              stddX.append(stdX)
                              lbl.append(names[count]+'-'+str(I)+str(J))
                              #ax1.errorbar(averageD, averageX, yerr=stdX, xerr=stdD, marker = 'o', markersize = 14., fmt='', color = clr, ecolor = clr, elinewidth=2., capsize=None, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None, data=None)
                              #ax1.text(averageD, averageX, names[count]+'_'+str(I)+str(J), color='grey', size=1, alpha=0.4)
                        else:
                              clr = couleur(color_min, color_max, clr)
                              ax1.scatter(NormalizedD, NormalizedX, c = clr, marker = 'o',  s = 80, alpha = 0.4)
            else:
                  clr = couleur(color_min, color_max, clr)
                  ax1.scatter(Ds, x1s, c = clr, marker = 'o',  s = 80, alpha = 0.4, label = names[count] + '_' + str(I) + str(J))
                  ax1.scatter(Ds, x2s, c = clr, marker = 'o',  s = 80, alpha = 0.4)
      count = count + 1

if (straight and error_bar):
      colr  = np.array(colr)
      avD   = np.array(avD)
      avX   = np.array(avX)
      stddD = np.array(stddD)
      stddX = np.array(stddX)
      lbl   = np.array(lbl)
      '''delta_min = min(avD) - 0.2
      delta_max = max(max(avD) + 0.2, 4.)
      normX_min = 0.
      normX_max = max(avX) + 0.2'''
      #color_min = min(colr)
      #color_max = max(colr)
      color_min = -4.81
      color_max = -3.39
      for i in range(len(avD)):
            if (colr[i] < color_max and colr[i] > color_min):
                  clr = couleur(color_min, color_max, colr[i])
                  ax1.errorbar(avD[i], avX[i], yerr=stddX[i], xerr=stddD[i], marker = 'o', markersize = 14., fmt='', color = clr, ecolor = clr, elinewidth=2., capsize=None, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None, data=None)
                  if (names_print):
                        ax1.text(avD[i], avX[i], lbl[i], color='black', fontsize=10, alpha=1.)

if (straight == 0):
      dlt = np.arange(-20., 50., 0.002)
      Xs  = np.arange(-14., 14., 0.01)
      Z   = np.zeros((2800,35000))
      for i in range(2800):
            for j in range(35000):
                  Z[i, j] = freqs[i, j]
      X,Y = np.meshgrid(dlt, Xs)

ax1.set_xlim(xmin = delta_min, xmax = delta_max)
ax1.set_ylim(ymin = normX_min, ymax = normX_max)
ax1.tick_params(axis='both', which='major', labelsize=25)
#artisanal_colorbar(color_min, color_max, 'Mass of the pair (stellar masses)', [-4.4, -4.6, -4.8], [r'$10^{-4.4}$', r'$10^{-4.6}$', r'$10^{-4.8}$'])
#artisanal_colorbar(color_min, color_max, 'Mass of the pair (stellar masses)', [-9, -8, -7, -6, -5, -4, -3, -2, -1], [r'$10^{-9}$', r'$10^{-8}$', r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$'])
artisanal_colorbar(color_min, color_max, 'Mass of the pair (stellar masses)', [-4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4], [r'$10^{-4.8}$', r'$10^{-4.6}$', r'$10^{-4.4}$', r'$10^{-4.2}$', r'$10^{-4}$', r'$10^{-3.8}$', r'$10^{-3.6}$', r'$10^{-3.4}$'])
if (straight):
      ax1.set_xlabel(xlabel="Normalized $\delta$", fontsize=25, labelpad=3)
      ax1.set_ylabel(ylabel="Normalized distance to main elliptic branch", fontsize=25, labelpad=4, rotation=90)
else:
      ax1.set_xlabel(xlabel="$\delta$", fontsize=25, labelpad=3)
      ax1.set_ylabel(ylabel="$X$",      fontsize=25, labelpad=4, rotation=0)
if (straight):
      ax1.axhline(y =  1., xmin = delta_min, xmax = delta_max, color = 'lightpink', linewidth = 6, linestyle = '-', alpha = 1, label = 'Continued separatrix')
      ax1.fill_between([1., delta_max], 0., 1., where=None, interpolate=False, alpha = 0.15, color = 'grey')
else:
      ax1.plot(delta_min - 1., 0., color = 'grey', linewidth = 2., linestyle = '-', alpha = 0.5, label = 'Frequency')
      ldn = [2.2, 2.6, 3.**(0.5)*2.**(2./3.), 3., 3.5, 4., 5., 6., 8., 12., 18., 26., 36., 48., 62., 78., 96.]
      cs = ax1.contour(X, Y, Z, ldn, linewidths = 2., colors = ['grey'], linestyles = 'solid', alpha = 0.5)
      #cs = ax1.contour(X, Y, Z, linewidths = 2., colors = ['grey'], linestyles = 'solid', alpha = 0.5)
      ax1.clabel(cs, inline=1, inline_spacing=10, fmt='%.1f', manual=True, fontsize=20)
      ax1.plot(delt[delt >= 1.], Xext[delt >= 1.], color = 'black',     linewidth = 4, linestyle = '-',  alpha = 1)
      ax1.plot(delt[delt >= 1.], Xhyp[delt >= 1.], color = 'lightpink', linewidth = 4, linestyle = '--', alpha = 1, label = 'Hyperbolic')
      ax1.plot(delt, Xint, color = 'black',     linewidth = 4, linestyle = '-',  alpha = 1, label = 'Elliptic')
      ax1.plot(delt, Xmin, color = 'lightpink', linewidth = 4, linestyle = '-',  alpha = 1, label = 'Separatrix')
      ax1.plot(delt, Xmax, color = 'lightpink', linewidth = 4, linestyle = '-',  alpha = 1)
#ax1.grid(linewidth=0.3, alpha = 0.45)
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




