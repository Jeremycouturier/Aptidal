import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np
import random

#Defining the system
path = "/home/jeremy/Documents/PosteriorLeleu/Kepler1705bc"
lbl  = 'Kepler-1705'
p    = 4 #Resonance p : p+1
delta_min = -3.3 #To be chosen by trial and error
delta_max = 17.4
X_min     = -3.2
X_max     = 7.4

eps_mass_ratio_bool = 1 #Set to 1 if masses are given as log10(m_j) and m1/(m1 + m2) and to 0 if masses are given as m1 and m2 directly

lbd1, T1, k1, h1, log10_eps, lbd2, T2, k2, h2, m12 = np.loadtxt(path, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
lbd1 = lbd1*m.pi/180.
lbd2 = lbd2*m.pi/180.
e1   = np.sqrt(k1**2 + h1**2)
vp1  = np.arctan2(h1, k1)
e2   = np.sqrt(k2**2 + h2**2)
vp2  = np.arctan2(h2, k2)
n10  = 2.*m.pi
n20  = p*n10/(p + 1)
if (eps_mass_ratio_bool):
      eps = 10.**(log10_eps)
      m1  = eps*m12
      m2  = eps - m1
else:
      m1  = log10_eps
      m2  = m12
      

f1s = [ 1.1904936978,  2.0252226899,  2.8404318567,  3.6496182441,  4.4561427851]
f2s = [-0.4283898341, -2.4840051833, -3.2832567218, -4.0837053718, -4.8847062975]

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
      while (sol < 2):
            X0 = random.random()
            X0 = 2.*D*X0 - D
            h  = 1.
            it = 0
            cv = 1
            while(abs(h) > 1.e-11):
                  num   = X0**4 - 6.*delta*X0**2 - 8.*X0 + 4.*H
                  denom = 4.*X0**3 - 12.*delta*X0 - 8.
                  h     = -num/denom
                  X0    = X0 + h
                  it    = it + 1
                  if (it > 10):
                        h  = 1.e-200
                        cv = 0
            new_sol = cv
            for i in range(sol):
                  if (abs(X0 - Sol[i]) < 1.e-9):
                        new_sol = 0
            sol = sol + new_sol
            if (new_sol):
                  Sol.append(X0)
      return [Sol[0], Sol[1]]


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
                  while(abs(h) > 1.e-11):
                        num   = X0**4 - 6.*delta*X0**2 - 8.*X0 + 4.*H
                        denom = 4.*X0**3 - 12.*delta*X0 - 8.
                        h     = -num/denom
                        X0    = X0 + h
                  new_sol = 1
                  for i in range(sol):
                        if (abs(X0 - Sol[i]) < 1.e-6):
                              new_sol = 0
                  sol = sol + new_sol
                  if (new_sol):
                        Sol.append(X0)
            Sol.sort()
            return [X1, X_ell, Sol[0], Sol[1], Sol[2]]



#Defining the topology of the phase space
Delta = np.linspace(delta_min, delta_max, 3000)
X1s   = []
X2s   = []
X3s   = []
X4s   = []
X5s   = []
for delta in Delta:
      if (delta <= 1):
            [X1]                 = topologie(delta)
      else:
            [X1, X2, X3, X4, X5] = topologie(delta)
      if (delta <= 1):
            X1s.append(X1)
            X2s.append(0.)
            X3s.append(0.)
            X4s.append(0.)
            X5s.append(0.)
      else:
            X1s.append(X1)
            X2s.append(X2)
            X3s.append(X3)
            X4s.append(X4)
            X5s.append(X5)
X1s = np.array(X1s)
X2s = np.array(X2s)
X3s = np.array(X3s)
X4s = np.array(X4s)
X5s = np.array(X5s)
Del = Delta[X2s != 0.]
X2s = X2s  [X2s != 0.]
X3s = X3s  [X3s != 0.]
X4s = X4s  [X4s != 0.]
X5s = X5s  [X5s != 0.]


#Getting the posterior points
x1s = []
x2s = []
N   = len(m1)
[X, Y, Ds] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2)
for i in range(N):
      [x1, x2] = levelLine2X(X[i], Y[i], Ds[i])
      x1s.append(x1)
      x2s.append(x2)
      print(i + 1, '/', N)
x1s = np.array(x1s)
x2s = np.array(x2s)
Ds  = np.array(Ds)


#Getting the level lines of (p+1)*n2/(p*n1)
m1_bar    = 10.**(np.average(np.log10(m1))) #The average of m1 and m2 over the posteriors are consider
m2_bar    = 10.**(np.average(np.log10(m2)))
C1_bar    = (p + 1)/p + m2_bar/m1_bar*((p + 1)/p)**(1./3.)
C2_bar    = 1. + m1_bar/m2_bar*((p + 1)/p)**(2./3.)
beta_bar  = 1.5*n10*p*(p*C1_bar + (p + 1)*C2_bar)
gamma_bar = m1_bar*n20/C2_bar*m.sqrt(f1s[p - 1]**2*C1_bar + f2s[p - 1]**2*C2_bar)
dlt = np.linspace(delta_min, delta_max, 2000)
XXX = np.linspace(X_min,     X_max,     2000)
Z   = np.zeros((2000,2000))
for i in range(2000):
      for j in range(2000):
            R   = 0.5*(2.*beta_bar/gamma_bar)**(-2./3.)*XXX[i]**2
            alp = (4./(27.*beta_bar*gamma_bar**2))**(-1./3.)*dlt[j]
            gpS = (C1_bar + C2_bar - alp/(3.*n10*p))/(p*C1_bar + (p + 1)*C2_bar)
            eps = gpS + R
            Z[i, j] = (p + 1)/p*(1. - 3.*(C1_bar + C2_bar)*(p*eps - 1.) - 3.*C2_bar*eps)**(-1.)
X,Y = np.meshgrid(dlt, XXX)

#Plotting
fig, ((ax1)) = py.subplots(1, 1, sharex=True, sharey=True, gridspec_kw={'width_ratios': [1]}, constrained_layout=False)
py.subplots_adjust(left=0.245, right=0.81, bottom=0.12, top=0.95)
ax1.plot(Delta, X1s, color = 'black', linewidth = 3, linestyle = '-',  alpha = 1,   label = 'Resonance center')
ax1.plot(Del,   X2s, color = 'black', linewidth = 3, linestyle = '--', alpha = 1,   label = 'External circulation')
ax1.plot(Del,   X3s, color = 'red',   linewidth = 3, linestyle = '--', alpha = 1,   label = 'Separatrix (fixed point)')
ax1.plot(Del,   X4s, color = 'red',   linewidth = 3, linestyle = '-',  alpha = 1,   label = 'Separatrix')
ax1.plot(Del,   X5s, color = 'red',   linewidth = 3, linestyle = '-',  alpha = 1)
ax1.scatter(Ds, x1s, c     = 'green', marker = "o",  s = 80,           alpha = 0.6, label = lbl)
ax1.scatter(Ds, x2s, c     = 'green', marker = "o",  s = 80,           alpha = 0.6)
ax1.set_xlim(xmin = delta_min, xmax = delta_max)
ax1.set_ylim(ymin = X_min,     ymax = X_max)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.set_xlabel(xlabel="$\delta$", fontsize=25, labelpad=3)
ax1.set_ylabel(ylabel="$X$",      fontsize=25, labelpad=4, rotation=0)
ax1.grid(linewidth=0.3, alpha = 0.45)
ax1.plot(delta_min - 1., 0., color = 'grey', linewidth = 2., linestyle = '-', alpha = 0.5, label = r'$T_2/T_1$')
ldn = [1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28]
cs = ax1.contour(X, Y, Z, ldn, linewidths = 2., colors = ['grey'], linestyles = 'solid', alpha = 0.5)
ax1.clabel(cs, inline=1, inline_spacing=10, fmt='%.2f', manual=True, fontsize=25)
py.legend(fontsize = 25)
py.show()
















