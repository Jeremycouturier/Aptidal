######## A script to convert elliptic elements to a level line of the Hamiltonian                       ########
######## 3*delta*Sigma - Sigma^2 + 2sqrt(2Sigma)cos(sigma) of the second fundamental model of resonance ########

import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np

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
      mu1   = G* (1. + m1)
      mu2   = G* (1. + m2) 
      n1    = 2.*m.pi/T1
      n2    = 2.*m.pi/T2
      a1    = (mu1/n1**2)**(1./3.)
      a2    = (mu2/n2**2)**(1./3.)
      Lbd1  = beta1*np.sqrt(mu1*a1)
      Lbd2  = beta2*np.sqrt(mu2*a2)
      D1    = Lbd1*(1. - np.sqrt(1. - e1**2))
      D2    = Lbd2*(1. - np.sqrt(1. - e2**2))

      #Defining the exact resonance
      n10   = 2.*m.pi
      n20   = p*n10/(p + 1)
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

e1  = 0.0382
e2  = 0.0603
vp1 = 2.719
vp2 = 0.370
lbd1= -2.055
lbd2= 0.262
m1  = 0.000271
m2  = 0.000403
p   = 2
T1  = 3.6718
T2  = 5.5087

[X, Y, delta] = ell2SFM(p, e1, e2, vp1, vp2, m1, m2, T1, T2, lbd1, lbd2)
print(X, Y)
print(delta)

'''
delta = 0.8  #Parameter of the model. Bifurcation is at delta = 1

Xmin = -3.26 #To be chosen by trial and error
Xmax = 3.9
Ymin = -3.64
Ymax = 3.64


if (delta < 1.):
      u      = (1. + m.sqrt(1. - delta**3))**(1./3.)
      v      =  1. - m.sqrt(1. - delta**3)
      if (v < 0.):
            v = -((-v)**(1./3.))
      else:
            v = v**(1./3.)
      X1     = u + v
      X2     = X1
      X3     = X1

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
      X2 = np.real(j*u +jb*v)
      X3 = np.real(jb*u +j*v)


######## Finding the hyperbolic fixed point ########
if   (delta >= 1 and X2**2 >= delta and X2**2 <= 3.*delta):
      X_ell = X3
      X_hyp = X2
elif (delta >= 1 and X3**2 >= delta and X3**2 <= 3.*delta):
      X_ell = X2
      X_hyp = X3
elif (delta >= 1):
      print("No hyperbolic equilibrium found")

######## Defining the Hamiltonian ########
def H(X,Y):
      return 2.*X - 0.25*(X**2 + Y**2 - 3.*delta)**2 + 9./4.*delta**2


######## Defining the level lines of interest of the Hamiltonian ########
if (delta < 1.):
      ldn = [H(X1, 0.), H(1.18*X1, 0.), H(1.36*X1, 0.), H(1.04*X1, 0.), H(1.5*X1, 0.)]
elif (delta >= 1. and delta <= 1.5):
      ldn = [H(X_hyp, 0.), H(1.18*X1, 0.), H(1.36*X1, 0.), H(1.45*X1, 0.), H(1.04*X1, 0.), H(1.56*X1, 0.)]
else:
      ldn = [H(X_hyp, 0.), H(1.18*X1, 0.), H(1.3*X1, 0.), H(1.04*X1, 0.), H(1.5*X1, 0.), H(0.5*(X_ell + X_hyp), 0.)]

ldn.sort()

def contour():

      x=np.linspace(Xmin,Xmax,2000)
      y=np.linspace(Ymin,Ymax,2000)
      Z=np.zeros((2000,2000))

      for i in range(2000):
            for j in range(2000):
                  #Z[i,j]=max(H(x[j],y[i]),H(0.99*Xmax,0.0))
                  Z[i,j]=max(H(x[j],y[i]),H(1.62*X1,0.0))
      
      py.xlim([Xmin,Xmax])
      py.ylim([Ymin,Ymax])
      ax = py.gca()
      ax.set_aspect('equal', adjustable='box')
      X,Y = np.meshgrid(x,y)
      ax.plot(X1, 0., 'ko')
      if (delta >= 1):
            ax.plot(X_ell, 0., 'ko')
      ax.grid(linewidth=0.3)
      py.pcolormesh(X,Y,Z,cmap='ocean_r')
      cbar=py.colorbar(orientation="vertical", pad=0.04)
      cbar.ax.tick_params(labelsize=25)
      py.xticks(fontsize=25)
      py.yticks(fontsize=25)
      cs=py.contour(X,Y,Z,ldn,linewidths=3,colors=['grey'],linestyles='solid')
      py.clabel(cs, inline=1, inline_spacing=10, fmt='%.3f', manual=True, fontsize=25)
      py.xlabel(r'$X$', fontsize=25)
      py.ylabel(r'$Y$', fontsize=25)
      py.show()
            
contour()'''  
            
            
            
