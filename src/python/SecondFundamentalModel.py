######## A script to plot the phase space of the one-parameter second fundamental model of resonance ########

import math as m
import cmath as cm
import matplotlib.pyplot as py
import matplotlib
import numpy as np

delta = 7  #Parameter of the model. Bifurcation is at delta = 1


Xmin = -7.1 #To be chosen by trial and error
Xmax = 7.1
Ymin = -7.5
Ymax = 7.5


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
      X2 = np.real(j*u + jb*v)
      X3 = np.real(jb*u + j*v)


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
      #ldn = [H(X_hyp, 0.), H(1.18*X1, 0.), H(1.3*X1, 0.), H(1.04*X1, 0.), H(1.5*X1, 0.), H(0.5*(X_ell + X_hyp), 0.)]
      ldn = [H(X_hyp, 0.), H(1.3*X1, 0.), H(1.04*X1, 0.), H(1.5*X1, 0.), H(0.5*(X_ell + X_hyp), 0.)]

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
            
contour()            
            
            
            
