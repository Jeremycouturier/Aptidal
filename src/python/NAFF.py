import math as m
import matplotlib.pyplot as py
import matplotlib
import numpy as np

########                    Trying to recover the fundamental frequencies and amplitudes of the function                     ########
######## f(t) = a_0 + 2*a_1*cos(nu_1*t) - 2*a_2*sin(nu_2*t) + 2*a_3*cos(nu_3*t) with the NAFF algorithm (Laskar et al. 1992) ########

a_0  = 1.0
a_1  = 2.0
a_2  = 3.0
a_3  = 1.5
nu_1 = 0.47267
nu_2 = 1.10372
nu_3 = 1.42289
T1   = 2.0*m.pi/(nu_1)
T2   = 2.0*m.pi/(nu_2)
T3   = 2.0*m.pi/(nu_3)
T    = 10*max(T1, T2, T3)
nu_0 = m.pi/T


N = 50000 #Number of points
omega_min = -2.0
omega_max = 2.0


def sinx(x): #Defining the function x -> sin(x)/x
      return np.sinc(x / np.pi)
      
def f(t):    #Defining the function f(t)
      return a_0 + 2.0*a_1*np.cos(nu_1*t) - 2.0*a_2*np.sin(nu_2*t) + 2.0*a_3*np.cos(nu_3*t)
      

######## Frequencies omega ########
X = np.linspace(omega_min,omega_max,N)


######## Defining the Fourier-Hanning transform ########
def hat_f(nu, a, omega):
      return a/(1.0-(nu-omega)**2/nu_0**2)*sinx((nu-omega)*T)


######## <f(t),exp(i*omega*t)> ########
Y = hat_f(0.0, a_0, X) + hat_f(nu_1, a_1, X) + hat_f(-nu_1, a_1, X) + hat_f(nu_2, a_2, X) + hat_f(-nu_2, -a_2, X) + hat_f(nu_3, a_3, X) + hat_f(-nu_3, a_3, X)


######## Naive determination of the average ########
A_naive = 1.0+1.0/(2.0*T)*(4.0/nu_1*np.sin(2.0*nu_1*T)+6.0/nu_2*np.cos(2.0*nu_2*T)+3.0/nu_3*np.sin(2.0*nu_3*T)-6.0/nu_2)
print("A_naive = ", A_naive)


######## Determination of the average thanks to the NAFF method                                              ########
######## I look for the maximum of <f(t),exp(i*omega*t)> in the vicinity of omega = 0 with a gradient ascent ########
def F(omg):
      return hat_f(0.0, a_0, omg) + hat_f(nu_1, a_1, omg) + hat_f(-nu_1, a_1, omg) + hat_f(nu_2, a_2, omg) + hat_f(-nu_2, -a_2, omg) + hat_f(nu_3, a_3, omg) + hat_f(-nu_3, a_3, omg)
step = 0.0001
omega = 0.0
gradient = 1.0
while (gradient > 0.00000001):
      gradient = (F(omega + 0.5*step) - F(omega - 0.5*step))/step
      omega    = omega + step*gradient

print("omega = ", omega)
A_naff = F(omega)
print("A_naff  = ", A_naff)
A_naff2 = F(0.0)
print("A_naff2  = ", A_naff2)


######## Precision of the average as a function of T ########
TT = np.linspace(20.0,1000.0,15000)
error_naive = []
error_naff  = []
error_naff2 = []
i = 0

for T in TT:
      nu_0 = m.pi/T
      
      ######## Naive determination ########
      A_naive = a_0+1.0/(2.0*T)*(4.0/nu_1*np.sin(2.0*nu_1*T)+6.0/nu_2*np.cos(2.0*nu_2*T)+3.0/nu_3*np.sin(2.0*nu_3*T)-6.0/nu_2)
      
      ######## Subtle determination with Fourier-Hanning transform (NAFF) ########
      step = 0.0001
      omega = 0.0
      gradient = 1.0
      while (gradient > 0.00000001):
            gradient = (F(omega + 0.5*step) - F(omega - 0.5*step))/step
            omega    = omega + step*gradient
      A_naff = F(omega)
      
      ######## Simply evaluating <f(t),exp(i*omega*t)> at omega = 0 ########
      A_naff2 = F(0.0)
      
      error_naive.append(abs(A_naive-a_0))
      error_naff.append(abs(A_naff-a_0))
      error_naff2.append(abs(A_naff2-a_0))
      print("progress = ", i/15000*100, "%")
      i = i+1
      
error_naive = np.log10(np.array(error_naive))
error_naff  = np.log10(np.array(error_naff))
error_naff2 = np.log10(np.array(error_naff2))
Tm3 = np.log10(250.0*TT**(-3))
Tm1 = np.log10(5*TT**(-1))
TT          = np.log10(TT)


my_yticks=[r'$10^{-8}$',r'$10^{-7}$',r'$10^{-6}$',r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',r'$1$']
py.yticks([-8,-7,-6,-5,-4,-3,-2,-1,0],my_yticks)
my_xticks=[r'$10^{1}$',r'$10^{1.2}$',r'$10^{1.4}$',r'$10^{1.6}$',r'$10^{1.8}$',r'$10^{2}$',r'$10^{2.2}$',r'$10^{2.4}$',r'$10^{2.6}$',r'$10^{2.8}$',r'$10^{3}$']
py.xticks([1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3],my_xticks)


py.plot(TT, error_naive, "-", color = "red",   linewidth = 5, label="Naive average")
py.plot(TT, error_naff,  "-", color = "orange",linewidth = 5, label=r"Maximum of $\left<f(t),e^{i\omega t}\right>$ around $0$")
py.plot(TT, error_naff2, "-", color = "blue",  linewidth = 5, label=r"$\left<f(t),1\right>$")
py.plot(TT, Tm3,         "-", color = "gold",  linewidth = 3, label=r"$250/T^3$")
py.plot(TT, Tm1,         "-", color = "green", linewidth = 3, label=r"$5/T$")
py.xticks(fontsize=25)
py.yticks(fontsize=25)
py.xlabel(r"$T$", fontsize=25)
py.ylabel("Relative error in the determination of the average", fontsize=25)
py.grid(linewidth=0.4)
py.legend(fontsize = 25)
py.show()

      
'''
py.plot(X, Y, "-", color = "blue", linewidth = 5)
py.xticks(fontsize=25)
py.yticks(fontsize=25)
py.xlabel(r"$\omega$", fontsize=25)
py.ylabel(r"$\left<f(t),e^{i\omega t}\right>$", fontsize=25)
py.grid(linewidth=0.4)
py.legend(fontsize = 25)
py.show()
'''






            
            
            
