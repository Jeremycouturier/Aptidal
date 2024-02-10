import math as m
import matplotlib.pyplot as py
import matplotlib
import numpy as np

######## Trying to compute the average of the function f(t) = a_0 + 2*a_1*cos(nu_1*t) - 2*a_2*sin(nu_2*t) + 2*a_3*cos(nu_3*t) ########
######## with the NAFF algorithm (Laskar et al. 1992) and different value for the power p of the weight function              ########

p    = 1

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
T    = 6   *max(T1, T2, T3)
dT   = 0.01*min(T1, T2, T3)
nu_0 = m.pi/T


N = 50000 #Number of points for plotting
omega_min = -2.0
omega_max = 2.0


def sinx(x): #Defining the function x -> sin(x)/x
      return np.sinc(x / np.pi)
      
def f(t):    #Defining the function f(t)
      return a_0 + 2.0*a_1*np.cos(nu_1*t) - 2.0*a_2*np.sin(nu_2*t) + 2.0*a_3*np.cos(nu_3*t)
      
def weight(t,p): #Defining the weight function
      Cp = 2**p*m.factorial(p)**2/m.factorial(2*p)
      return Cp*(1-np.cos(m.pi*t/T))**p
      
      
######## Numerical check of the formulae for < cos(nu t),1 > and < sin(nu t),1 > ########
def hp(x,p):
      prod = 1
      for i in range(1, p+1):
            prod = prod * (x**2 - i**2*m.pi**2)
      return (-1)**p*np.pi**(2*p)*m.factorial(p)**2*sinx(x)/prod
      
nu = 1.37162027
cos_nut_scalar_1 = m.cos(nu*T)*hp(nu*T,p) # < cos(nu*t),1 >
sin_nut_scalar_1 = m.sin(nu*T)*hp(nu*T,p) # < sin(nu*t),1 >

cos_nut_scalar_1_from_integral = 0.0
sin_nut_scalar_1_from_integral = 0.0
n = int(2*T/dT)
t = 0
cos_scalar_t = m.cos(nu*t)*weight(t,p)
sin_scalar_t = m.sin(nu*t)*weight(t,p)
for i in range(0,n):
      t = (i+1)*dT
      cos_scalar_t_dT = m.cos(nu*t)*weight(t,p)
      sin_scalar_t_dT = m.sin(nu*t)*weight(t,p)
      cos_nut_scalar_1_from_integral = cos_nut_scalar_1_from_integral + 0.5*dT*(cos_scalar_t + cos_scalar_t_dT)
      sin_nut_scalar_1_from_integral = sin_nut_scalar_1_from_integral + 0.5*dT*(sin_scalar_t + sin_scalar_t_dT)
      cos_scalar_t = cos_scalar_t_dT
      sin_scalar_t = sin_scalar_t_dT

dt = 2*T-n*dT
cos_scalar_t_dT = m.cos(nu*2*T)*weight(2*T,p)
sin_scalar_t_dT = m.sin(nu*2*T)*weight(2*T,p)
cos_nut_scalar_1_from_integral = cos_nut_scalar_1_from_integral + 0.5*dt*(cos_scalar_t + cos_scalar_t_dT)
sin_nut_scalar_1_from_integral = sin_nut_scalar_1_from_integral + 0.5*dt*(sin_scalar_t + sin_scalar_t_dT)
cos_nut_scalar_1_from_integral = cos_nut_scalar_1_from_integral/(2*T)
sin_nut_scalar_1_from_integral = sin_nut_scalar_1_from_integral/(2*T)
print("cos_nut_scalar_1               =", cos_nut_scalar_1)
print("cos_nut_scalar_1_from_integral =", cos_nut_scalar_1_from_integral)
print("difference                     =", abs(cos_nut_scalar_1_from_integral-cos_nut_scalar_1))
print("sin_nut_scalar_1               =", sin_nut_scalar_1)
print("sin_nut_scalar_1_from_integral =", sin_nut_scalar_1_from_integral)
print("difference                     =", abs(sin_nut_scalar_1_from_integral-sin_nut_scalar_1))

############################################
######## Computing the average of f ########
############################################

#Naively analytically:
A_naive_analytical = a_0+1.0/(2.0*T)*(2.0*a_1/nu_1*np.sin(2.0*nu_1*T)+2.0*a_2/nu_2*np.cos(2.0*nu_2*T)+2.0*a_3/nu_3*np.sin(2.0*nu_3*T)-2.0*a_2/nu_2)
print("A_naive_analytical = ", A_naive_analytical)

#Naively numerically:
A_naive_numerical = 0.0
n = int(2*T/dT)
t = 0
ft = f(t)
for i in range(0,n):
      t = (i+1)*dT
      ft_dT = f(t)
      A_naive_numerical = A_naive_numerical + 0.5*dT*(ft + ft_dT)
      ft = ft_dT

dt = 2*T-n*dT
ft_dT = f(2*T)
A_naive_numerical = A_naive_numerical + 0.5*dt*(ft + ft_dT)
A_naive_numerical = A_naive_numerical/(2*T)
print("A_naive_numerical =", A_naive_numerical)
print("difference        =", abs(A_naive_analytical-A_naive_numerical))
print("T                 =", T)
print("dT                =", dT)

'''
#########################################################################################
######## Precision of the numerical computation of integrals as a function of dT ########
#########################################################################################
N_points  = 2500
log10_dTT = np.linspace(-3,0,N_points)
dTT       = 10**log10_dTT * min(T1, T2, T3)
trapeze             = []
simpson_3pt         = []
gauss_2pt           = []
gauss_3pt           = []
gauss_4pt           = []
incr = 0
for dT in dTT:
      A_naive_trapeze     = 0.0
      A_naive_simpson_3pt = 0.0
      A_naive_gauss_2pt   = 0.0
      A_naive_gauss_3pt   = 0.0
      A_naive_gauss_4pt   = 0.0
      n  = int(2*T/dT)
      dT = 2*T/n #Slightly changing the value of dT so that 2T/dT is an integer
      ######## Trapeze ########
      t  = 0
      ft = f(t)
      for i in range(0,n):
            t = (i+1)*dT
            ft_dT = f(t)
            A_naive_trapeze = A_naive_trapeze + 0.5*dT*(ft + ft_dT)
            ft = ft_dT
      A_naive_trapeze = A_naive_trapeze/(2*T)
      trapeze.append(A_naive_trapeze)
      ######## Simpson_3pt ########
      t  = 0
      ft = f(t)
      for i in range(0,n):
            t1  = (i+0.5)*dT
            t2  = (i+1)  *dT
            ft1 = f(t1)
            ft2 = f(t2)
            A_naive_simpson_3pt = A_naive_simpson_3pt + dT*(ft + 4*ft1 + ft2)/6
            ft = ft2
      A_naive_simpson_3pt = A_naive_simpson_3pt/(2*T)
      simpson_3pt.append(A_naive_simpson_3pt)
      ######## Gauss_2pt ########
      w1  = 1
      w2  = 1
      xi1 = -1/m.sqrt(3)
      xi2 =  1/m.sqrt(3)
      for i in range(0,n):
            t1  = dT*(i+0.5+xi1/2)
            t2  = dT*(i+0.5+xi2/2)
            ft1 = f(t1)
            ft2 = f(t2)
            A_naive_gauss_2pt = A_naive_gauss_2pt + 0.5*dT*w1*ft1 + 0.5*dT*w2*ft2
      A_naive_gauss_2pt = A_naive_gauss_2pt/(2*T)
      gauss_2pt.append(A_naive_gauss_2pt)
      ######## Gauss_3pt ########
      w1  = 5/9
      w2  = 8/9
      w3  = 5/9
      xi1 = -m.sqrt(3/5)
      xi2 = 0
      xi3 = m.sqrt(3/5)
      for i in range(0,n):
            t1  = dT*(i+0.5+xi1/2)
            t2  = dT*(i+0.5+xi2/2)
            t3  = dT*(i+0.5+xi3/2)
            ft1 = f(t1)
            ft2 = f(t2)
            ft3 = f(t3)
            A_naive_gauss_3pt = A_naive_gauss_3pt + 0.5*dT*w1*ft1 + 0.5*dT*w2*ft2 + 0.5*dT*w3*ft3
      A_naive_gauss_3pt = A_naive_gauss_3pt/(2*T)
      gauss_3pt.append(A_naive_gauss_3pt)
      ######## Gauss_4pt ########
      w1  = (18+m.sqrt(30))/36
      w2  = (18+m.sqrt(30))/36
      w3  = (18-m.sqrt(30))/36
      w4  = (18-m.sqrt(30))/36
      xi1 = -m.sqrt(3/7-2/7*m.sqrt(6/5))
      xi2 =  m.sqrt(3/7-2/7*m.sqrt(6/5))
      xi3 = -m.sqrt(3/7+2/7*m.sqrt(6/5))
      xi4 =  m.sqrt(3/7+2/7*m.sqrt(6/5))
      for i in range(0,n):
            t1  = dT*(i+0.5+xi1/2)
            t2  = dT*(i+0.5+xi2/2)
            t3  = dT*(i+0.5+xi3/2)
            t4  = dT*(i+0.5+xi4/2)
            ft1 = f(t1)
            ft2 = f(t2)
            ft3 = f(t3)
            ft4 = f(t4)
            A_naive_gauss_4pt = A_naive_gauss_4pt + 0.5*dT*w1*ft1 + 0.5*dT*w2*ft2 + 0.5*dT*w3*ft3 + 0.5*dT*w4*ft4
      A_naive_gauss_4pt = A_naive_gauss_4pt/(2*T)
      gauss_4pt.append(A_naive_gauss_4pt)
      
      incr = incr + 1
      print("Progress =", incr/N_points*100, " %")
      
trapeze           = np.array(trapeze)
simpson_3pt       = np.array(simpson_3pt)
gauss_2pt         = np.array(gauss_2pt)
gauss_3pt         = np.array(gauss_3pt)
gauss_4pt         = np.array(gauss_4pt)
error_trapeze     = np.log10(abs(trapeze-A_naive_analytical))
error_simpson_3pt = np.log10(abs(simpson_3pt-A_naive_analytical))
error_gauss_2pt   = np.log10(abs(gauss_2pt-A_naive_analytical))
error_gauss_3pt   = np.log10(abs(gauss_3pt-A_naive_analytical))
error_gauss_4pt   = np.log10(abs(gauss_4pt-A_naive_analytical))
dTT               = np.log10(dTT/min(T1, T2, T3))
my_yticks=[r'$10^{-16}$',r'$10^{-14}$',r'$10^{-12}$',r'$10^{-10}$',r'$10^{-8}$',r'$10^{-6}$',r'$10^{-4}$',r'$10^{-2}$',r'$1$']
py.yticks([-16,-14,-12,-10,-8,-6,-4,-2,0],my_yticks)
my_xticks=[r'$10^{-3}$',r'$10^{-2.5}$',r'$10^{-2}$',r'$10^{-1.5}$',r'$10^{-1}$',r'$10^{-0.5}$',r'$1$']
py.xticks([-3,-2.5,-2,-1.5,-1,-0.5,0],my_xticks)
py.plot(dTT, error_trapeze,     "-", color = "red",    linewidth = 5, label="Trapeze")
py.plot(dTT, error_simpson_3pt, "-", color = "blue",   linewidth = 5, label="Simpson 3pt")
py.plot(dTT, error_gauss_2pt,   "-", color = "green",  linewidth = 5, label="Gauss 2pt")
py.plot(dTT, error_gauss_3pt,   "-", color = "orange", linewidth = 5, label="Gauss 3pt")
py.plot(dTT, error_gauss_4pt,   "-", color = "purple", linewidth = 5, label="Gauss 4pt")
py.xticks(fontsize=25)
py.yticks(fontsize=25)
py.xlabel(r"$dT\,/\,\min(T_1,T_2,T_3)$", fontsize=25)
py.ylabel(r"Relative error in the computation of $\frac{1}{2T}\int_0^{2T}f(t)dt$", fontsize=25)
py.grid(linewidth=0.4)
py.legend(fontsize = 25)
py.show()
'''
'''
#######################################################################
######## Plot of Re(<f(t),1>) as a function of p for a given T ########
#######################################################################
N_points = 100001
om = np.linspace(omega_min,omega_max,N_points)

p1_numerical = []
n  = int(2*T/dT)
dT = 2*T/n #Slightly changing the value of dT so that 2T/dT is an integer
w1  = 5/9
w2  = 8/9
w3  = 5/9
xi1 = -m.sqrt(3/5)
xi2 = 0
xi3 = m.sqrt(3/5)
for omeg in om:
      A  = 0
      for i in range(0,n):
            t1  = dT*(i+0.5+xi1/2)
            t2  = dT*(i+0.5+xi2/2)
            t3  = dT*(i+0.5+xi3/2)
            ft1 = f(t1)*m.cos(omeg*t1)*weight(t1,1)
            ft2 = f(t2)*m.cos(omeg*t2)*weight(t2,1)
            ft3 = f(t3)*m.cos(omeg*t3)*weight(t3,1)
            A = A + 0.5*dT*w1*ft1 + 0.5*dT*w2*ft2 + 0.5*dT*w3*ft3
      A = A/(2*T)
      p1_numerical.append(A)
p1_numerical = np.array(p1_numerical)

p0 = 0.5*a_0*(np.cos(om*T)*hp(om*T,0)+np.cos(-om*T)*hp(-om*T,0)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,0)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,0))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,0)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,0)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,0)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,0))
p1 = 0.5*a_0*(np.cos(om*T)*hp(om*T,1)+np.cos(-om*T)*hp(-om*T,1)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,1)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,1))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,1)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,1)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,1)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,1))
p2 = 0.5*a_0*(np.cos(om*T)*hp(om*T,2)+np.cos(-om*T)*hp(-om*T,2)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,2)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,2))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,2)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,2)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,2)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,2))
p3 = 0.5*a_0*(np.cos(om*T)*hp(om*T,3)+np.cos(-om*T)*hp(-om*T,3)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,3)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,3))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,3)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,3)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,3)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,3))
p4 = 0.5*a_0*(np.cos(om*T)*hp(om*T,4)+np.cos(-om*T)*hp(-om*T,4)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,4)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,4))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,4)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,4)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,4)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,4))
p5 = 0.5*a_0*(np.cos(om*T)*hp(om*T,5)+np.cos(-om*T)*hp(-om*T,5)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,5)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,5))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,5)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,5)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,5)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,5))
py.plot(om, p5, "-", color = "red",    linewidth = 3, label=r'$p=5$')
py.plot(om, p4, "-", color = "blue",   linewidth = 3, label=r'$p=4$')
py.plot(om, p3, "-", color = "gold",   linewidth = 3, label=r'$p=3$')
py.plot(om, p2, "-", color = "orange", linewidth = 3, label=r'$p=2$')
py.plot(om, p1, "-", color = "green",  linewidth = 3, label=r'$p=1$')
py.plot(om, p0, "-", color = "purple", linewidth = 3, label=r'$p=0$')
py.vlines(x=0,    ymin=-3, ymax=3, colors='black', ls='-', lw=2)
py.vlines(x=nu_1, ymin=-3, ymax=3, colors='black', ls='-', lw=2)
py.vlines(x=nu_2, ymin=-3, ymax=3, colors='black', ls='-', lw=2)
py.vlines(x=nu_3, ymin=-3, ymax=3, colors='black', ls='-', lw=2)
py.xticks(fontsize=25)
py.yticks(fontsize=25)
py.xlabel(r"$\omega$", fontsize=30)
py.ylabel(r"$\mathfrak{R}_e\left<f(t),e^{i\omega t}\right>$", fontsize=30)
py.grid(linewidth=0.4)
py.legend(fontsize = 25)
py.show()
'''

######################################################################################
######## Precision of the computation of the average as a function of T and p ########
######################################################################################
N_points = 25000
log10_TT = np.linspace(0,3,N_points)
TT       = 10**log10_TT * max(T1, T2, T3)
error_0 = []
error_1 = []
error_2 = []
error_3 = []
error_4 = []
error_5 = []
error_6 = []
incr = 0
om = 0
for T in TT:
      ######## p = 0 ########
      error = abs(0.5*a_0*(np.cos(om*T)*hp(om*T,0)+np.cos(-om*T)*hp(-om*T,0)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,0)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,0))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,0)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,0)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,0)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,0))-a_0)
      error_0.append(error)
      error = abs(0.5*a_0*(np.cos(om*T)*hp(om*T,1)+np.cos(-om*T)*hp(-om*T,1)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,1)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,1))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,1)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,1)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,1)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,1))-a_0)
      error_1.append(error)
      error = abs(0.5*a_0*(np.cos(om*T)*hp(om*T,2)+np.cos(-om*T)*hp(-om*T,2)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,2)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,2))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,2)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,2)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,2)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,2))-a_0)
      error_2.append(error)
      error = abs(0.5*a_0*(np.cos(om*T)*hp(om*T,3)+np.cos(-om*T)*hp(-om*T,3)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,3)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,3))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,3)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,3)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,3)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,3))-a_0)
      error_3.append(error)
      error = abs(0.5*a_0*(np.cos(om*T)*hp(om*T,4)+np.cos(-om*T)*hp(-om*T,4)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,4)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,4))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,4)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,4)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,4)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,4))-a_0)
      error_4.append(error)
      error = abs(0.5*a_0*(np.cos(om*T)*hp(om*T,5)+np.cos(-om*T)*hp(-om*T,5)) + a_1*(np.cos((nu_1+om)*T)*hp((nu_1+om)*T,5)+np.cos((nu_1-om)*T)*hp((nu_1-om)*T,5))-a_2*(np.sin((nu_2+om)*T)*hp((nu_2+om)*T,5)+np.sin((nu_2-om)*T)*hp((nu_2-om)*T,5)) + a_3*(np.cos((nu_3+om)*T)*hp((nu_3+om)*T,5)+np.cos((nu_3-om)*T)*hp((nu_3-om)*T,5))-a_0)
      error_5.append(error)
      incr = incr + 1
      print("Progress =", incr/N_points*100, " %")

error_0 = np.log10(np.array(error_0))
error_1 = np.log10(np.array(error_1))
error_2 = np.log10(np.array(error_2))
error_3 = np.log10(np.array(error_3))
error_4 = np.log10(np.array(error_4))
error_5 = np.log10(np.array(error_5))
TT      = np.log10(TT/max(T1, T2, T3))
my_yticks=[r'$10^{-16}$',r'$10^{-14}$',r'$10^{-12}$',r'$10^{-10}$',r'$10^{-8}$',r'$10^{-6}$',r'$10^{-4}$',r'$10^{-2}$',r'$1$']
py.yticks([-16,-14,-12,-10,-8,-6,-4,-2,0],my_yticks)
my_xticks=[r'$1$',r'$10^{0.5}$',r'$10^{1}$',r'$10^{1.5}$',r'$10^{2}$',r'$10^{2.5}$',r'$10^{3}$']
py.xticks([0,0.5,1,1.5,2,2.5,3],my_xticks)
py.plot(TT, error_0, "-", color = "purple", linewidth = 3, label=r"$p=0$")
py.plot(TT, error_1, "-", color = "green",  linewidth = 3, label=r"$p=1$")
py.plot(TT, error_2, "-", color = "orange", linewidth = 3, label=r"$p=2$")
py.plot(TT, error_3, "-", color = "gold",   linewidth = 3, label=r"$p=3$")
py.plot(TT, error_4, "-", color = "blue",   linewidth = 3, label=r"$p=4$")
py.plot(TT, error_5, "-", color = "red",    linewidth = 3, label=r"$p=5$")
py.xticks(fontsize=25)
py.yticks(fontsize=25)
py.xlabel(r"$T\,/\,\max(T_1,T_2,T_3)$", fontsize=30)
py.ylabel(r"$\left|a_0-\left<f(t),1\right>\right|$", fontsize=30)
py.grid(linewidth=0.4)
py.legend(fontsize = 25)
py.show()

'''
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
TT = np.linspace(20.0,1000.0,5000)
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
      print("progress = ", i/5000*100, "%")
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






            
            
            
