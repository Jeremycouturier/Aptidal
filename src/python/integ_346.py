import math as m
import matplotlib.pyplot as py
import matplotlib
import numpy as np

tau   = 0.0025
T     = 3000
os    = 10

m0    = 1
m1    = 0.0000164
m2    = 0.0000232
m3    = 0.0000191
G     = 39.478417604357434475337964
mu1   = G*(m0+m1)
mu2   = G*(m0+m2)
mu3   = G*(m0+m3)
beta1 = m0*m1/(m0+m1)
beta2 = m0*m2/(m0+m2)
beta3 = m0*m3/(m0+m3)
a1    = 1
a2    = 1.2114137
a3    = 1.5874011
e1    = 0.078
e2    = 0.15
e3    = 0.037
lbd_1 = -1.2678
lbd_2 = 0.27182
lbd_3 = -2.8035
vrp_1 = -1.0614
vrp_2 = 0.81578
vrp_3 = -3.07126
Lbd_1 = beta1*m.sqrt(mu1*a1)
Lbd_2 = beta2*m.sqrt(mu2*a2)
Lbd_3 = beta3*m.sqrt(mu3*a3)
Lbd_10= m1*m.sqrt(G*m0*a1)
Lbd_20= m2*m.sqrt(G*m0*a2)
Lbd_30= m3*m.sqrt(G*m0*a3)
D_1   = Lbd_1*(1-m.sqrt(1-e1**2))
D_2   = Lbd_2*(1-m.sqrt(1-e2**2))
D_3   = Lbd_3*(1-m.sqrt(1-e3**2))
sq2D_1= m.sqrt(2*D_1/Lbd_10)
sq2D_2= m.sqrt(2*D_2/Lbd_20)
sq2D_3= m.sqrt(2*D_3/Lbd_30)
g_1   = -vrp_1
g_2   = -vrp_2
g_3   = -vrp_3
K_12  = G*m1*m2/a2
K_13  = G*m1*m3/a3
K_23  = G*m2*m3/a3
C_1   =  2.8404322879410663916
C_2   = -3.2832571441534739165
C_3   =  1.1904935981362627651
C_4   = -0.4283690611689001670
C_5   =  2.0252221881953218485
C_6   = -2.4840046993671887066
phi_1 = lbd_1-2*lbd_2+lbd_3
phi_2 = lbd_2-lbd_3
phi_3 = -2*lbd_2+3*lbd_3
sig_1 = phi_3-vrp_1
sig_2 = phi_3-vrp_2
sig_3 = phi_3-vrp_3
u_1   = sq2D_1*m.cos(sig_1)
u_2   = sq2D_2*m.cos(sig_2)
u_3   = sq2D_3*m.cos(sig_3)
v_1   = sq2D_1*m.sin(sig_1)
v_2   = sq2D_2*m.sin(sig_2)
v_3   = sq2D_3*m.sin(sig_3)
Phi_1 = Lbd_1
Phi_2 = 4*Lbd_1+3*Lbd_2+2*Lbd_3
Phi_3 = (Lbd_1-D_1)+(Lbd_2-D_2)+(Lbd_3-D_3)
eta   = m.sqrt(mu1/a1**3)

def old2new(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3):
      sq2D_1= m.sqrt(2*D_1/Lbd_10)
      sq2D_2= m.sqrt(2*D_2/Lbd_20)
      sq2D_3= m.sqrt(2*D_3/Lbd_30)
      vrp_1 = -g_1
      vrp_2 = -g_2
      vrp_3 = -g_3
      phi_1 = lbd_1-2*lbd_2+lbd_3
      phi_2 = lbd_2-lbd_3
      phi_3 = -2*lbd_2+3*lbd_3
      sig_1 = phi_3-vrp_1
      sig_2 = phi_3-vrp_2
      sig_3 = phi_3-vrp_3
      Phi_1 = Lbd_1
      Phi_2 = 4*Lbd_1+3*Lbd_2+2*Lbd_3
      Phi_3 = (Lbd_1-D_1)+(Lbd_2-D_2)+(Lbd_3-D_3)
      u_1   = sq2D_1*m.cos(sig_1)
      u_2   = sq2D_2*m.cos(sig_2)
      u_3   = sq2D_3*m.cos(sig_3)
      v_1   = sq2D_1*m.sin(sig_1)
      v_2   = sq2D_2*m.sin(sig_2)
      v_3   = sq2D_3*m.sin(sig_3)
      return [phi_1, v_1, v_2, v_3, Phi_1, u_1, u_2, u_3]
      

def gradH(phi_1, v_1, v_2, v_3, Phi_1, u_1, u_2, u_3):
      dHkdPhi1 = (beta3**3*mu3**2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*beta2**3*mu2**2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3+(beta1**3*mu1**2)/Phi_1**3
      dHpdphi1 = K_13*(m.cos(phi_1)*(C_4*v_3+C_3*v_1)-m.sin(phi_1)*(C_4*u_3+C_3*u_1))+K_12*(3*m.cos(3*phi_1)*(C_2*v_2+C_1*v_1)-3*m.sin(3*phi_1)*(C_2*u_2+C_1*u_1))
      dHkdu1   = (3*Lbd_10*beta3**3*mu3**2*u_1)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_10*beta2**3*mu2**2*u_1)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
      dHkdu2   = (3*Lbd_20*beta3**3*mu3**2*u_2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_20*beta2**3*mu2**2*u_2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
      dHkdu3   = (3*Lbd_30*beta3**3*mu3**2*u_3)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_30*beta2**3*mu2**2*u_3)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
      dHkdv1   = (3*Lbd_10*beta3**3*mu3**2*v_1)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_10*beta2**3*mu2**2*v_1)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
      dHkdv2   = (3*Lbd_20*beta3**3*mu3**2*v_2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_20*beta2**3*mu2**2*v_2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
      dHkdv3   = (3*Lbd_30*beta3**3*mu3**2*v_3)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_30*beta2**3*mu2**2*v_3)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
      dHpdu1   = C_1*K_12*m.cos(3*phi_1)+C_3*K_13*m.cos(phi_1)
      dHpdu2   = C_2*K_12*m.cos(3*phi_1)+C_5*K_23
      dHpdu3   = C_4*K_13*m.cos(phi_1)+C_6*K_23
      dHpdv1   = C_1*K_12*m.sin(3*phi_1)+C_3*K_13*m.sin(phi_1)
      dHpdv2   = C_2*K_12*m.sin(3*phi_1)
      dHpdv3   = C_4*K_13*m.sin(phi_1)
      return [dHpdphi1, dHkdv1+dHpdv1, dHkdv2+dHpdv2, dHkdv3+dHpdv3, dHkdPhi1, dHkdu1+dHpdu1, dHkdu2+dHpdu2, dHkdu3+dHpdu3]
      
def gradH_original(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3):
      sq2D_1   = m.sqrt(2*D_1/Lbd_10)
      sq2D_2   = m.sqrt(2*D_2/Lbd_20)
      sq2D_3   = m.sqrt(2*D_3/Lbd_30)
      vrp_1    = -g_1
      vrp_2    = -g_2
      vrp_3    = -g_3
      dHkdLbd1 = mu1**2*beta1**3/Lbd_1**3
      dHkdLbd2 = mu2**2*beta2**3/Lbd_2**3
      dHkdLbd3 = mu3**2*beta3**3/Lbd_3**3
      dHpdlbd1 = K_13*(-(C_4*sq2D_3*m.sin(vrp_3-2*lbd_3+lbd_1))-C_3*sq2D_1*m.sin(vrp_1-2*lbd_3+lbd_1))+K_12*(-(3*C_2*sq2D_2*m.sin(vrp_2-4*lbd_2+3*lbd_1))-3*C_1*sq2D_1*m.sin(vrp_1-4*lbd_2+3*lbd_1))
      dHpdlbd2 = K_23*(-(2*C_6*sq2D_3*m.sin(vrp_3-3*lbd_3+2*lbd_2))-2*C_5*sq2D_2*m.sin(vrp_2-3*lbd_3+2*lbd_2))+K_12*(4*C_2*sq2D_2*m.sin(vrp_2-4*lbd_2+3*lbd_1)+4*C_1*sq2D_1*m.sin(vrp_1-4*lbd_2+3*lbd_1))
      dHpdlbd3 = K_13*(2*C_4*sq2D_3*m.sin(vrp_3-2*lbd_3+lbd_1)+2*C_3*sq2D_1*m.sin(vrp_1-2*lbd_3+lbd_1))+K_23*(3*C_6*sq2D_3*m.sin(vrp_3-3*lbd_3+2*lbd_2)+3*C_5*sq2D_2*m.sin(vrp_2-3*lbd_3+2*lbd_2))
      dHpdg1   = C_3*K_13*sq2D_1*m.sin(vrp_1-2*lbd_3+lbd_1)+C_1*K_12*sq2D_1*m.sin(vrp_1-4*lbd_2+3*lbd_1)
      dHpdg2   = C_5*K_23*sq2D_2*m.sin(vrp_2-3*lbd_3+2*lbd_2)+C_2*K_12*sq2D_2*m.sin(vrp_2-4*lbd_2+3*lbd_1)
      dHpdg3   = C_4*K_13*sq2D_3*m.sin(vrp_3-2*lbd_3+lbd_1)+C_6*K_23*sq2D_3*m.sin(vrp_3-3*lbd_3+2*lbd_2)
      dHpdsq2D1= C_3*K_13*m.cos(vrp_1-2*lbd_3+lbd_1)+C_1*K_12*m.cos(vrp_1-4*lbd_2+3*lbd_1)
      dHpdsq2D2= C_5*K_23*m.cos(vrp_2-3*lbd_3+2*lbd_2)+C_2*K_12*m.cos(vrp_2-4*lbd_2+3*lbd_1)
      dHpdsq2D3= C_4*K_13*m.cos(vrp_3-2*lbd_3+lbd_1)+C_6*K_23*m.cos(vrp_3-3*lbd_3+2*lbd_2)
      dHpdD1   = dHpdsq2D1/m.sqrt(2*Lbd_10*D_1)
      dHpdD2   = dHpdsq2D2/m.sqrt(2*Lbd_20*D_2)
      dHpdD3   = dHpdsq2D3/m.sqrt(2*Lbd_30*D_3)
      return [dHpdlbd1, dHpdlbd2, dHpdlbd3, dHpdg1, dHpdg2, dHpdg3, dHkdLbd1, dHkdLbd2, dHkdLbd3, dHpdD1, dHpdD2, dHpdD3]
      
def gradH_expanded(phi_1, v_1, v_2, v_3, Phi_1, u_1, u_2, u_3):
      dHkdPhi1 = eta*(-((9*Lbd_30*v_3**2)/(2*Lbd_20))-(9*v_3**2)/4-(9*Lbd_20*v_2**2)/(4*Lbd_30)-(9*v_2**2)/2-(9*Lbd_10*v_1**2)/(4*Lbd_30)-(9*Lbd_10*v_1**2)/(2*Lbd_20)-(9*Lbd_30*u_3**2)/(2*Lbd_20)-(9*u_3**2)/4-(9*Lbd_20*u_2**2)/(4*Lbd_30)-(9*u_2**2)/2-(9*Lbd_10*u_1**2)/(4*Lbd_30)-(9*Lbd_10*u_1**2)/(2*Lbd_20)-(9*Phi_3)/(2*Lbd_30)-(9*Phi_3)/Lbd_20+(3*Phi_2)/(2*Lbd_30)+(9*Phi_2)/(2*Lbd_20)-(3*Phi_1)/(2*Lbd_30)-(9*Phi_1)/Lbd_20-(3*Phi_1)/Lbd_10)
      dHkdu1   = eta*(-((9*Lbd_10*Lbd_30*u_1*v_3**2)/(2*Lbd_20))-(27*Lbd_10*u_1*v_3**2)/4-(27*Lbd_10*Lbd_20*u_1*v_2**2)/(4*Lbd_30)-(9*Lbd_10*u_1*v_2**2)/2-(27*Lbd_10**2*u_1*v_1**2)/(4*Lbd_30)-(9*Lbd_10**2*u_1*v_1**2)/(2*Lbd_20)-(9*Lbd_10*Lbd_30*u_1*u_3**2)/(2*Lbd_20)-(27*Lbd_10*u_1*u_3**2)/4-(27*Lbd_10*Lbd_20*u_1*u_2**2)/(4*Lbd_30)-(9*Lbd_10*u_1*u_2**2)/2-(27*Lbd_10**2*u_1**3)/(4*Lbd_30)-(9*Lbd_10**2*u_1**3)/(2*Lbd_20)-(27*Lbd_10*Phi_3*u_1)/(2*Lbd_30)-(9*Lbd_10*Phi_3*u_1)/Lbd_20+(9*Lbd_10*Phi_2*u_1)/(2*Lbd_30)+(9*Lbd_10*Phi_2*u_1)/(2*Lbd_20)-(9*Lbd_10*Phi_1*u_1)/(2*Lbd_30)-(9*Lbd_10*Phi_1*u_1)/Lbd_20)
      dHkdu2   = eta*(-((9*Lbd_30*u_2*v_3**2)/2)-(27*Lbd_20*u_2*v_3**2)/4-(27*Lbd_20**2*u_2*v_2**2)/(4*Lbd_30)-(9*Lbd_20*u_2*v_2**2)/2-(27*Lbd_10*Lbd_20*u_2*v_1**2)/(4*Lbd_30)-(9*Lbd_10*u_2*v_1**2)/2-(9*Lbd_30*u_2*u_3**2)/2-(27*Lbd_20*u_2*u_3**2)/4-(27*Lbd_20**2*u_2**3)/(4*Lbd_30)-(9*Lbd_20*u_2**3)/2-(27*Lbd_10*Lbd_20*u_1**2*u_2)/(4*Lbd_30)-(9*Lbd_10*u_1**2*u_2)/2-(27*Lbd_20*Phi_3*u_2)/(2*Lbd_30)-9*Phi_3*u_2+(9*Lbd_20*Phi_2*u_2)/(2*Lbd_30)+(9*Phi_2*u_2)/2-(9*Lbd_20*Phi_1*u_2)/(2*Lbd_30)-9*Phi_1*u_2)
      dHkdu3   = eta*(-((9*Lbd_30**2*u_3*v_3**2)/(2*Lbd_20))-(27*Lbd_30*u_3*v_3**2)/4-(9*Lbd_30*u_3*v_2**2)/2-(27*Lbd_20*u_3*v_2**2)/4-(9*Lbd_10*Lbd_30*u_3*v_1**2)/(2*Lbd_20)-(27*Lbd_10*u_3*v_1**2)/4-(9*Lbd_30**2*u_3**3)/(2*Lbd_20)-(27*Lbd_30*u_3**3)/4-(9*Lbd_30*u_2**2*u_3)/2-(27*Lbd_20*u_2**2*u_3)/4-(9*Lbd_10*Lbd_30*u_1**2*u_3)/(2*Lbd_20)-(27*Lbd_10*u_1**2*u_3)/4-(9*Lbd_30*Phi_3*u_3)/Lbd_20-(27*Phi_3*u_3)/2+(9*Lbd_30*Phi_2*u_3)/(2*Lbd_20)+(9*Phi_2*u_3)/2-(9*Lbd_30*Phi_1*u_3)/Lbd_20-(9*Phi_1*u_3)/2)
      dHkdv1   = eta*(-((9*Lbd_10*Lbd_30*v_1*v_3**2)/(2*Lbd_20))-(27*Lbd_10*v_1*v_3**2)/4-(27*Lbd_10*Lbd_20*v_1*v_2**2)/(4*Lbd_30)-(9*Lbd_10*v_1*v_2**2)/2-(27*Lbd_10**2*v_1**3)/(4*Lbd_30)-(9*Lbd_10**2*v_1**3)/(2*Lbd_20)-(9*Lbd_10*Lbd_30*u_3**2*v_1)/(2*Lbd_20)-(27*Lbd_10*u_3**2*v_1)/4-(27*Lbd_10*Lbd_20*u_2**2*v_1)/(4*Lbd_30)-(9*Lbd_10*u_2**2*v_1)/2-(27*Lbd_10**2*u_1**2*v_1)/(4*Lbd_30)-(9*Lbd_10**2*u_1**2*v_1)/(2*Lbd_20)-(27*Lbd_10*Phi_3*v_1)/(2*Lbd_30)-(9*Lbd_10*Phi_3*v_1)/Lbd_20+(9*Lbd_10*Phi_2*v_1)/(2*Lbd_30)+(9*Lbd_10*Phi_2*v_1)/(2*Lbd_20)-(9*Lbd_10*Phi_1*v_1)/(2*Lbd_30)-(9*Lbd_10*Phi_1*v_1)/Lbd_20)
      dHkdv2   = eta*(-((9*Lbd_30*v_2*v_3**2)/2)-(27*Lbd_20*v_2*v_3**2)/4-(27*Lbd_20**2*v_2**3)/(4*Lbd_30)-(9*Lbd_20*v_2**3)/2-(27*Lbd_10*Lbd_20*v_1**2*v_2)/(4*Lbd_30)-(9*Lbd_10*v_1**2*v_2)/2-(9*Lbd_30*u_3**2*v_2)/2-(27*Lbd_20*u_3**2*v_2)/4-(27*Lbd_20**2*u_2**2*v_2)/(4*Lbd_30)-(9*Lbd_20*u_2**2*v_2)/2-(27*Lbd_10*Lbd_20*u_1**2*v_2)/(4*Lbd_30)-(9*Lbd_10*u_1**2*v_2)/2-(27*Lbd_20*Phi_3*v_2)/(2*Lbd_30)-9*Phi_3*v_2+(9*Lbd_20*Phi_2*v_2)/(2*Lbd_30)+(9*Phi_2*v_2)/2-(9*Lbd_20*Phi_1*v_2)/(2*Lbd_30)-9*Phi_1*v_2)
      dHkdv3   = eta*(-((9*Lbd_30**2*v_3**3)/(2*Lbd_20))-(27*Lbd_30*v_3**3)/4-(9*Lbd_30*v_2**2*v_3)/2-(27*Lbd_20*v_2**2*v_3)/4-(9*Lbd_10*Lbd_30*v_1**2*v_3)/(2*Lbd_20)-(27*Lbd_10*v_1**2*v_3)/4-(9*Lbd_30**2*u_3**2*v_3)/(2*Lbd_20)-(27*Lbd_30*u_3**2*v_3)/4-(9*Lbd_30*u_2**2*v_3)/2-(27*Lbd_20*u_2**2*v_3)/4-(9*Lbd_10*Lbd_30*u_1**2*v_3)/(2*Lbd_20)-(27*Lbd_10*u_1**2*v_3)/4-(9*Lbd_30*Phi_3*v_3)/Lbd_20-(27*Phi_3*v_3)/2+(9*Lbd_30*Phi_2*v_3)/(2*Lbd_20)+(9*Phi_2*v_3)/2-(9*Lbd_30*Phi_1*v_3)/Lbd_20-(9*Phi_1*v_3)/2)
      dHpdphi1 = K_13*(m.cos(phi_1)*(C_4*v_3+C_3*v_1)-m.sin(phi_1)*(C_4*u_3+C_3*u_1))+K_12*(3*m.cos(3*phi_1)*(C_2*v_2+C_1*v_1)-3*m.sin(3*phi_1)*(C_2*u_2+C_1*u_1))
      dHpdu1   = C_1*K_12*m.cos(3*phi_1)+C_3*K_13*m.cos(phi_1)
      dHpdu2   = C_2*K_12*m.cos(3*phi_1)+C_5*K_23
      dHpdu3   = C_4*K_13*m.cos(phi_1)+C_6*K_23
      dHpdv1   = C_1*K_12*m.sin(3*phi_1)+C_3*K_13*m.sin(phi_1)
      dHpdv2   = C_2*K_12*m.sin(3*phi_1)
      dHpdv3   = C_4*K_13*m.sin(phi_1)
      return [dHpdphi1, dHkdv1+dHpdv1, dHkdv2+dHpdv2, dHkdv3+dHpdv3, dHkdPhi1, dHkdu1+dHpdu1, dHkdu2+dHpdu2, dHkdu3+dHpdu3]

def timestep(phi_1, v_1, v_2, v_3, Phi_1, u_1, u_2, u_3):
      tT = 2./3.
      k1 = gradH(phi_1, v_1, v_2, v_3, Phi_1, u_1, u_2, u_3)
      k2 = gradH(phi_1+tT*tau*k1[4], v_1+tT*tau*k1[5]/Lbd_10, v_2+tT*tau*k1[6]/Lbd_20, v_3+tT*tau*k1[7]/Lbd_30, Phi_1-tT*tau*k1[0], u_1-tT*tau*k1[1]/Lbd_10, u_2-tT*tau*k1[2]/Lbd_20, u_3-tT*tau*k1[3]/Lbd_30)
      return [phi_1+0.25*tau*(k1[4]+3*k2[4]), v_1+0.25*tau*(k1[5]+3*k2[5])/Lbd_10, v_2+0.25*tau*(k1[6]+3*k2[6])/Lbd_20, v_3+0.25*tau*(k1[7]+3*k2[7])/Lbd_30, Phi_1-0.25*tau*(k1[0]+3*k2[0]), u_1-0.25*tau*(k1[1]+3*k2[1])/Lbd_10, u_2-0.25*tau*(k1[2]+3*k2[2])/Lbd_20, u_3-0.25*tau*(k1[3]+3*k2[3])/Lbd_30]      

def timestep_original(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3):
      tT = 2./3.
      k1 = gradH_original(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      k2 = gradH_original(lbd_1+tT*tau*k1[6], lbd_2+tT*tau*k1[7], lbd_3+tT*tau*k1[8], g_1+tT*tau*k1[9], g_2+tT*tau*k1[10], g_3+tT*tau*k1[11], Lbd_1-tT*tau*k1[0], Lbd_2-tT*tau*k1[1], Lbd_3-tT*tau*k1[2], D_1-tT*tau*k1[3], D_2-tT*tau*k1[4], D_3-tT*tau*k1[5])
      return [lbd_1+0.25*tau*(k1[6]+3*k2[6]), lbd_2+0.25*tau*(k1[7]+3*k2[7]), lbd_3+0.25*tau*(k1[8]+3*k2[8]), g_1+0.25*tau*(k1[9]+3*k2[9]), g_2+0.25*tau*(k1[10]+3*k2[10]), g_3+0.25*tau*(k1[11]+3*k2[11]), Lbd_1-0.25*tau*(k1[0]+3*k2[0]), Lbd_2-0.25*tau*(k1[1]+3*k2[1]), Lbd_3-0.25*tau*(k1[2]+3*k2[2]), D_1-0.25*tau*(k1[3]+3*k2[3]), D_2-0.25*tau*(k1[4]+3*k2[4]), D_3-0.25*tau*(k1[5]+3*k2[5])]

    
N_step = T/tau
N_step = int(N_step)

file = open("/home/jeremy/Documents/Aptidal_simulation/test/RK2_uv.txt", "w")
for i in range (N_step):
      
      if (i%os == 0):
            file.write(str(tau*i)+" "+str(phi_1)+" "+str(v_1)+" "+str(v_2)+" "+str(v_3)+" "+str(Phi_1)+" "+str(u_1)+" "+str(u_2)+" "+str(u_3)+"\n")
      tmsp  = timestep(phi_1, v_1, v_2, v_3, Phi_1, u_1, u_2, u_3)
      phi_1 = tmsp[0]
      v_1   = tmsp[1]
      v_2   = tmsp[2]
      v_3   = tmsp[3]
      Phi_1 = tmsp[4]
      u_1   = tmsp[5]
      u_2   = tmsp[6]
      u_3   = tmsp[7]
file.close()

file = open("/home/jeremy/Documents/Aptidal_simulation/test/RK2_old.txt", "w")
for i in range (N_step):
      
      uv = old2new(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      if (i%os == 0):
            file.write(str(tau*i)+" "+str(uv[0])+" "+str(uv[1])+" "+str(uv[2])+" "+str(uv[3])+" "+str(uv[4])+" "+str(uv[5])+" "+str(uv[6])+" "+str(uv[7])+"\n")
      tmsp  = timestep_original(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      lbd_1 = tmsp[0]
      lbd_2 = tmsp[1]
      lbd_3 = tmsp[2]
      g_1   = tmsp[3]
      g_2   = tmsp[4]
      g_3   = tmsp[5]
      Lbd_1 = tmsp[6]
      Lbd_2 = tmsp[7]
      Lbd_3 = tmsp[8]
      D_1   = tmsp[9]
      D_2   = tmsp[10]
      D_3   = tmsp[11]
file.close()























