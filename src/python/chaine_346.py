import math as m
import matplotlib.pyplot as py
import matplotlib
import numpy as np

tau   = 0.0078125
T     = 5000
os    = 20

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
e1    = 0.02#0.078
e2    = 0.02#0.15
e3    = 0.02#0.037
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
      Hk    = -beta1**3*mu1**2/(2*Lbd_1**2)-beta2**3*mu2**2/(2*Lbd_2**2)-beta3**3*mu3**2/(2*Lbd_3**2)
      Hp    = K_12*(C_1*sq2D_1*m.cos(3*lbd_1-4*lbd_2+vrp_1)+C_2*sq2D_2*m.cos(3*lbd_1-4*lbd_2+vrp_2))+K_13*(C_3*sq2D_1*m.cos(lbd_1-2*lbd_3+vrp_1)+C_4*sq2D_3*m.cos(lbd_1-2*lbd_3+vrp_3))+K_23*(C_5*sq2D_2*m.cos(2*lbd_2-3*lbd_3+vrp_2)+C_6*sq2D_3*m.cos(2*lbd_2-3*lbd_3+vrp_3))
      e1    = m.sqrt(1-(1-0.5*Lbd_10/Lbd_1*(u_1**2+v_1**2))**2)
      e2    = m.sqrt(1-(1-0.5*Lbd_20/Lbd_2*(u_2**2+v_2**2))**2)
      e3    = m.sqrt(1-(1-0.5*Lbd_30/Lbd_3*(u_3**2+v_3**2))**2)
      return [phi_1, v_1, v_2, v_3, Phi_1, u_1, u_2, u_3, Phi_2, Phi_3, Hk+Hp, e1, e2, e3]
      
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
            

def exp_tauLA(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3):
      k1 = gradH_original(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      return [lbd_1+0.5*tau*k1[6], lbd_2+0.5*tau*k1[7], lbd_3+0.5*tau*k1[8], g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3]

def exp_tauLB(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3):
      tT = 2./3.
      k1 = gradH_original(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      k2 = gradH_original(lbd_1, lbd_2, lbd_3, g_1+tT*tau*k1[9], g_2+tT*tau*k1[10], g_3+tT*tau*k1[11], Lbd_1-tT*tau*k1[0], Lbd_2-tT*tau*k1[1], Lbd_3-tT*tau*k1[2], D_1-tT*tau*k1[3], D_2-tT*tau*k1[4], D_3-tT*tau*k1[5])
      return [lbd_1, lbd_2, lbd_3, g_1+0.25*tau*(k1[9]+3*k2[9]), g_2+0.25*tau*(k1[10]+3*k2[10]), g_3+0.25*tau*(k1[11]+3*k2[11]), Lbd_1-0.25*tau*(k1[0]+3*k2[0]), Lbd_2-0.25*tau*(k1[1]+3*k2[1]), Lbd_3-0.25*tau*(k1[2]+3*k2[2]), D_1-0.25*tau*(k1[3]+3*k2[3]), D_2-0.25*tau*(k1[4]+3*k2[4]), D_3-0.25*tau*(k1[5]+3*k2[5])]

    
N_step = T/tau
N_step = int(N_step)

file = open("/home/jeremy/Documents/Aptidal_simulation/test/SABA1_python.txt", "w")
uv = old2new(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
file.write(str(0.0)+" "+str(uv[0])+" "+str(uv[1])+" "+str(uv[2])+" "+str(uv[3])+" "+str(uv[4])+" "+str(uv[5])+" "+str(uv[6])+" "+str(uv[7])+" "+str(uv[8])+" "+str(uv[9])+" "+str(uv[10])+" "+str(uv[11])+" "+str(uv[12])+" "+str(uv[13])+"\n")

for i in range (N_step):
      
      tmsp  = exp_tauLA(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      lbd_1 = tmsp[0]
      lbd_2 = tmsp[1]
      lbd_3 = tmsp[2]
      tmsp  = exp_tauLB(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      g_1   = tmsp[3]
      g_2   = tmsp[4]
      g_3   = tmsp[5]
      Lbd_1 = tmsp[6]
      Lbd_2 = tmsp[7]
      Lbd_3 = tmsp[8]
      D_1   = tmsp[9]
      D_2   = tmsp[10]
      D_3   = tmsp[11]
      tmsp  = exp_tauLA(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      lbd_1 = tmsp[0]
      lbd_2 = tmsp[1]
      lbd_3 = tmsp[2]
      uv = old2new(lbd_1, lbd_2, lbd_3, g_1, g_2, g_3, Lbd_1, Lbd_2, Lbd_3, D_1, D_2, D_3)
      if (i%os == 0 and i != 0):
            file.write(str(tau*i)+" "+str(uv[0])+" "+str(uv[1])+" "+str(uv[2])+" "+str(uv[3])+" "+str(uv[4])+" "+str(uv[5])+" "+str(uv[6])+" "+str(uv[7])+" "+str(uv[8])+" "+str(uv[9])+" "+str(uv[10])+" "+str(uv[11])+" "+str(uv[12])+" "+str(uv[13])+"\n")
      
file.close()























