import math as m
import matplotlib.pyplot as py
import matplotlib
import numpy as np

m0    = 1
m1    = 0.0000164
m2    = 0.0000232
G     = 39.478417604357434475337964
mu1   = G*(m0+m1)
mu2   = G*(m0+m2)
beta1 = m0*m1/(m0+m1)
beta2 = m0*m2/(m0+m2)
a1    = 1
a2    = 1.2114137
e1    = 0.078
e2    = 0.15
lbd_1 = -1.2678
lbd_2 = 0.27182
vrp_1 = -1.0614
vrp_2 = 0.81578
Lbd_1 = beta1*m.sqrt(mu1*a1)
Lbd_2 = beta2*m.sqrt(mu2*a2)
Lbd_10= m1*m.sqrt(G*m0*a1)
Lbd_20= m2*m.sqrt(G*m0*a2)
D1    = Lbd_1*(1-m.sqrt(1-e1**2))
D2    = Lbd_2*(1-m.sqrt(1-e2**2))
sq2D_1= m.sqrt(2*D1/Lbd_10)
sq2D_2= m.sqrt(2*D2/Lbd_20)
K_12  = G*m1*m2/a2
C_1   =  2.8404322879
C_2   = -3.2832571442
phi_1 = lbd_1-lbd_2
phi_2 = -3*lbd_1+4*lbd_2
sig_1 = phi_2-vrp_1
sig_2 = phi_2-vrp_2
u_1   = sq2D_1*m.cos(sig_1)
u_2   = sq2D_2*m.cos(sig_2)
v_1   = sq2D_1*m.sin(sig_1)
v_2   = sq2D_2*m.sin(sig_2)
Phi_1 = 4*Lbd_1+3*Lbd_2
Phi_2 = Lbd_1-D1+Lbd_2-D2

dHkdLbd1 = mu1**2*beta1**3/Lbd_1**3
dHkdLbd2 = mu2**2*beta2**3/Lbd_2**3
dHkdLbd3 = mu3**2*beta3**3/Lbd_3**3
dHpdlbd1 = K_13*(-(C_4*sq2D_3*m.sin(vrp_3-2*lbd_3+lbd_1))-C_3*sq2D_1*m.sin(vrp_1-2*lbd_3+lbd_1))+K_12*(-(3*C_2*sq2D_2*m.sin(vrp_2-4*lbd_2+3*lbd_1))-3*C_1*sq2D_1*m.sin(vrp_1-4*lbd_2+3*lbd_1))
dHpdlbd2 = K_23*(-(2*C_6*sq2D_3*m.sin(vrp_3-3*lbd_3+2*lbd_2))-2*C_5*sq2D_2*m.sin(vrp_2-3*lbd_3+2*lbd_2))+K_12*(4*C_2*sq2D_2*m.sin(vrp_2-4*lbd_2+3*lbd_1)+4*C_1*sq2D_1*m.sin(vrp_1-4*lbd_2+3*lbd_1))
dHpdlbd3 = K_13*(2*C_4*sq2D_3*m.sin(vrp_3-2*lbd_3+lbd_1)+2*C_3*sq2D_1*m.sin(vrp_1-2*lbd_3+lbd_1))+K_23*(3*C_6*sq2D_3*m.sin(vrp_3-3*lbd_3+2*lbd_2)+3*C_5*sq2D_2*m.sin(vrp_2-3*lbd_3+2*lbd_2))
dHpdvrp1 = C_3*K_13*sq2D_1*m.sin(vrp_1-2*lbd_3+lbd_1)+C_1*K_12*sq2D_1*m.sin(vrp_1-4*lbd_2+3*lbd_1)
dHpdvrp2 = C_5*K_23*sq2D_2*m.sin(vrp_2-3*lbd_3+2*lbd_2)+C_2*K_12*sq2D_2*m.sin(vrp_2-4*lbd_2+3*lbd_1)
dHpdvrp3 = C_4*K_13*sq2D_3*m.sin(vrp_3-2*lbd_3+lbd_1)+C_6*K_23*sq2D_3*m.sin(vrp_3-3*lbd_3+2*lbd_2)
dHpdsq2D1= C_3*K_13*m.cos(vrp_1-2*lbd_3+lbd_1)+C_1*K_12*m.cos(vrp_1-4*lbd_2+3*lbd_1)
dHpdsq2D2= C_5*K_23*m.cos(vrp_2-3*lbd_3+2*lbd_2)+C_2*K_12*m.cos(vrp_2-4*lbd_2+3*lbd_1)
dHpdsq2D3= C_4*K_13*m.cos(vrp_3-2*lbd_3+lbd_1)+C_6*K_23*m.cos(vrp_3-3*lbd_3+2*lbd_2)

print("dHkdLbd1 = ", dHkdLbd1)
print("dHkdLbd2 = ", dHkdLbd2)
print("dHkdLbd3 = ", dHkdLbd3)
print("dHpdlbd1 = ", dHpdlbd1)
print("dHpdlbd2 = ", dHpdlbd2)
print("dHpdlbd3 = ", dHpdlbd3)
print("dHpdvrp1 = ", dHpdvrp1)
print("dHpdvrp2 = ", dHpdvrp2)
print("dHpdvrp3 = ", dHpdvrp3)
print("dHpdsq2D1= ", dHpdsq2D1)
print("dHpdsq2D2= ", dHpdsq2D2)
print("dHpdsq2D3= ", dHpdsq2D3)

dHkdPhi1 = (beta3**3*mu3**2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*beta2**3*mu2**2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3+(beta1**3*mu1**2)/Phi_1**3
dHkdPhi2 = (beta2**3*mu2**2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3-(beta3**3*mu3**2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3
dHkdPhi3 = (3*beta3**3*mu3**2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*beta2**3*mu2**2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
dHkdu1   = (3*Lbd_10*beta3**3*mu3**2*u_1)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_10*beta2**3*mu2**2*u_1)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
dHkdu2   = (3*Lbd_20*beta3**3*mu3**2*u_2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_20*beta2**3*mu2**2*u_2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
dHkdu3   = (3*Lbd_30*beta3**3*mu3**2*u_3)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_30*beta2**3*mu2**2*u_3)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
dHkdv1   = (3*Lbd_10*beta3**3*mu3**2*v_1)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_10*beta2**3*mu2**2*v_1)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
dHkdv2   = (3*Lbd_20*beta3**3*mu3**2*v_2)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_20*beta2**3*mu2**2*v_2)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
dHkdv3   = (3*Lbd_30*beta3**3*mu3**2*v_3)/(3*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2)+3*Phi_3-Phi_2+Phi_1)**3-(2*Lbd_30*beta2**3*mu2**2*v_3)/(-(2*((Lbd_30*(v_3**2+u_3**2))/2+(Lbd_20*(v_2**2+u_2**2))/2+(Lbd_10*(v_1**2+u_1**2))/2))-2*Phi_3+Phi_2-2*Phi_1)**3
dHpdphi1 = K_13*(m.cos(phi_1)*(C_4*v_3+C_3*v_1)-m.sin(phi_1)*(C_4*u_3+C_3*u_1))+K_12*(3*m.cos(3*phi_1)*(C_2*v_2+C_1*v_1)-3*m.sin(3*phi_1)*(C_2*u_2+C_1*u_1))
dHpdu1   = C_1*K_12*m.cos(3*phi_1)+C_3*K_13*m.cos(phi_1)
dHpdu2   = C_2*K_12*m.cos(3*phi_1)+C_5*K_23
dHpdu3   = C_4*K_13*m.cos(phi_1)+C_6*K_23
dHpdv1   = C_1*K_12*m.sin(3*phi_1)+C_3*K_13*m.sin(phi_1)
dHpdv2   = C_2*K_12*m.sin(3*phi_1)
dHpdv3   = C_4*K_13*m.sin(phi_1)

print("\n")
print("dHkdPhi1 = ", dHkdPhi1)
print("dHkdPhi2 = ", dHkdPhi2)
print("dHkdPhi3 = ", dHkdPhi3)
print("dHkdu1   = ", dHkdu1)
print("dHkdu2   = ", dHkdu2)
print("dHkdu3   = ", dHkdu3)
print("dHkdv1   = ", dHkdv1)
print("dHkdv2   = ", dHkdv2)
print("dHkdv3   = ", dHkdv3)
print("dHpdphi1 = ", dHpdphi1)
print("dHpdu1   = ", dHpdu1)
print("dHpdu2   = ", dHpdu2)
print("dHpdu3   = ", dHpdu3)
print("dHpdv1   = ", dHpdv1)
print("dHpdv2   = ", dHpdv2)
print("dHpdv3   = ", dHpdv3)






