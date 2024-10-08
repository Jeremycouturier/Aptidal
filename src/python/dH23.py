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
a2    = 1.310371
e1    = 0.078
e2    = 0.15
lbd1  = -1.2678
lbd2  = 0.27182
vrp1  = -1.0614
vrp2  = 0.81578
Lbd1  = beta1*m.sqrt(mu1*a1)
Lbd2  = beta2*m.sqrt(mu2*a2)
Lbd10 = m1*m.sqrt(G*m0*a1)
Lbd20 = m2*m.sqrt(G*m0*a2)
D1    = Lbd1*(1-m.sqrt(1-e1**2))
D2    = Lbd2*(1-m.sqrt(1-e2**2))
sq2D1 = m.sqrt(2*D1/Lbd10)
sq2D2 = m.sqrt(2*D2/Lbd20)
K     = G*m1*m2/a2
C_1   = 2.025220535376242
C_2   = -2.484003105205423

dHdLbd1 = mu1**2*beta1**3/Lbd1**3
dHdLbd2 = mu2**2*beta2**3/Lbd2**3
dHdlbd1 = K*(-(2*C_2*sq2D2*m.sin(vrp2-3*lbd2+2*lbd1))-2*C_1*sq2D1*m.sin(vrp1-3*lbd2+2*lbd1))
dHdlbd2 = K*(3*C_2*sq2D2*m.sin(vrp2-3*lbd2+2*lbd1)+3*C_1*sq2D1*m.sin(vrp1-3*lbd2+2*lbd1))
dHdg1   = C_1*K*sq2D1*m.sin(vrp1-3*lbd2+2*lbd1)
dHdg2   = C_2*K*sq2D2*m.sin(vrp2-3*lbd2+2*lbd1)
dHdsq2D1= C_1*K*m.cos(vrp1-3*lbd2+2*lbd1)
dHdsq2D2= C_2*K*m.cos(vrp2-3*lbd2+2*lbd1)

'''
print("dHdLbd1 = ", dHdLbd1)
print("dHdLbd2 = ", dHdLbd2)
print("dHdlbd1 = ", dHdlbd1)
print("dHdlbd2 = ", dHdlbd2)
print("dHdg1   = ", dHdg1)
print("dHdg2   = ", dHdg2)
print("dHdsq2D1= ", dHdsq2D1)
print("dHdsq2D2= ", dHdsq2D2)'''

sig1  = 3*lbd2-2*lbd1-vrp1
sig2  = 3*lbd2-2*lbd1-vrp2
Phi1  = 3*Lbd1+2*Lbd2
Phi2  = Lbd1+Lbd2-D1-D2

dHdsig1 = -(C_1*K*sq2D1*m.sin(sig1))
dHdsig2 = -(C_2*K*sq2D2*m.sin(sig2))
dHdPhi1 = (mu1**2*beta1**3)/(-(2*Phi2)+Phi1-sq2D2**2*Lbd20-sq2D1**2*Lbd10)**3-(mu2**2*beta2**3)/(3*Phi2-Phi1+(3*sq2D2**2*Lbd20)/2+(3*sq2D1**2*Lbd10)/2)**3
dHdPhi2 = (3*mu2**2*beta2**3)/(3*Phi2-Phi1+(3*sq2D2**2*Lbd20)/2+(3*sq2D1**2*Lbd10)/2)**3-(2*mu1**2*beta1**3)/(-(2*Phi2)+Phi1-sq2D2**2*Lbd20-sq2D1**2*Lbd10)**3
dHdsq2D1= C_1*K*m.cos(sig1)+(3*sq2D1*mu2**2*Lbd10*beta2**3)/(3*Phi2-Phi1+(3*sq2D2**2*Lbd20)/2+(3*sq2D1**2*Lbd10)/2)**3-(2*sq2D1*mu1**2*Lbd10*beta1**3)/(-(2*Phi2)+Phi1-sq2D2**2*Lbd20-sq2D1**2*Lbd10)**3
dHdsq2D2= C_2*K*m.cos(sig2)+(3*sq2D2*mu2**2*Lbd20*beta2**3)/(3*Phi2-Phi1+(3*sq2D2**2*Lbd20)/2+(3*sq2D1**2*Lbd10)/2)**3-(2*sq2D2*mu1**2*Lbd20*beta1**3)/(-(2*Phi2)+Phi1-sq2D2**2*Lbd20-sq2D1**2*Lbd10)**3

'''
print("dHdsig1 = ", dHdsig1)
print("dHdsig2 = ", dHdsig2)
print("dHdPhi1 = ", dHdPhi1)
print("dHdPhi2 = ", dHdPhi2)
print("dHdsq2D1= ", dHdsq2D1)
print("dHdsq2D2= ", dHdsq2D2)'''

u1   = sq2D1*m.cos(sig1)
u2   = sq2D2*m.cos(sig2)
v1   = sq2D1*m.sin(sig1)
v2   = sq2D2*m.sin(sig2)

dHKdu1 = (3*u1*mu2**2*Lbd10*beta2**3)/(3*Phi2-Phi1+(3*(v2**2+u2**2)*Lbd20)/2+(3*(v1**2+u1**2)*Lbd10)/2)**3-(2*u1*mu1**2*Lbd10*beta1**3)/(-(2*Phi2)+Phi1-(v2**2+u2**2)*Lbd20-(v1**2+u1**2)*Lbd10)**3
dHKdu2 = (3*u2*mu2**2*Lbd20*beta2**3)/(3*Phi2-Phi1+(3*(v2**2+u2**2)*Lbd20)/2+(3*(v1**2+u1**2)*Lbd10)/2)**3-(2*u2*mu1**2*Lbd20*beta1**3)/(-(2*Phi2)+Phi1-(v2**2+u2**2)*Lbd20-(v1**2+u1**2)*Lbd10)**3
dHKdv1 = (3*v1*mu2**2*Lbd10*beta2**3)/(3*Phi2-Phi1+(3*(v2**2+u2**2)*Lbd20)/2+(3*(v1**2+u1**2)*Lbd10)/2)**3-(2*v1*mu1**2*Lbd10*beta1**3)/(-(2*Phi2)+Phi1-(v2**2+u2**2)*Lbd20-(v1**2+u1**2)*Lbd10)**3
dHKdv2 = (3*v2*mu2**2*Lbd20*beta2**3)/(3*Phi2-Phi1+(3*(v2**2+u2**2)*Lbd20)/2+(3*(v1**2+u1**2)*Lbd10)/2)**3-(2*v2*mu1**2*Lbd20*beta1**3)/(-(2*Phi2)+Phi1-(v2**2+u2**2)*Lbd20-(v1**2+u1**2)*Lbd10)**3

dHPdu1 = C_1*K
dHPdu2 = C_2*K

dHKdPhi1 = (mu1**2*beta1**3)/(-(2*Phi2)+Phi1-sq2D2**2*Lbd20-sq2D1**2*Lbd10)**3-(mu2**2*beta2**3)/(3*Phi2-Phi1+(3*sq2D2**2*Lbd20)/2+(3*sq2D1**2*Lbd10)/2)**3
dHKdPhi2 = (3*mu2**2*beta2**3)/(3*Phi2-Phi1+(3*sq2D2**2*Lbd20)/2+(3*sq2D1**2*Lbd10)/2)**3-(2*mu1**2*beta1**3)/(-(2*Phi2)+Phi1-sq2D2**2*Lbd20-sq2D1**2*Lbd10)**3

print("dHKdPhi1 = ", dHKdPhi1)
print("dHKdPhi2 = ", dHKdPhi2)
print("dHKdu1   = ", dHKdu1)
print("dHKdu2   = ", dHKdu2)
print("dHKdv1   = ", dHKdv1)
print("dHKdv2   = ", dHKdv2)
print("dHPdu1   = ", dHPdu1)
print("dHPdu2   = ", dHPdu2)


























