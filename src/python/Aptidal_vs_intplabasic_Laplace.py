import math as m
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np

intplabasic = "/home/jeremy/Documents/Aptidal_simulation/intplabasic/To_be_kept/ABAH1_chain_346_ABAH1.ell" #Change 235 by 112 or 3468 or 346

Aptidal_O1  = "/home/jeremy/Documents/Aptidal_simulation/test/To_be_kept/SABA1_346_O1.txt" #Change 235 by 112 or 3468 or 346
Aptidal_O2  = "/home/jeremy/Documents/Aptidal_simulation/test/To_be_kept/SABA1_346_O2.txt" #Change 235 by 112 or 3468 or 346
Aptidal_O3  = "/home/jeremy/Documents/Aptidal_simulation/test/To_be_kept/SABA1_346_O3.txt" #Change 235 by 112 or 3468 or 346
Aptidal_O3p = "/home/jeremy/Documents/Aptidal_simulation/test/To_be_kept/SABA1_346_O3+.txt"#Change 235 by 112 or 3468 or 346


tIPB, lbd1_IPB, k1_IPB, h1_IPB, lbd2_IPB, k2_IPB, h2_IPB, lbd3_IPB, k3_IPB, h3_IPB = np.loadtxt(intplabasic, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 3, 4, 8, 9, 10, 14, 15, 16]))
tATD_O1,  phi_1_ATD_O1,  e1_ATD_O1,  e2_ATD_O1,  e3_ATD_O1                         = np.loadtxt(Aptidal_O1,  dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 6, 12, 18]))
tATD_O2,  phi_1_ATD_O2,  e1_ATD_O2,  e2_ATD_O2,  e3_ATD_O2                         = np.loadtxt(Aptidal_O2,  dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 6, 12, 18]))
tATD_O3,  phi_1_ATD_O3,  e1_ATD_O3,  e2_ATD_O3,  e3_ATD_O3                         = np.loadtxt(Aptidal_O3,  dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 6, 12, 18]))
tATD_O3p, phi_1_ATD_O3p, e1_ATD_O3p, e2_ATD_O3p, e3_ATD_O3p                        = np.loadtxt(Aptidal_O3p, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 6, 12, 18]))

lbd1_IPB     = np.fmod(lbd1_IPB, 2.*m.pi)
lbd2_IPB     = np.fmod(lbd2_IPB, 2.*m.pi)
lbd3_IPB     = np.fmod(lbd3_IPB, 2.*m.pi)
#phi_1_IPB    = 2.*lbd1_IPB-9./2.*lbd2_IPB+5./2.*lbd3_IPB
phi_1_IPB    = lbd1_IPB-2.*lbd2_IPB+lbd3_IPB
#phi_1_IPB    = lbd1_IPB-lbd2_IPB
phi_1_ATD_O1 = 180./m.pi*phi_1_ATD_O1 + 360
phi_1_ATD_O2 = 180./m.pi*phi_1_ATD_O2 + 360
phi_1_ATD_O3 = 180./m.pi*phi_1_ATD_O3 + 360

#Fixing phi_1_ATD_O1 and phi_1_ATD_O2
phi_1_ATD_O1 = phi_1_ATD_O1 + 180
phi_1_ATD_O1 = np.fmod(phi_1_ATD_O1, 360)
phi_1_ATD_O1 = phi_1_ATD_O1 - 180
for i in range (0, len(phi_1_ATD_O2)):
      if(phi_1_ATD_O2[i] < -180):
            phi_1_ATD_O2[i] = phi_1_ATD_O2[i] + 360

e1_IPB       = np.sqrt(k1_IPB**2 + h1_IPB**2)
e2_IPB       = np.sqrt(k2_IPB**2 + h2_IPB**2)
e3_IPB       = np.sqrt(k3_IPB**2 + h3_IPB**2)

def Laplace_fix(phi):
      phi[0]   = np.fmod(phi[0] + 600.*np.pi, 2.*np.pi)
      max_jump = 0.8*m.pi
      jump     = 0
      for i in range (0, len(phi) - 1):
            jump = phi[i+1] - phi[i]
            while (jump > max_jump):
                  phi[i+1] = phi[i+1] - 2.*m.pi
                  jump = phi[i+1] - phi[i]
            while (jump < -max_jump):
                  phi[i+1] = phi[i+1] + 2.*m.pi
                  jump = phi[i+1] - phi[i]
      return phi           

phi_1_IPB = 180./m.pi*Laplace_fix(phi_1_IPB)
#phi_1_IPB = 180./m.pi*phi_1_IPB
   
#fig, (ax1, ax2, ax3) = py.subplots(1, 3, sharex=True, sharey=True, gridspec_kw={'width_ratios': [1,1,1]}, constrained_layout=False)
fig, ax1 = py.subplots()
py.subplots_adjust(left=0.2, right=0.7, bottom=0.14, top=0.9)

ax1.scatter(tIPB,    phi_1_IPB,     c='black', marker="o", s=25, alpha=1, label='Complete system')
ax1.scatter(tATD_O1, phi_1_ATD_O1,  c='blue',  marker="o", s=25, alpha=1, label='Averaged system at first order in eccentricity')
ax1.scatter(tATD_O2, phi_1_ATD_O2,  c='red',   marker="o", s=25, alpha=1, label='Averaged system at second order in eccentricity')
ax1.scatter(tATD_O3, phi_1_ATD_O3,  c='green', marker="o", s=25, alpha=1, label='Averaged system at third order in eccentricity')


ax1.set_xlim(xmin=0, xmax=999)
ax1.set_ylim(ymin=-180, ymax=180)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.tick_params(axis='y', which='major', labelsize=35)
ax1.set_xlabel(xlabel="Time (inner orbital period)", fontsize=25, labelpad=3)
#ax1.set_ylabel(ylabel=r"$\varphi_1=2\,\lambda_1-9/2\,\lambda_2+5/2\,\lambda_3$ ($^\circ$)", fontsize=25, labelpad=-4, rotation=90)
ax1.set_ylabel(ylabel=r"$\varphi_1=\lambda_1-2\,\lambda_2+\lambda_3$", fontsize=25, labelpad=-7, rotation=90)
#ax1.set_ylabel(ylabel=r"$\varphi_1=\lambda_1-\lambda_2$ ($^\circ$)", fontsize=25, labelpad=0, rotation=90)
#ax1.set_title(label=r"$$", loc='center', pad=6.0, y=0.96, color='black', fontsize=25, rotation=0)
ax1.grid(linewidth=0.3)

h1,labels = ax1.get_legend_handles_labels()
h1 = ['Complete system', 'Averaged system at first order in eccentricity', 'Averaged system at second order in eccentricity', 'Averaged system at third order in eccentricity']
ax1.legend(h1,fontsize = 20, loc='lower left')

'''# y-ticks
yticks=[-60, 0, 60, 120, 180, 240, 300]
ax1.set_yticks(yticks)
dic_y = {-60 : r"$-\frac{\pi}{3}$", 0 : r"$0$", 60 : r"$\frac{\pi}{3}$", 120 : r"$\frac{2\pi}{3}$", 180 : r"$\pi$", 240 : r"$\frac{4\pi}{3}$", 300 : r"$\frac{5\pi}{3}$"}
labels_y = [dic_y.get(T, yticks[i]) for i,T in enumerate(yticks)]
ax1.set_yticklabels(labels_y)'''

# y-ticks
yticks=[-180, -120, -60, 0, 60, 120, 180]
ax1.set_yticks(yticks)
dic_y = {-180 : r"$-\pi$", -120 : r"$-\frac{2\pi}{3}$", -60 : r"$-\frac{\pi}{3}$", 0 : r"$0$", 60 : r"$\frac{\pi}{3}$", 120 : r"$\frac{2\pi}{3}$", 180 : r"$\pi$"}
labels_y = [dic_y.get(T, yticks[i]) for i,T in enumerate(yticks)]
ax1.set_yticklabels(labels_y)

# x-ticks
'''xticks=[np.log10(0.5), 0, 1, 2, 3, 4]
ax1.set_xticks(xticks)
dic_x = {np.log10(0.5) : r"$\frac{1}{2}$", 0 : r"$1$", 1 : r"$10^{1}$", 2 : r"$10^{2}$", 3 : r"$10^{3}$", 4 : r"$10^{4}$"}
labels_x = [dic_x.get(T, xticks[i]) for i,T in enumerate(xticks)]
ax1.set_xticklabels(labels_x)'''

#py.legend(fontsize = 25)
#figure = py.gcf() ##get current figure
#figure.set_size_inches(16, 5)
#py.savefig("/home/jeremy/Documents/git/moon-formation/latex/image/error_distribution_theta_min=0.5.png", dpi=160)
py.show()








