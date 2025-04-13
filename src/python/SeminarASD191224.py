import math as m
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np

FixedPointAveraged   = "/home/jeremy/Documents/Aptidal_simulation/SABA4_3468_fixedPoint.txt"
FixedPointUnaveraged = "/home/jeremy/Documents/Aptidal_simulation/UnaveragedSABA4_3468_fixedPoint.txt"
LibrationCenter      = "/home/jeremy/Documents/Aptidal_simulation/UnaveragedSABA4_3468_LibrationCenter.txt"



tFPA, phi1FPA, phi2FPA = np.loadtxt(FixedPointAveraged,   dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 9]))
tFPU, phi1FPU, phi2FPU = np.loadtxt(FixedPointUnaveraged, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 9]))
tLC,  phi1LC,  phi2LC  = np.loadtxt(LibrationCenter,      dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2, 9]))
phi1FPA = 180./m.pi*(phi1FPA + 4*m.pi)
phi1FPU = 180./m.pi*(phi1FPU + 4*m.pi)
phi2FPA = 180./m.pi*(phi2FPA + 2*m.pi)
phi2FPU = 180./m.pi*(phi2FPU + 2*m.pi)
phi1LC  = 180./m.pi*(phi1LC  + 4*m.pi)
phi2LC  = 180./m.pi*(phi2LC  + 4*m.pi)

fig, ax1 = py.subplots()
py.subplots_adjust(left=0.1, right=0.9, bottom=0.22, top=0.8)

'''ax1.scatter(tFPA, phi1FPA, c='red',  marker="o", s=15, alpha=0.01)
ax1.scatter(tFPA, phi2FPA, c='blue', marker="o", s=15, alpha=0.01)
ax1.scatter(tFPU, phi1FPU, c='red',  marker="o", s=25, alpha=0.6, label=r'$\phi_1$')
ax1.scatter(tFPU, phi2FPU, c='blue', marker="o", s=25, alpha=0.6, label=r'$\phi_2$')'''
#ax1.scatter(tLC,  phi1LC,  c='red',  marker="o", s=75, alpha=0.6, label=r'$\phi_1$')
ax1.scatter(tLC,  phi2LC,  c='blue', marker="o", s=75, alpha=0.6, label=r'$\phi_2$')
#ax1.scatter(tI3, phiI3,  c='green', marker="o", s=25, alpha=0.6, label=r'After $3$ iterations of NAFF reduction')


#ax1.set_xlim(xmin=0, xmax=999)
#ax1.set_ylim(ymin=-180, ymax=180)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.tick_params(axis='y', which='major', labelsize=25)
ax1.set_xlabel(xlabel="Time (inner orbital period)", fontsize=25, labelpad=3)
ax1.set_ylabel(ylabel=r"$\phi_1$ ($^\circ$)", fontsize=30, labelpad=-2, rotation=90)
#ax1.grid(linewidth=0.3, alpha = 0.5)

'''h1,labels = ax1.get_legend_handles_labels()
h1 = [r'Fixed point of the model at $\mathcal{O}(e_j^3)$', r'After $1$ iteration  of NAFF reduction', r'After $2$ iterations of NAFF reduction', r'After $3$ iterations of NAFF reduction']
ax1.legend(h1,fontsize = 20, loc='lower left')'''

'''# y-ticks
yticks=[-60, 0, 60, 120, 180, 240, 300]
ax1.set_yticks(yticks)
dic_y = {-60 : r"$-\frac{\pi}{3}$", 0 : r"$0$", 60 : r"$\frac{\pi}{3}$", 120 : r"$\frac{2\pi}{3}$", 180 : r"$\pi$", 240 : r"$\frac{4\pi}{3}$", 300 : r"$\frac{5\pi}{3}$"}
labels_y = [dic_y.get(T, yticks[i]) for i,T in enumerate(yticks)]
ax1.set_yticklabels(labels_y)'''

'''# y-ticks
yticks=[-180, -120, -60, 0, 60, 120, 180]
ax1.set_yticks(yticks)
dic_y = {-180 : r"$-\pi$", -120 : r"$-\frac{2\pi}{3}$", -60 : r"$-\frac{\pi}{3}$", 0 : r"$0$", 60 : r"$\frac{\pi}{3}$", 120 : r"$\frac{2\pi}{3}$", 180 : r"$\pi$"}
labels_y = [dic_y.get(T, yticks[i]) for i,T in enumerate(yticks)]
ax1.set_yticklabels(labels_y)'''

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








