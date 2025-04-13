import math as m
import matplotlib.pyplot as py
import matplotlib as mpl
import numpy as np

FixedPoint = "/home/jeremy/Documents/Aptidal_simulation/GroupMeeting/UnaveragedSABA4_346_fixedPoint.txt"
Iteration1 = "/home/jeremy/Documents/Aptidal_simulation/GroupMeeting/UnaveragedSABA1_346.txt"
Iteration2 = "/home/jeremy/Documents/Aptidal_simulation/GroupMeeting/UnaveragedSABA2_346.txt"
Iteration3 = "/home/jeremy/Documents/Aptidal_simulation/GroupMeeting/UnaveragedSABA3_346.txt"


tPF, phiPF = np.loadtxt(FixedPoint, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2]))
tI1, phiI1 = np.loadtxt(Iteration1, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2]))
tI2, phiI2 = np.loadtxt(Iteration2, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2]))
tI3, phiI3 = np.loadtxt(Iteration3, dtype = np.float64, delimiter=' ', unpack=True, usecols=np.array([0, 2]))
phiPF = 180./m.pi*phiPF
phiI1 = 180./m.pi*phiI1
phiI2 = 180./m.pi*phiI2
phiI3 = 180./m.pi*phiI3

fig, ax1 = py.subplots()
py.subplots_adjust(left=0.2, right=0.7, bottom=0.14, top=0.9)

ax1.scatter(tPF, phiPF,  c='black', marker="o", s=125, alpha=0.6, label=r'Fixed point of the model at $\mathcal{O}(e_j^3)$')
ax1.scatter(tI1, phiI1,  c='blue',  marker="o", s=125, alpha=0.6, label=r'After $1$ iteration  of NAFF reduction')
ax1.scatter(tI2, phiI2,  c='red',   marker="o", s=125, alpha=0.6, label=r'After $2$ iterations of NAFF reduction')
#ax1.scatter(tI3, phiI3,  c='green', marker="o", s=25, alpha=0.6, label=r'After $3$ iterations of NAFF reduction')


#ax1.set_xlim(xmin=0, xmax=999)
#ax1.set_ylim(ymin=-180, ymax=180)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax1.tick_params(axis='y', which='major', labelsize=25)
ax1.set_xlabel(xlabel="Time (inner orbital period)", fontsize=25, labelpad=3)
ax1.set_ylabel(ylabel=r"$\varphi_1=\lambda_1-2\,\lambda_2+\lambda_3$ ($^\circ$)", fontsize=30, labelpad=-2, rotation=90)
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








