import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('font',**{'size':16})

### R_MOY ###

#GET DATA
xfc=np.loadtxt('FC/fort.7',usecols=0)
xfc=1.E15*xfc*2.4188843265857E-17
yfc=np.loadtxt('FC/fort.7',usecols=1)
ylfc=np.loadtxt('FC/fort.30',usecols=1)

#plt.plot(x,y,label=r'$R_{moy}^{tot}$')
plt.plot(xfc,ylfc,label=r'$R_{moy}^{bound} (FC)$')
#plt.xlim(np.min(xfc),np.max(xfc))


xnfc=np.loadtxt('nonFC/fort.7',usecols=0)
xnfc=1.E15*xnfc*2.4188843265857E-17
ynfc=np.loadtxt('nonFC/fort.7',usecols=1)
ylnfc=np.loadtxt('nonFC/fort.30',usecols=1)

#plt.plot(x,y,label=r'$R_{moy}^{tot}$')
plt.plot(xnfc,ylnfc,label=r'$R_{moy}^{bound} (no FC)$')
plt.xlim(np.min(xnfc),np.max(xnfc))






plt.grid()
plt.legend()
#plt.title(r'Average R')
plt.xlabel(r'$t\ (fs)$')
plt.ylabel(r'$R\ (a_0)$')
plt.tight_layout()
#plt.show()
plt.savefig('rmoy.png')
plt.close()
