import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('font',**{'size':16})

### HHG ###

#GET DATA
x=np.loadtxt('../HHG.dat',usecols=0)
y1=np.loadtxt('../HHG.dat',usecols=1)
y2=np.loadtxt('../HHG.dat',usecols=2)

"""
#LINEAR PLOTS
plt.plot(x,y1,label='Single Atom')
plt.plot(x,y2,label='Many Atoms')
plt.xlim(np.min(x),np.max(x))
plt.grid()
plt.legend()
#plt.title(r'Linear HHG Spectra')
plt.xlabel(r'$\omega\ (\omega_0)$')
plt.ylabel(r'$P$')
plt.tight_layout()
plt.savefig('hhg_lin.png')
plt.close()
"""

#LOG PLOTS
plt.semilogy(x,y1,label='Single Atom')
plt.semilogy(x,y2,label='Many Atoms')
plt.xlim(np.min(x),np.max(x))
plt.grid()
plt.legend()
#plt.title(r'Logarithmic HHG Spectra')
plt.xlabel(r'$\omega\ (\omega_0)$')
plt.ylabel(r'$P$')
plt.tight_layout()
plt.savefig('hhg_log.png')
plt.close()

### POPULATIONS ###

#GET DATA
x=np.loadtxt('../fort.100',usecols=0)
x=1.E15*x*2.4188843265857E-17
ydiss=np.ones(len(x))
y=[]
for i in range(19):
    yt=np.loadtxt('../fort.%03d'%(100+i),usecols=1)
    y.append(yt)
    ydiss=ydiss-yt
y=np.array(y)

plt.figure(figsize=(12.8,4.8))
for i in range(19):
    plt.plot(x,y[i],label=r'$v=%d$'%(i))
plt.plot(x,ydiss,label=r'Dissociative')
plt.xlim(np.min(x),np.max(x))
plt.grid()
plt.legend(ncol=2,bbox_to_anchor=(1.04,1),borderaxespad=0)
#plt.title(r'Bound States Population')
plt.xlabel(r'$t\ (fs)$')
plt.ylabel(r'$P$')
plt.tight_layout()
plt.savefig('pop.png')
plt.close()

### R_MOY ###

#GET DATA
x=np.loadtxt('../fort.7',usecols=0)
x=1.E15*x*2.4188843265857E-17
y=np.loadtxt('../fort.7',usecols=1)
yl=np.loadtxt('../fort.30',usecols=1)

plt.plot(x,y,label=r'$R_{moy}^{tot}$')
plt.plot(x,yl,label=r'$R_{moy}^{bound}$')
plt.xlim(np.min(x),np.max(x))
plt.grid()
plt.legend()
#plt.title(r'Average R')
plt.xlabel(r'$t\ (fs)$')
plt.ylabel(r'$R\ (a_0)$')
plt.tight_layout()
plt.savefig('rmoy.png')
plt.close()
