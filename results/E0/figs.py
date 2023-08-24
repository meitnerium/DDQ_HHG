import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('font',**{'size':16})

### R_MOY ###

#GET DATA
xfc=np.loadtxt('FC/fort.7',usecols=0)
xfcfs=1.E15*xfc*2.4188843265857E-17
yfc=np.loadtxt('FC/fort.7',usecols=1)
ylfc=np.loadtxt('FC/fort.30',usecols=1)

#plt.plot(x,y,label=r'$R_{moy}^{tot}$')
plt.plot(xfcfs,ylfc,label=r'$R_{moy}^{bound} (FC)$')
#plt.xlim(np.min(xfc),np.max(xfc))
plt.xlim(0,3000)

plt.grid()
plt.legend()
#plt.title(r'Average R')
plt.xlabel(r'$t\ (fs)$')
plt.ylabel(r'$R\ (a_0)$')
plt.tight_layout()
#plt.show()
#plt.savefig('fc.png')
#plt.close()


xnfc=np.loadtxt('nonFC/fort.7',usecols=0)
xnfcfs=1.E15*xnfc*2.4188843265857E-17
#ynfcfs=np.loadtxt('nonFC/fort.7',usecols=1)
ylnfc=np.loadtxt('nonFC/fort.30',usecols=1)
plt.plot(xnfcfs,ylnfc,label=r'$R_{moy}^{bound} (no FC)$')
plt.xlim(0,250)

plt.grid()
plt.legend()
#plt.title(r'Average R')
plt.xlabel(r'$t\ (fs)$')
plt.ylabel(r'$R\ (a_0)$')
plt.tight_layout()
#plt.show()
plt.savefig('fcvsnon-fc.png')
plt.close()

#plt.plot(x,y,label=r'$R_{moy}^{tot}$')




#plt.grid()
#plt.legend()
#plt.title(r'Average R')
#plt.xlabel(r'$t\ (fs)$')
#plt.ylabel(r'$R\ (a_0)$')
#plt.tight_layout()
#plt.show()
#plt.savefig('rmoy.png')
#plt.close()



#w800=0.057       
#t4field = np.arange(0,xfc[-1],10000)
#../nonFC-0_phase-0/fort.2
field=np.loadtxt('../nonFC-0_phase-0/fort.2',usecols=1)
#field800np_phase0 = np.cos(w800*t4field)*np.max(np.max(ylfc))#,np.max(ylnfc))
#plt.plot(ylfc,field800np_phase0,label="Field 800 nm")
a=1
b=10
plt.plot(xfc,b*field+a,label="Field")

#plt.xlim(np.min(xfc),np.max(xfc))
plt.xlim(0,3000)

plt.grid()
plt.legend()
#plt.title(r'Average R')
plt.xlabel(r'$t\ (fs)$')
plt.ylabel(r'$R\ (a_0)$')
plt.tight_layout()
#plt.show()
plt.savefig('field.png')
plt.close()
