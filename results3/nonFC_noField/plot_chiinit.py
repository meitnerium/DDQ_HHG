import numpy as np
import matplotlib.pyplot as plt

chi1init = np.genfromtxt("chi1init.dat")
x=[]
y=[]
for row in chi1init:
    x.append(row[0])
    y.append(np.abs(row[1]+1j*row[2])**2)
    #print(x,y)


chi1initFC = np.genfromtxt("../../results2/FC_field-0/chi1init.dat")
xFC=[]
yFC=[]
for row in chi1initFC:
    xFC.append(row[0])
    yFC.append(np.abs(row[1]+1j*row[2])**2)
    #print(x,y)


plt.plot(xFC,yFC, label="FC")
plt.plot(x,y, label="non FC")
plt.legend()
plt.xlabel("R(u.a.)")
plt.ylabel(r"$|\Psi|^2$")
plt.xlim(1,3)
plt.savefig("chi1init.png")
plt.show()

