import numpy as np
import sys
sys.path.insert(0, '../src/')
from fmt1d import FMT1D
import matplotlib.pyplot as plt
import timeit

starttime = timeit.default_timer()
plt.style.use(['science'])

d = np.array([1.0,3.0])
L = 10*d[1]
fmt = FMT1D(L=L,d=d,method='RF')
fmt.Set_External_Potential(extpotmodel='hardwall')
fmt1 = FMT1D(L=L,d=d,method='WBI')
fmt1.Set_External_Potential(extpotmodel='hardwall')
fmt2 = FMT1D(L=L,d=d,method='WBII')
fmt2.Set_External_Potential(extpotmodel='hardwall')

eta = 0.39
x1 = 0.25
r = d[0]/d[1]
rhob = np.array([eta/(np.pi*d[0]**3*(1+(1-x1)/x1/r**3)/6), eta/(np.pi*d[0]**3*(1+(1-x1)/x1/r**3)/6)*(1-x1)/x1])

fmt.Set_BulkDensities(rhob)
fmt1.Set_BulkDensities(rhob)
fmt2.Set_BulkDensities(rhob)

print('mu =',fmt.mu)
print('rhob=',fmt.rhob)

fmt.Calculate_Equilibrium(logoutput=False)
fmt1.Calculate_Equilibrium(logoutput=False)
fmt2.Calculate_Equilibrium(logoutput=False)

# data from: NOWORYTA, JERZY P., et al. “Hard Sphere Mixtures near a Hard Wall.” Molecular Physics, vol. 95, no. 3, Oct. 1998, pp. 415–24, doi:10.1080/00268979809483175.
folder = 'data/'
name = 'hardwall-mixture-eta=0.39-x1=0.25-ratio3'
data = np.loadtxt(folder+name+'.dat')
[xdata, rhodata,xdata1, rhodata1] = data.T
plt.scatter(xdata,rhodata*fmt.rhob[0],marker='o',edgecolors='C0',facecolors='none',label='MC')
plt.plot(fmt.z,fmt.rho[0],':k',label='RF')
plt.plot(fmt1.z,fmt1.rho[0],'--k',label='WBI')
plt.plot(fmt2.z,fmt2.rho[0],'k',label='WBII')
plt.xlabel(r'$z/d_s$')
plt.ylabel(r'$\rho_s(z) d_s^3$')
plt.xlim(0.5,8)
plt.ylim(0.005,0.025)
plt.text(4,0.017,r'$\eta=$'+str(eta))
plt.text(4,0.015,r'$x_s=$'+str(x1))
plt.text(4,0.013,r'$d_b/d_s=$'+str(d[1]/d[0]))
plt.legend(loc='upper right')
plt.savefig(name+'-small.png',dpi=200)
plt.show()
plt.close()

plt.scatter(xdata1,rhodata1*fmt.rhob[1]*d[1]**3,marker='s',edgecolors='C3',facecolors='none',label='MC')
plt.plot(fmt.z,fmt.rho[1]*d[1]**3,':k',label='RF')
plt.plot(fmt1.z,fmt1.rho[1]*d[1]**3,'--k',label='WBI')
plt.plot(fmt2.z,fmt2.rho[1]*d[1]**3,'k',label='WBII')
plt.xlabel(r'$z/d_s$')
plt.ylabel(r'$\rho_b(z) d_b^3$')
plt.xlim(1.5,8)
plt.ylim(0.0,5.2)
plt.text(4,3.0,r'$\eta=$'+str(eta))
plt.text(4,2.6,r'$x_s=$'+str(x1))
plt.text(4,2.2,r'$d_b/d_s=$'+str(d[1]/d[0]))
plt.legend(loc='upper right')
plt.savefig(name+'-big.png',dpi=200)
plt.show()
plt.close()
