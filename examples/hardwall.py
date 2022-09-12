import numpy as np
import sys
sys.path.insert(0, '../src/')
from fmt1d import FMT1D
import matplotlib.pyplot as plt

# plt.style.use(['science'])
    
fmt = FMT1D(L=6.0,method='WBI')
fmt.Set_External_Potential(extpotmodel='hardwall')

etaarray = np.array([0.4257,0.4783])

for eta in etaarray:

    rhob = eta/(np.pi/6.0)
    fmt.Set_BulkDensities(np.array([rhob]))

    fmt.Calculate_Equilibrium(logoutput=False)

    data = np.loadtxt('data/hardwall-eta'+str(eta)+'.dat')
    [xdata, rhodata] = [data[:,0],data[:,1]]
    plt.scatter(xdata,rhodata,marker='o',edgecolors='C0',facecolors='none',label='MC')
    plt.plot(fmt.z,fmt.rho[0],'-',color='k',label=fmt.method)
    plt.xlabel(r'$z/d$')
    plt.ylabel(r'$\rho(z)d^3$')
    plt.xlim(0.5,3)
    # plt.ylim(0.0,7)
    plt.text(2,rhodata.max()/2,r'$\eta=$'+str(eta))
    plt.legend(loc='best')
    plt.savefig('hardwall-eta'+str(eta)+'.png',dpi=200)
    plt.show()
    plt.close()