import numpy as np
import sys
sys.path.insert(0, '../src/')
from fmt1d import FMT1D
import matplotlib.pyplot as plt
import timeit

starttime = timeit.default_timer()
plt.style.use(['science'])
    
fmtmethod = ['RF','WBI','WBII']
fmt = FMT1D(L=6.0,method='RF')
fmt2 = FMT1D(L=6.0,method='WBI')
fmt3 = FMT1D(L=6.0,method='WBII')
fmt.Set_External_Potential(extpotmodel='hardwall')
fmt2.Set_External_Potential(extpotmodel='hardwall')
fmt3.Set_External_Potential(extpotmodel='hardwall')

etaarray = np.array([0.4257,0.4783])

for eta in etaarray:

    rhob = eta/(np.pi/6.0)
    fmt.Set_BulkDensities(np.array([rhob]))
    fmt2.Set_BulkDensities(np.array([rhob]))
    fmt3.Set_BulkDensities(np.array([rhob]))

    print('rhob = ',rhob)
    print('muid = ',fmt.muid[0])
    print('muexc = ',fmt.muexc[0],fmt2.muexc[0],fmt3.muexc[0])

    fmt.Calculate_Equilibrium(logoutput=False)
    fmt2.Calculate_Equilibrium(logoutput=False)
    fmt3.Calculate_Equilibrium(logoutput=False)

    data = np.loadtxt('data/hardwall-eta'+str(eta)+'.dat')
    [xdata, rhodata] = [data[:,0],data[:,1]]
    plt.scatter(xdata,rhodata,marker='o',edgecolors='C0',facecolors='none',label='MC')
    plt.plot(fmt.z,fmt.rho[0],':',color='k',label='RF')
    plt.plot(fmt2.z,fmt2.rho[0],'--',color='k',label='WBI')
    plt.plot(fmt3.z,fmt3.rho[0],'-',color='k',label='WBIII')
    plt.xlabel(r'$z/d$')
    plt.ylabel(r'$\rho(z)d^3$')
    plt.xlim(0.5,3)
    # plt.ylim(0.0,7)
    plt.text(2,rhodata.max()/2,r'$\eta=$'+str(eta))
    plt.legend(loc='best')
    plt.savefig('hardwall-eta'+str(eta)+'.png',dpi=200)
    # plt.show()
    plt.close()