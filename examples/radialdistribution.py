import numpy as np
import sys
sys.path.insert(0, '../src/')
from fmt1d import FMT1D
import matplotlib.pyplot as plt
import timeit
import pandas as pd

starttime = timeit.default_timer()
plt.style.use(['science'])
    
fmt = FMT1D(L=9.0,method='RF',geometry='Spherical')
fmt2 = FMT1D(L=9.0,method='WBI',geometry='Spherical')
fmt3 = FMT1D(L=9.0,method='WBII',geometry='Spherical')
fmt.Set_External_Potential(extpotmodel='hardsphere')
fmt2.Set_External_Potential(extpotmodel='hardsphere')
fmt3.Set_External_Potential(extpotmodel='hardsphere')

rhobarray = np.array([0.2,0.9])

df = pd.read_excel('data/MCdata-radialdistribution-hardsphere-Barker1971.xls',sheet_name='radialdistributionfunction') 

for rhob in rhobarray:

    fmt.Set_BulkDensities(np.array([rhob]))
    fmt2.Set_BulkDensities(np.array([rhob]))
    fmt3.Set_BulkDensities(np.array([rhob]))

    print('rhob = ',rhob)
    print('muid = ',fmt.muid[0])
    print('muexc = ',fmt.muexc[0],fmt2.muexc[0],fmt3.muexc[0])

    fmt.Calculate_Equilibrium(logoutput=False)
    fmt2.Calculate_Equilibrium(logoutput=False)
    fmt3.Calculate_Equilibrium(logoutput=False)

    plt.scatter(df['r'],df['rhob='+str(rhob)],marker='o',edgecolors='C0',facecolors='none',label='MC')
    plt.plot(fmt.z,fmt.rho[0]/fmt.rhob[0],':',color='k',label='RF')
    plt.plot(fmt2.z,fmt2.rho[0]/fmt.rhob[0],'--',color='k',label='WBI')
    plt.plot(fmt3.z,fmt3.rho[0]/fmt.rhob[0],'-',color='k',label='WBIII')
    plt.xlabel(r'$r/d$')
    plt.ylabel(r'$g(r)$')
    plt.xlim(1,3)
    # plt.ylim(0.0,5.5)
    plt.text(2,df['rhob='+str(rhob)].max()/2,r'$\rho_b=$'+str(rhob))
    plt.legend(loc='best')
    plt.savefig('radialdistribution-rhob'+str(rhob)+'.png',dpi=200)
    # plt.show()
    plt.close()