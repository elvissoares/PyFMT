import numpy as np
import sys
sys.path.insert(0, '../src/')
from fmt1d import FMT1D
import matplotlib.pyplot as plt
import pandas as pd

# plt.style.use(['science'])
    
fmt = FMT1D(L=10.0,method='WBI',geometry='Spherical')
fmt.Set_External_Potential(extpotmodel='hardsphere')

rhobarray = np.array([0.2,0.9])

df = pd.read_excel('data/MCdata-radialdistribution-hardsphere-Barker1971.xls',sheet_name='radialdistributionfunction') 

for rhob in rhobarray:

    fmt.Set_BulkDensities(np.array([rhob]))

    fmt.Calculate_Equilibrium(logoutput=False)

    plt.scatter(df['r'],df['rhob='+str(rhob)],marker='o',edgecolors='C0',facecolors='none',label='MC')
    plt.plot(fmt.z,fmt.rho[0]/fmt.rhob[0],color='k',label='WBI')
    plt.xlabel(r'$r/d$')
    plt.ylabel(r'$g(r)$')
    plt.xlim(1,3)
    # plt.ylim(0.0,5.5)
    plt.text(2,df['rhob='+str(rhob)].max()/2,r'$\rho_b d^3=$'+str(rhob))
    plt.legend(loc='best')
    plt.savefig('radialdistribution-rhob'+str(rhob)+'.png',dpi=200)
    # plt.show()
    plt.close()