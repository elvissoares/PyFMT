import numpy as np
from scipy.ndimage import convolve1d
# Author: Elvis do A. Soares
# Github: @elvissoares
# Date: 2020-06-16
# Updated: 2022-09-12

twopi = 2*np.pi

def phi2func(eta):
    return np.piecewise(eta,[eta<=1e-3,eta>1e-3],[lambda eta: 1+eta**2/9,lambda eta: 1+(2*eta-eta**2+2*np.log(1-eta)*(1-eta))/(3*eta)])

def phi3func(eta):
    return np.piecewise(eta,[eta<=1e-3,eta>1e-3],[lambda eta: 1-4*eta/9,lambda eta: 1-(2*eta-3*eta**2+2*eta**3+2*np.log(1-eta)*(1-eta)**2)/(3*eta**2)])

def phi1func(eta):
    return np.piecewise(eta,[eta<=1e-3,eta>1e-3],[lambda eta: 1-2*eta/9-eta**2/18,lambda eta: 2*(eta+np.log(1-eta)*(1-eta)**2)/(3*eta**2)])

def dphi1dnfunc(eta):
    return np.piecewise(eta,[eta<=1e-3,eta>1e-3],[lambda eta: -2/9-eta/9-eta**2/15.0,lambda eta: (2*(eta-2)*eta+4*(eta-1)*np.log(1-eta))/(3*eta**3)])

def dphi2dnfunc(eta):
    return np.piecewise(eta,[eta<=1e-3,eta>1e-3],[lambda eta: 2*eta/9+eta**2/6.0,lambda eta: -(2*eta+eta**2+2*np.log(1-eta))/(3*eta**2)])

def dphi3dnfunc(eta):
    return np.piecewise(eta,[eta<=1e-3,eta>1e-3],[lambda eta: -4.0/9+eta/9,lambda eta: -2*(1-eta)*(eta*(2+eta)+2*np.log(1-eta))/(3*eta**3)])

def convolve1dplanar(rho,w,z):
    return convolve1d(rho,w, mode='nearest')

def convolve1dspherical(rho,w,r):
    return convolve1d(rho*r,w, mode='nearest')/r

def integrate1dplanar(f,z,dz):
    return np.sum(f)*dz

def integrate1dspherical(f,r,dr):
    return np.sum(f*4*np.pi*r**2)*dr

####################################################
####### The FMT Functional on 1d geometries   ######
####################################################
class FMT1D():
    def __init__(self,L,d=np.array([1.0]),method='WBI',geometry='Planar'):
        self.geometry = geometry
        self.method = method
        self.L = L
        self.d = d
        self.species = d.size

        self.delta = 0.01*self.d.min()
        self.z = np.arange(0.0,self.L,self.delta)+0.5*self.delta
        if self.geometry == 'Planar':
            self.convolve = convolve1dplanar
            self.integrate = integrate1dplanar
        elif self.geometry == 'Spherical':
            self.convolve = convolve1dspherical
            self.integrate = integrate1dspherical
        self.N = self.z.size

        self.rho = np.empty((self.species,self.N),dtype=np.float32)
        self.Vext = np.empty((self.species,self.N),dtype=np.float32)

        self.n3 = np.empty(self.N,dtype=np.float32)
        self.n2 = np.empty(self.N,dtype=np.float32)
        self.n2vec = np.empty(self.N,dtype=np.float32)

        self.w3 = np.empty(self.species,dtype=object)
        self.w2 = np.empty(self.species,dtype=object)
        self.w2vec = np.empty(self.species,dtype=object)
        self.c1hs = np.empty((self.species,self.N),dtype=np.float32)
        self.c1 = np.empty((self.species,self.N),dtype=np.float32)

        for i in range(self.species):
            nd = int(self.d[i]/self.delta)+1
            x = np.linspace(-0.5*self.d[i],0.5*self.d[i],nd)

            self.w3[i] = np.pi*((0.5*self.d[i])**2-x**2)
            self.w2[i] = self.d[i]*np.pi*np.ones(nd)
            self.w2vec[i] = twopi*x

    def Set_BulkDensities(self,rhob):
        self.rhob = rhob
        self.Calculate_mu()
        self.Set_InitialCondition()
    
    def Set_External_Potential(self,extpotmodel='hardwall',params='None'):
        self.extpotmodel = extpotmodel
        self.params = params
        if self.extpotmodel  == 'hardwall' or self.extpotmodel  == 'hardpore':
            self.Vext[:] = 0.0

    def Set_InitialCondition(self):
        nsig = (0.5*self.d/self.delta).astype(int)
        if self.extpotmodel  == 'hardwall':
            for i in range(self.species):
                self.rho[i,:] = self.rhob[i]
                self.rho[i,:nsig[i]] = 1.0e-16
        elif self.extpotmodel  == 'hardpore':
            for i in range(self.species):
                self.rho[i,:] = self.rhob[i]
                self.rho[i,:nsig[i]] = 1.0e-16
                self.rho[i,-nsig[i]:] = 1.0e-16
        self.Update_System()

    def Update_System(self):
        self.Calculate_weighted_densities()
        self.Calculate_c1()
        self.Calculate_Omega()

    def Calculate_weighted_densities(self):
        self.n3[:] = self.convolve(self.rho[0],self.w3[0],self.z)*self.delta
        self.n2[:] = self.convolve(self.rho[0],self.w2[0],self.z)*self.delta
        self.n2vec[:] = self.convolve(self.rho[0],self.w2vec[0],self.z)*self.delta
        if self.geometry == 'Spherical':
            self.n2vec[:] += self.n3/self.z
        self.n1vec = self.n2vec/(twopi*self.d[0])
        self.n0 = self.n2/(np.pi*self.d[0]**2)
        self.n1 = self.n2/(twopi*self.d[0])
        
        for i in range(1,self.species):
            n3 = self.convolve(self.rho[i],self.w3[i],self.z)*self.delta
            n2 = self.convolve(self.rho[i],self.w2[i],self.z)*self.delta
            n2vec = self.convolve(self.rho[i],self.w2vec[i],self.z)*self.delta
            self.n3[:] += n3
            self.n2[:] += n2
            self.n2vec[:] += n2vec
            if self.geometry == 'Spherical':
                self.n2vec[:] += n3/self.z
            self.n1vec[:] += n2vec/(twopi*self.d[i])
            self.n0[:] += n2/(np.pi*self.d[i]**2)
            self.n1[:] += n2/(twopi*self.d[i])
            
        self.oneminusn3 = 1-self.n3

        if self.method == 'RF' or self.method == 'WBI': 
            self.phi2 = 1.0
            self.dphi2dn3 = 0.0
        elif self.method == 'WBII': 
            self.phi2 = phi2func(self.n3)
            self.dphi2dn3 = dphi2dnfunc(self.n3)

        if self.method == 'WBI': 
            self.phi3 = phi1func(self.n3)
            self.dphi3dn3 = dphi1dnfunc(self.n3)
        elif self.method == 'WBII': 
            self.phi3 = phi3func(self.n3)
            self.dphi3dn3 = dphi3dnfunc(self.n3)
        else:
            self.phi3 = 1.0
            self.dphi3dn3 = 0.0

    def Calculate_Free_energy(self):
        self.Fid = self.integrate(self.rho*(np.log(self.rho)-1.0),self.z,self.delta)

        aux = -self.n0*np.log(self.oneminusn3)+(self.phi2/self.oneminusn3)*(self.n1*self.n2-(self.n1vec*self.n2vec)) + (self.phi3/(24*np.pi*self.oneminusn3**2))*(self.n2*self.n2*self.n2-3*self.n2*(self.n2vec*self.n2vec))

        self.Fhs = self.integrate(aux,self.z,self.delta)

        self.F = self.Fid+self.Fhs

    def Calculate_Omega(self):
        self.Calculate_Free_energy()
        self.Omega = (self.F + self.integrate((self.Vext-self.mu[:,np.newaxis])*self.rho,self.z,self.delta))/self.L

    def Calculate_c1(self):
        self.dPhidn0 = -np.log(self.oneminusn3 )
        self.dPhidn1 = self.n2*self.phi2/self.oneminusn3
        self.dPhidn2 = self.n1*self.phi2/self.oneminusn3  + (3*self.n2*self.n2-3*(self.n2vec*self.n2vec))*self.phi3/(24*np.pi*self.oneminusn3**2)

        self.dPhidn3 = self.n0/self.oneminusn3 +(self.n1*self.n2-(self.n1vec*self.n2vec))*(self.dphi2dn3 + self.phi2/self.oneminusn3)/self.oneminusn3 + (self.n2*self.n2*self.n2-3*self.n2*(self.n2vec*self.n2vec))*(self.dphi3dn3+2*self.phi3/self.oneminusn3)/(24*np.pi*self.oneminusn3**2)

        self.dPhidn1vec0 = -self.n2vec*self.phi2/self.oneminusn3 
        self.dPhidn2vec0 = -self.n1vec*self.phi2/self.oneminusn3  - self.n2*self.n2vec*self.phi3/(4*np.pi*self.oneminusn3**2)

        for i in range(self.species):
            self.c1hs[i,:] = -self.convolve(self.dPhidn2 + self.dPhidn1/(twopi*self.d[i]) + self.dPhidn0/(np.pi*self.d[i]**2),self.w2[i],self.z)*self.delta - self.convolve(self.dPhidn3,self.w3[i],self.z)*self.delta + self.convolve(self.dPhidn2vec0+self.dPhidn1vec0/(twopi*self.d[i]),self.w2vec[i],self.z)*self.delta
            if self.geometry == 'Spherical':
                self.c1hs[i,:] += -self.convolve(self.dPhidn2vec0+self.dPhidn1vec0/(twopi*self.d),self.w3,self.z)*self.delta/self.z

        del self.dPhidn0,self.dPhidn1,self.dPhidn2,self.dPhidn3,self.dPhidn1vec0,self.dPhidn2vec0,
        
        self.c1[:] = self.c1hs

    def Calculate_mu(self):
        self.muid = np.log(self.rhob)

        n3 = np.sum(self.rhob*np.pi*self.d**3/6)
        n2 = np.sum(self.rhob*np.pi*self.d**2)
        n1 = np.sum(self.rhob*self.d/2)
        n0 = np.sum(self.rhob)

        if self.method == 'RF' or self.method == 'WBI': 
            phi2 = 1.0
            dphi2dn3 = 0.0
        elif self.method == 'WBII': 
            phi2 = phi2func(n3)
            dphi2dn3 = dphi2dnfunc(n3)

        if self.method == 'WBI': 
            phi3 = phi1func(n3)
            dphi3dn3 = dphi1dnfunc(n3)
        elif self.method == 'WBII': 
            phi3 = phi3func(n3)
            dphi3dn3 = dphi3dnfunc(n3)
        else: 
            phi3 = 1.0
            dphi3dn3 = 0.0

        dPhidn0 = -np.log(1-n3)
        dPhidn1 = n2*phi2/(1-n3)
        dPhidn2 = n1*phi2/(1-n3) + (3*n2**2)*phi3/(24*np.pi*(1-n3)**2)
        dPhidn3 = n0/(1-n3) +(n1*n2)*(dphi2dn3 + phi2/(1-n3))/(1-n3) + (n2**3)*(dphi3dn3+2*phi3/(1-n3))/(24*np.pi*(1-n3)**2)

        self.muhs = dPhidn0+dPhidn1*self.d/2+dPhidn2*np.pi*self.d**2+dPhidn3*np.pi*self.d**3/6

        self.muexc = self.muhs

        self.mu = self.muid + self.muexc

    def Calculate_Equilibrium(self,alpha0=0.19,dt=0.01,atol=1e-8,rtol=1e-6,logoutput=False):
        " Global variables for the FIRE algorithm"
        Ndelay = 20
        Nmax = 10000
        finc = 1.1
        fdec = 0.5
        fa = 0.99
        Nnegmax = 2000
        dtmax = 10*dt
        dtmin = 0.02*dt
        alpha = alpha0
        Npos = 0
        Nneg = 0

        lnrho = np.log(self.rho)
        V = np.zeros_like(self.rho)
        F = -self.rho*(lnrho -self.c1 - self.mu[:,np.newaxis]+self.Vext)

        error0 = max(np.abs(F.min()),F.max())

        for i in range(Nmax):

            P = (F*V).sum() # dissipated power
            if (P>0):
                Npos = Npos + 1
                if Npos>Ndelay:
                    dt = min(dt*finc,dtmax)
                    alpha = alpha*fa
            else:
                Npos = 0
                Nneg = Nneg + 1
                if Nneg > Nnegmax: break
                if i> Ndelay:
                    dt = max(dt*fdec,dtmin)
                    alpha = alpha0
                lnrho[:] += - 0.5*dt*V
                V[:] = 0.0
                self.rho[:] = np.exp(lnrho)
                self.Update_System()

            V[:] += 0.5*dt*F
            V[:] = (1-alpha)*V + alpha*F*np.linalg.norm(V)/np.linalg.norm(F)
            lnrho[:] += dt*V
            self.rho[:] = np.exp(lnrho)
            self.Update_System()
            F[:] = -self.rho*(lnrho -self.c1 - self.mu[:,np.newaxis]+self.Vext)
            V[:] += 0.5*dt*F

            error = max(np.abs(F.min()),F.max())
            if error/error0 < rtol or error < atol: break

            if logoutput: print(i,self.Omega,error)
        self.Niter = i

        del V, F  


####################################################
class FMTspherical():
    def __init__(self,N,delta,d=1.0,method='WBI'):
        self.method = method
        self.N = N
        self.delta = delta
        self.L = N*delta
        self.d = d

        self.n3 = np.empty(self.N,dtype=np.float32)
        self.n2 = np.empty(self.N,dtype=np.float32)
        self.n2vec = np.empty(self.N,dtype=np.float32)

        self.r = np.arange(0,self.L,self.delta)+ 0.5*self.delta

        self.nsig = int(self.d/self.delta)

        self.w3 = np.zeros(self.nsig,dtype=np.float32)
        self.w2 = np.zeros(self.nsig,dtype=np.float32)
        self.w2vec = np.zeros(self.nsig,dtype=np.float32)
        
        r = np.linspace(-0.5*self.d,0.5*self.d,self.nsig)
        
        self.w3[:] = np.pi*((0.5*self.d)**2-r**2)
        self.w2[:] = self.d*np.pi
        self.w2vec[:] = twopi*r

    def weighted_densities(self,rho):
        self.n3[:] = convolve1d(rho*self.r, weights=self.w3, mode='nearest')*self.delta/self.r
        self.n2[:] = convolve1d(rho*self.r, weights=self.w2, mode='nearest')*self.delta/self.r
        self.n2vec[:] = self.n3/self.r + convolve1d(rho*self.r, weights=self.w2vec, mode='nearest')*self.delta/self.r

        self.n1vec = self.n2vec/(twopi*self.d)
        self.n0 = self.n2/(np.pi*self.d**2)
        self.n1 = self.n2/(twopi*self.d)
        self.oneminusn3 = 1-self.n3

        if self.method == 'RF' or self.method == 'WBI': 
            self.phi2 = 1.0
            self.dphi2dn3 = 0.0
        elif self.method == 'WBII': 
            self.phi2 = phi2func(self.n3)
            self.dphi2dn3 = dphi2dnfunc(self.n3)

        if self.method == 'WBI': 
            self.phi3 = phi1func(self.n3)
            self.dphi3dn3 = dphi1dnfunc(self.n3)
        elif self.method == 'WBII': 
            self.phi3 = phi3func(self.n3)
            self.dphi3dn3 = dphi3dnfunc(self.n3)
        else:
            self.phi3 = 1.0
            self.dphi3dn3 = 0.0

    def Phi(self,rho):
        self.weighted_densities(rho)
        return -self.n0*np.log(self.oneminusn3)+(self.phi2/self.oneminusn3)*(self.n1*self.n2-(self.n1vec*self.n2vec)) + (self.phi3/(24*np.pi*self.oneminusn3**2))*(self.n2*self.n2*self.n2-3*self.n2*(self.n2vec*self.n2vec))

    def F(self,rho):
        return np.sum(self.Phi(rho)*4*np.pi*self.r**2)*self.delta

    def dPhidn(self,rho):
        self.weighted_densities(rho)

        self.dPhidn0 = -np.log(self.oneminusn3 )
        self.dPhidn1 = self.n2*self.phi2/self.oneminusn3
        self.dPhidn2 = self.n1*self.phi2/self.oneminusn3  + (3*self.n2*self.n2-3*(self.n2vec*self.n2vec))*self.phi3/(24*np.pi*self.oneminusn3**2)

        self.dPhidn3 = self.n0/self.oneminusn3 +(self.n1*self.n2-(self.n1vec*self.n2vec))*(self.dphi2dn3 + self.phi2/self.oneminusn3)/self.oneminusn3 + (self.n2*self.n2*self.n2-3*self.n2*(self.n2vec*self.n2vec))*(self.dphi3dn3+2*self.phi3/self.oneminusn3)/(24*np.pi*self.oneminusn3**2)

        self.dPhidn1vec0 = -self.n2vec*self.phi2/self.oneminusn3 
        self.dPhidn2vec0 = -self.n1vec*self.phi2/self.oneminusn3  - self.n2*self.n2vec*self.phi3/(4*np.pi*self.oneminusn3**2)

        dPhidn = convolve1d((self.dPhidn2 + self.dPhidn1/(twopi*self.d) + self.dPhidn0/(np.pi*self.d**2))*self.r, weights=self.w2, mode='nearest')*self.delta/self.r
        dPhidn += convolve1d(self.dPhidn3*self.r, weights=self.w3, mode='nearest')*self.delta/self.r
        dPhidn += convolve1d((self.dPhidn2vec0+self.dPhidn1vec0/(twopi*self.d))*self.r, weights=self.w3, mode='nearest')*self.delta/(self.r)**2 - convolve1d((self.dPhidn2vec0+self.dPhidn1vec0/(twopi*self.d))*self.r, weights=self.w2vec, mode='nearest')*self.delta/self.r

        del self.dPhidn0,self.dPhidn1,self.dPhidn2,self.dPhidn3,self.dPhidn1vec0,self.dPhidn2vec0
        
        return dPhidn

    def c1(self,rho):
        return -self.dPhidn(rho)

    def mu(self,rhob):
        n3 = rhob*np.pi*self.d**3/6
        n2 = rhob*np.pi*self.d**2
        n1 = rhob*self.d/2
        n0 = rhob

        if self.method == 'RF' or self.method == 'WBI': 
            phi2 = 1.0
            dphi2dn3 = 0.0
        elif self.method == 'WBII': 
            phi2 = phi2func(n3)
            dphi2dn3 = dphi2dnfunc(n3)

        if self.method == 'WBI': 
            phi3 = phi1func(n3)
            dphi3dn3 = dphi1dnfunc(n3)
        elif self.method == 'WBII': 
            phi3 = phi3func(n3)
            dphi3dn3 = dphi3dnfunc(n3)
        else: 
            phi3 = 1.0
            dphi3dn3 = 0.0

        dPhidn0 = -np.log(1-n3)
        dPhidn1 = n2*phi2/(1-n3)
        dPhidn2 = n1*phi2/(1-n3) + (3*n2**2)*phi3/(24*np.pi*(1-n3)**2)
        dPhidn3 = n0/(1-n3) +(n1*n2)*(dphi2dn3 + phi2/(1-n3))/(1-n3) + (n2**3)*(dphi3dn3+2*phi3/(1-n3))/(24*np.pi*(1-n3)**2)

        return (dPhidn0+dPhidn1*self.d/2+dPhidn2*np.pi*self.d**2+dPhidn3*np.pi*self.d**3/6)

##### Take a example using FMT ######
if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import pandas as pd
    import timeit

    starttime = timeit.default_timer()
    plt.style.use(['science'])
    
    #############################
    test0 = False # optimizing alpha and delta
    test1 = False # hard wall (1D-planar)
    test2 = True # hard wall mixture (1D-planar)
    test3 = False # hard-sphere RDF (1D-spherical)
    test4 = False

    ######################################################
    if test0:
        fmt = FMT1D(L=6.0,method='WBI')
        fmt.Set_External_Potential(extpotmodel='hardwall')
        eta = 0.4783
        rhob = eta/(np.pi/6.0)
        fmt.Set_BulkDensities(np.array([rhob]))

        print('rhob = ',rhob)
        print('muid = ',fmt.muid[0])
        print('muexc = ',fmt.muexc[0])

        # alphaarray = np.linspace(0.1,0.5,41)
        # Niterarray = np.empty_like(alphaarray)

        # for i in range(alphaarray.size):
        #     fmt.Set_InitialCondition()
        #     fmt.Calculate_Equilibrium(alpha0=alphaarray[i],logoutput=False)
        #     Niterarray[i] = fmt.Niter

        # print('alpha=',alphaarray[Niterarray==Niterarray.min()])
        # plt.plot(alphaarray,Niterarray)
        # plt.show()

        deltaarray = np.power(10.0,np.linspace(-4,-2,10))
        Niterarray = np.empty_like(deltaarray)

        for i in range(deltaarray.size):
            fmt.Set_InitialCondition()
            fmt.Calculate_Equilibrium(dt=deltaarray[i],logoutput=False)
            Niterarray[i] = fmt.Niter

        print('delta=',deltaarray[Niterarray==Niterarray.min()])
        plt.plot(deltaarray,Niterarray)
        plt.show()

    ######################################################
    if test1:
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

    if test2:
        d = np.array([1.0,3.0])
        L = 10*d[1]
        fmt = FMT1D(L=L,d=d,method='RF')
        fmt.Set_External_Potential(extpotmodel='hardwall')
        fmt1 = FMT1D(L=L,d=d,method='WBI')
        fmt1.Set_External_Potential(extpotmodel='hardwall')
        fmt2 = FMT1D(L=L,d=d,method='WBII')
        fmt2.Set_External_Potential(extpotmodel='hardwall')

        eta = 0.12
        x1 = 0.83
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
        name = 'hardwall-mixture-eta=0.12-x1=0.83-ratio3'
        # name = 'hardwall-mixture-eta=0.39-x1=0.25-ratio3'
        data = np.loadtxt(folder+name+'.dat')
        [xdata, rhodata,xdata1, rhodata1] = data.T
        plt.scatter(xdata,rhodata*fmt.rhob[0],marker='o',edgecolors='C0',facecolors='none',label='MC')
        plt.plot(fmt.z,fmt.rho[0],':k',label='RF')
        plt.plot(fmt1.z,fmt1.rho[0],'--k',label='WBI')
        plt.plot(fmt2.z,fmt2.rho[0],'k',label='WBII')
        plt.xlabel(r'$z/d_s$')
        plt.ylabel(r'$\rho_s(z) d_s^3$')
        plt.xlim(0.5,8)
        plt.ylim(0.03,0.05)
        plt.text(4,0.047,r'$\eta=$'+str(eta))
        plt.text(4,0.045,r'$x_s=$'+str(x1))
        plt.text(4,0.043,r'$d_b/d_s=$'+str(d[1]/d[0]))
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
        plt.ylim(0.15,0.4)
        plt.text(4,0.3,r'$\eta=$'+str(eta))
        plt.text(4,0.27,r'$x_s=$'+str(x1))
        plt.text(4,0.24,r'$d_b/d_s=$'+str(d[1]/d[0]))
        plt.legend(loc='upper right')
        plt.savefig(name+'-big.png',dpi=200)
        plt.show()
        plt.close()

    ################################################################
    if test3:
        sigma = 1.0 
        delta = 0.01*sigma
        N = 900
        L = N*delta
        fmt = FMTspherical(N,delta,d=sigma)
        rhobarray = np.arange(0.1,0.95,0.01)

        nsig = int(sigma/delta)

        n0 = np.ones(N,dtype=np.float32)
        n0[:nsig] = 1.0e-16
        # n[N-nsig:] = 1.0e-12

        r = np.arange(0,L,delta)
        Vol = 4*np.pi*L**3/3
            
        def Omega(lnn,mu):
            n[nsig:-2*nsig] = np.exp(lnn)
            phi = fmt.Phi(n)
            Omegak = n*(np.log(n)-1.0) + phi - mu*n
            return np.sum(4*np.pi*r**2*Omegak*delta)/Vol

        def dOmegadnR(lnn,mu):
            n[nsig:-2*nsig] = np.exp(lnn)
            dOmegadrho = n*4*np.pi*r**2*(np.log(n) -fmt.c1(n)- mu)*delta
            return dOmegadrho[nsig:-2*nsig]/Vol

        n = n0.copy()

        for rhob in rhobarray:

            mu = np.log(rhob) + fmt.mu(rhob)

            lnn = np.log(n0[nsig:-2*nsig])+np.log(rhob)
        
            [lnnsol,Omegasol,Niter] = optimize_fire2(lnn,Omega,dOmegadnR,mu,alpha0=0.62,rtol=1.0e-6,dt=2.0,logoutput=True)

            n[nsig:-2*nsig] = np.exp(lnnsol)

            nmean = np.sum(n*4*np.pi*r**2*delta)/Vol
            print('rhob=',rhob,'\n nmean = ',nmean,'\n Omega/N =',Omegasol)

            np.save('results/radialdistribution-fmt-wbi-rhob'+str(rhob)+'-N'+str(N)+'-delta'+str(delta)+'.npy',[r,n])

            # data = np.loadtxt('../MCdataHS/hardwall-eta'+str(eta)+'.dat')
            # [xdata, rhodata] = [data[:,0],data[:,1]]
            # plt.scatter(xdata,rhodata,marker='o',edgecolors='C0',facecolors='none',label='MC')
            # plt.plot(r,n/rhob,label='DFT')
            # plt.xlabel(r'$r/\d$')
            # plt.ylabel(r'$g(r)$')
            # plt.xlim(1.0,2.2)
            # plt.ylim(0.5,6)
            # plt.legend(loc='best')
            # plt.savefig('hardsphere-eta'+str(eta)+'.png', bbox_inches='tight')
            # plt.show()
            # plt.close()

    ################################################################
    if test4:

        starttime = timeit.default_timer()

        rhob = 0.84
        kT = 0.71

        sigma = 1.0
        d = sigma*(1+0.2977*kT)/(1+0.33163*kT+0.0010477*kT**2)
        delta = 0.01*d
        L = 9.0*d
        N = int(L/delta)
        fmt = FMTspherical(N,delta,d=sigma)

        def Vext(r,eps,sigma):
            return 4*eps*((sigma/r)**12-(sigma/r)**6)

        nsig = int(0.5*sigma/delta)
        r = np.arange(0,L,delta)+0.5*delta

        n0 = np.ones(N,dtype=np.float32)
        n0[:nsig] = 1.0e-16
        # n[N-nsig:] = 1.0e-12
        eps = 0.71
        beta = 1.0/kT
        V = beta*Vext(r,eps,sigma)
        V[V>beta*1e6] = beta*1e6
        
        Vol = 4*np.pi*L**3/3
            
        def Omega(lnn,mu):
            n = np.exp(lnn)
            phi = fmt.Phi(n)
            Omegak = n*(np.log(n)-1.0) + phi - mu*n + V*n
            return np.sum(4*np.pi*r**2*Omegak*delta)/Vol

        def dOmegadnR(lnn,mu):
            n = np.exp(lnn)
            dOmegadrho = n*4*np.pi*r**2*(np.log(n) -fmt.c1(n)- mu + V)*delta
            return dOmegadrho/Vol

        n = n0.copy()

        mu = np.log(rhob) + fmt.mu(rhob)

        lnn = np.log(n0)+np.log(rhob)
    
        [lnnsol,Omegasol,Niter] = optimize_fire2(lnn,Omega,dOmegadnR,mu,alpha0=0.62,rtol=1.0e-5,dt=2.0,logoutput=True)

        n = np.exp(lnnsol)

        nmean = np.sum(n*4*np.pi*r**2*delta)/Vol
        print('rhob=',rhob)
        print('rhomean = ',nmean)

        print("time elapsed:", timeit.default_timer() - starttime, 'sec')

        df = pd.read_excel('MCdata/MCdata-radialdistribution-lennardjones-Verlet1968.xls',sheet_name='Argon') # Shukla2000

        plt.scatter(df['r']/3.405,df['KT=0.71-rhob=0.84'],marker='o',edgecolors='C0',facecolors='none',label=r'${}^{36}$Ar @ 85 K')
        plt.plot(r,np.exp(-V),'--',color='C3',label='Ideal')
        plt.plot(r,n/rhob,'k',label='DFT')
        plt.xlabel(r'$r/\sigma$')
        plt.ylabel(r'$g(r)$')
        plt.text(3.55,2.0,r'$k_B T/\epsilon = 0.75$')
        plt.text(3.6,1.7,r'$\rho_b \sigma^3 = 0.84$')
        plt.xlim(0.0,8.0)
        plt.ylim(0,3.5)
        plt.legend(loc='best')
        # plt.savefig('hardsphere-eta'+str(eta)+'.png', bbox_inches='tight')
        plt.show()
        # plt.close()
    
    print("time :", timeit.default_timer() - starttime, 'sec')