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
        n2sig = (self.d/self.delta).astype(int)
        if self.extpotmodel  == 'hardwall':
            for i in range(self.species):
                self.rho[i,:] = self.rhob[i]
                self.rho[i,:nsig[i]] = 1.0e-16
        elif self.extpotmodel  == 'hardpore':
            for i in range(self.species):
                self.rho[i,:] = self.rhob[i]
                self.rho[i,:nsig[i]] = 1.0e-16
                self.rho[i,-nsig[i]:] = 1.0e-16
        elif self.extpotmodel  == 'hardsphere':
            for i in range(self.species):
                self.rho[i,:] = self.rhob[i]
                self.rho[i,:n2sig[i]] = 1.0e-16
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
                self.c1hs[i,:] += -self.convolve(self.dPhidn2vec0+self.dPhidn1vec0/(twopi*self.d[i]),self.w3[i],self.z)*self.delta/self.z

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