# PyFMT
A python implementation of the Fundamental Measure Theory for hard-sphere mixture in classical Density Functional Theory

For a fluid composed by hard-spheres with temperature T, total volume V, and chemical potential of each species $\mu_i$  specified, the grand potential, $\Omega$, is written as

$$\Omega[\{\rho_i(\boldsymbol{r})\}] = F^\text{id}[\{\rho_i(\boldsymbol{r})\}] + F^\text{hs}[\{\rho_i(\boldsymbol{r})\}]+ \sum_i \int_{V} d \boldsymbol{r} [V_i^{(\text{ext})}(\boldsymbol{r})-\mu_i] \rho_i(\boldsymbol{r})$$

where $\rho_i(\boldsymbol{r})$ is the local density of the component i, and $V^\text{ext}_{i}$ is the external potential. 

The ideal-gas contribution $F^\text{id}$ is given by the exact expression

$$F^\text{id}[\{\rho_i (\boldsymbol{r})\}] = k_B T\sum_i \int_{V} d\boldsymbol{r}\ \rho_i(\boldsymbol{r})[\ln(\rho_i (\boldsymbol{r})\Lambda_i^3)-1]$$

where $k_B$, and $\Lambda_i$ is the well-known thermal de Broglie wavelength of each component.

The hard-sphere contribution, $F^{\textrm{hs}}$, represents the hard-sphere exclusion volume correlation described by the fundamental measure theory (FMT) as

$$F^\text{hs}[\{\rho_i (\boldsymbol{r})\}] = k_B T\int_{V} d \boldsymbol{r}\ \Phi_\textrm{FMT}(\{ n_\alpha(\boldsymbol{r})\})$$

where $ n_\alpha(\boldsymbol{r}) = \sum_i \int_{V} d \boldsymbol{s}\ \rho_i (\boldsymbol{s})\omega^{(\alpha)}_i(\boldsymbol{r}-\boldsymbol{s})$ are the weigthed densities given by the convolution with the weigth function $\omega^{(\alpha)}_i(\boldsymbol{r})$. The function $\Phi$ can be described using different formulations of the fundamental measure theory (FMT) as

- [x] **R**osenfeld **F**unctional (**RF**) - [Rosenfeld, Y., Phys. Rev. Lett. 63, 980–983 (1989)](https://link.aps.org/doi/10.1103/PhysRevLett.63.980)
- [x] **W**hite **B**ear version **I** (**WBI**) - [Yu, Y.-X. & Wu, J., J. Chem. Phys. 117, 10156–10164 (2002)](http://aip.scitation.org/doi/10.1063/1.1520530); [Roth, R., Evans, R., Lang, A. & Kahl, G., J. Phys. Condens. Matter 14, 12063–12078 (2002)](https://iopscience.iop.org/article/10.1088/0953-8984/14/46/313)
- [x] **W**hite **B**ear version **II** (**WBII**) - [Hansen-Goos, H. & Roth, R. J., Phys. Condens. Matter 18, 8413–8425 (2006)](https://iopscience.iop.org/article/10.1088/0953-8984/18/37/002)

where [x] represents the implemented functionals.

The thermodynamic equilibrium is given by the functional derivative of the grand potential in the form 

$$ \frac{\partial \Omega}{\partial \rho_i(\boldsymbol{r})} = k_B T \ln(\rho_i(\boldsymbol{r}) \Lambda_i^3) + \frac{\partial F^\text{hs}[\rho_j]}{\partial \rho_i(\boldsymbol{r})}  +V_i^{(\text{ext})}(\boldsymbol{r})-\mu_i = 0$$

# Examples

On the folder 'examples' you can find different applications of the FMT. 

## Hard-Sphere near a Hardwall
|![Figure1](https://github.com/elvissoares/PyFMT/blob/main/examples/hardwall-eta0.4257.png)|![Figure2](https://github.com/elvissoares/PyFMT/blob/main/examples/hardwall-eta0.4783.png)|
|:--:|:--:|
| <b>Fig.1 - The density profiles of a pure hard-sphere fluid at a planar hard wall with bulk packing fraction of η = 0.4257. </b>| <b>Fig.2 - The density profiles of a pure hard-sphere fluid at a planar hard wall with bulk packing fraction of η = 0.4783. </b>|

## Hard-Sphere Mixture near a Hardwall

|![Figure3](https://github.com/elvissoares/PyFMT/blob/main/examples/hardwall-mixture-eta%3D0.39-x1%3D0.25-ratio3-small.png)|![Figure4](https://github.com/elvissoares/PyFMT/blob/main/examples/hardwall-mixture-eta%3D0.39-x1%3D0.25-ratio3-big.png)|
|:--:|:--:|
| <b>Fig.3 - Density profiles at a planar hard wall of the *small spheres* of a binary mixture with size ratio σ_b = 3σ_s and packing η = 0.39 and x1 = 0.25. </b>| <b>Fig.4 - Density profiles at a planar hard wall of the *big spheres* of a binary mixture with size ratio σ_b = 3σ_s and packing η = 0.39 and x1 = 0.25. </b>|

## Hard-Sphere Radial Distribution Function

|![Figure5](https://github.com/elvissoares/PyFMT/blob/main/examples/radialdistribution-rhob0.2.png)|![Figure6](https://github.com/elvissoares/PyFMT/blob/main/examples/radialdistribution-rhob0.9.png)|
|:--:|:--:|
| <b>Fig.5 - The radial distribution function of a pure hard-sphere fluid with bulk density ρ_b = 0.2. </b>| <b>Fig.6 - The radial distribution function of a pure hard-sphere fluid with bulk density ρ_b = 0.9.  </b>|