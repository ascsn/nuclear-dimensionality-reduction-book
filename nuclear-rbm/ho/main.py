import numpy as np
from scipy import special
import scipy
import matplotlib.pyplot as plt
import time
from scipy.optimize import minimize
from scipy.linalg import eigh

def exact_wavefunction(n,alpha):
    """
    Exact wavefunction, from 
    https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#Hamiltonian_and_energy_eigenstates
    modified for our notation.

    Parameters
    ----------
    n : int
        Order of the solution.
    alpha : float
        The HO parameter.

    Returns
    -------
    wf : function
        The exact wavefunction.

    """
    norm = 1/np.sqrt(2.**n * np.math.factorial(n))*(np.sqrt(alpha)/np.pi)**(1/4)
    herm = special.hermite(n)
    def wf(x):
        return norm*np.exp(-np.sqrt(alpha)*x**2/2)*herm(x)
    return wf

def exact_eigenvalue(n,alpha):
    """
    Exact eigenvalue $\lambda$ of the ODE
    
    $$ -f''(x) + alpha x^2f(x) = \lambda f(x) $$
    
    We assume throughout that $m = \hbar = 1$. Then, $\lambda$ is related to the
    energy as $\lambda = 2E$.

    Parameters
    ----------
    n : int
        Order of the solution.
    alpha : float
        The HO parameter.

    Returns
    -------
    float
        The eigenvalue.

    """
    
    return np.sqrt(alpha)*(n+0.5)

class FiniteDifferenceSolver:
    def __init__(self,alpha,regularGrid,h):
        self.alpha = alpha
        self.grid = regularGrid
        
        self.H = self.make_Hamiltonian_matrix()
        
    def make_kinetic_term(self):
        size = len(self.grid)
        offDiag = np.zeros(size)
        offDiag[1] = 1
        
        #Kinetic term
        H = -1*(-2*np.identity(size) + scipy.linalg.toeplitz(offDiag))/h**2
        return H
        
    def make_potential_term(self):
        H = np.diag(self.alpha*self.grid**2)
        return H
        
    def make_Hamiltonian_matrix(self):
        return self.make_kinetic_term() + self.make_potential_term()
    
    def solve(self):
        t0 = time.time()
        evals, evects = np.linalg.eigh(self.H)
                    
        t1 = time.time()
        
        solveTime = t1 - t0
        return evals, evects, solveTime
    
class StandardRBM:
    def __init__(self,psiMat,kineticTerm,potentialTerm):
        #potentialTerm should not be multiplied by alpha
        self.psiMat = psiMat
        self.kineticTerm = kineticTerm
        self.potentialTerm = potentialTerm
        self.matrixDim = psiMat.shape[0]
        
        self.kinProj, self.potProj = self.compute_projections()
        self.overlaps = self.compute_overlaps()
        
    def compute_projections(self):
        kinProj = np.zeros(2*(self.matrixDim,))
        potProj = np.zeros(2*(self.matrixDim,))
        for i in range(self.matrixDim):
            for j in range(i,self.matrixDim):
                kinProj[i,j] = self.psiMat[i] @ self.kineticTerm @ self.psiMat[j]
                kinProj[j,i] = kinProj[i,j]
                
                potProj[i,j] = self.psiMat[i] @ self.potentialTerm @ self.psiMat[j]
                potProj[j,i] = potProj[i,j]
                    
        return kinProj, potProj
    
    def compute_overlaps(self):
        overlapMatrix = np.zeros(2*(self.matrixDim,))
        for i in range(self.matrixDim):
            for j in range(i,self.matrixDim):
                overlapMatrix[i,j] = self.psiMat[i] @ self.psiMat[j]
                overlapMatrix[j,i] = overlapMatrix[i,j]
        return overlapMatrix
        
    def solve(self,alpha):
        mat = self.kinProj + alpha*self.potProj
            
        vals, vecs = eigh(mat,b=self.overlaps)
        
        a = vecs[:,0]
        
        wf = a @ self.psiMat
        wf /= np.linalg.norm(wf)
        
        return vals[0], wf, a
    
class VariationalRBM:
    def __init__(self,psiMat,h,regularGrid):
        self.psiMat = psiMat
        self.h = h
        self.grid = regularGrid
        
        self.m1, self.m2, self.m3 = self.matrix_components()
        
    def matrix_components(self):
        nComps, vecSize = self.psiMat.shape
        
        offDiag = np.zeros(vecSize)
        offDiag[1] = 1
        
        kinTerm = -1*(-2*np.identity(vecSize) + scipy.linalg.toeplitz(offDiag))/self.h**2
        
        m1 = np.zeros((nComps,nComps))
        m2 = np.zeros((nComps,nComps))
        m3 = np.zeros((nComps,nComps))
        for i in range(nComps):
            for j in range(i+1):
                #Kinetic energy
                m1[i,j] = self.psiMat[i] @ kinTerm @ self.psiMat[j]
                m1[j,i] = m1[i,j]
                
                #Potential energy
                m2[i,j] = np.sum(self.psiMat[i]*self.psiMat[j]*self.grid**2)
                m2[j,i] = m2[i,j]
                
                #Normalization
                m3[i,j] = self.psiMat[i] @ self.psiMat[j]
                m3[j,i] = m3[i,j]
        return m1, m2, m3
    
    def energy(self,a,alpha):
        return a @ (self.m1+alpha*self.m2) @ a/(a @ self.m3 @ a)

    def energy_grad(self,a,alpha):
        f1Eval = a @ (self.m1+alpha*self.m2) @ a
        f2Eval = a @ self.m3 @ a
        return 2*(self.m1+alpha*self.m2) @ a/f2Eval - 2*f1Eval/f2Eval**2 * self.m3 @ a
    
    def solve(self,alpha,a0):
        t0 = time.time()
        opt = minimize(self.energy,a0,args=(alpha,),jac=self.energy_grad)
        
        #Normalize total wavefunction
        a = opt.x
        wf = a@T
        wf /= np.linalg.norm(wf)
        
        eigenvalue = self.energy(a,alpha)
        
        t1 = time.time()
        
        return eigenvalue, wf, t1-t0, a
    
class InferHamiltonian:
    def __init__(self,psiMat,arrOfHamiltonians,alphaVals):
        self.psiMat = psiMat
        self.arrOfHamiltonians = arrOfHamiltonians
        self.alphaVals = alphaVals
        self.matrixDim = len(alphaVals)
        
        self.projections = self.compute_projections()
        self.overlaps = self.compute_overlaps()
        
    def compute_projections(self):
        projections = np.zeros(3*(self.matrixDim,))
        for k in range(self.matrixDim):
            for i in range(self.matrixDim):
                for j in range(i,self.matrixDim):
                    projections[k,i,j] = self.psiMat[i] @ self.arrOfHamiltonians[k] @ self.psiMat[j]
                    projections[k,j,i] = projections[k,i,j]
                    
        return projections
    
    def compute_overlaps(self):
        overlapMatrix = np.zeros(2*(self.matrixDim,))
        for i in range(self.matrixDim):
            for j in range(i,self.matrixDim):
                overlapMatrix[i,j] = self.psiMat[i] @ self.psiMat[j]
                overlapMatrix[j,i] = overlapMatrix[i,j]
        return overlapMatrix
    
    def fit_matrices(self,polyDegree=1):
        #Note that minCoeffs is indexed backwards - the first element is
        #the coefficient of the highest degree term in the polynomial
        minCoeffs = np.zeros((polyDegree+1,self.matrixDim,self.matrixDim,),)
        for i in range(self.matrixDim):
            for j in range(i,self.matrixDim):
                minCoeffs[:,i,j] = np.polyfit(self.alphaVals,
                                              self.projections[:,i,j],
                                              polyDegree)
                minCoeffs[:,j,i] = minCoeffs[:,i,j]
        return minCoeffs
        
    def solve(self,alpha,coeffs):
        mat = np.zeros(2*(self.matrixDim,))
        maxPolyDim = coeffs.shape[0] - 1
        for (coeffIter,coeff) in enumerate(coeffs):
            mat += alpha**(maxPolyDim - coeffIter) * coeff
            
        vals, vecs = eigh(mat,b=self.overlaps)
        
        a = vecs[:,0]
        
        wf = a @ self.psiMat
        wf /= np.linalg.norm(wf)
        
        return vals[0], wf, a

#First define global variables
h = 10**(-2) ### grid spacing for domain (Warning around 10**(-3) it starts to get slow).
### HO global parameters 
n = 0 # principle quantum number to solve in HO
# define the domain boundaries
xLims = [-10,10]
grid = np.arange(xLims[0],xLims[1]+h,h)
m = len(grid)
print('Number of grid points: ',m)

# Select alpha values to use to solve SE exactly.
alphaVals = [.5,5,10,15]
#Here, we choose values of alpha to solve exactly. The solutions are the basis functions
#to be used with the RBM
T = np.zeros((len(alphaVals),m)) 

listOfHamiltonians = []

for i,alpha in enumerate(alphaVals):
    fdSolver = FiniteDifferenceSolver(alpha,grid,h)
    evals, evecs, runTime = fdSolver.solve()
    print('FD time: %.3e s'%runTime)
    listOfHamiltonians.append(fdSolver.make_Hamiltonian_matrix())
    kinTerm = fdSolver.make_kinetic_term()
    potTerm = fdSolver.make_potential_term()/alpha
    T[i] = evecs[:,n]
    
a0 = np.array([1.,1,1,1])
alphaTest = 0.3

fig, ax = plt.subplots()
fdSol = FiniteDifferenceSolver(alphaTest,grid,h)
t0 = time.time()
evals, evecs, _ = fdSol.solve()
t1 = time.time()
ax.plot(grid,evecs[:,0],label='FD')
print('Finite Difference time: %.3e s'%(t1-t0))

rbm = StandardRBM(T,kinTerm,potTerm)
t0 = time.time()
val, wf, a = rbm.solve(alphaTest)
t1 = time.time()
ax.plot(grid,wf,label='RBM')
print('Standard RBM time: %.3e s'%(t1-t0))

rbm = VariationalRBM(T,h,grid)
eigenvalue, wfRBM, runTime, a = rbm.solve(alphaTest,a0)
ax.plot(grid,wfRBM,label='Variational RBM',ls='--')
print('Variational RBM time: %.3e s'%runTime)

inferObj = InferHamiltonian(T,listOfHamiltonians,alphaVals)
coeffs = inferObj.fit_matrices()
t0 = time.time()
val, wf, a = inferObj.solve(alphaTest,coeffs)
t1 = time.time()
ax.plot(grid,wf,label='Infer Matrix',ls=':')
print('Inferred Matrix time: %.3e s'%(t1-t0))


ax.legend()

# exactEig = 2*exact_eigenvalue(n,alphaTest) #Convert from energy to eigenvalue
# print('Eigenvalue relative difference: %.3e'%((eigenvalue-exactEig)/exactEig))

