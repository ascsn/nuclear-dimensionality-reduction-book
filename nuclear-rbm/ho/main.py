import numpy as np
from scipy import special
import scipy
import matplotlib.pyplot as plt
import time
from scipy.optimize import minimize

def make_wavefunction(n,k):
    '''
    Definition of exact HO wavefunction taken from Zettili page 240.
    
    Parameters
    ----------
    n : TYPE
        principle quantum number for SE equation
    k : TYPE
        harmonic oscillator parameter (mass*omega)^2 from the potential term \mu^2 \omega^{2} x^{2} in the SE.

    Returns
    -------
    wf : function
        1-d wavefunction for the 1d-harmonic oscillator as a function of position x. 
    '''
    herm = special.hermite(n)
    def wf(x):
        result = (1/np.sqrt(np.sqrt(np.pi)*2**(n)*np.math.factorial(n)*(k)**(.25)))*np.exp(-x**(2)*np.sqrt(k)/2)*herm(x*k**(.25))
        return result
    return wf

def get_eigenvalue(n,mass,alpha):
    '''
    Exact eigenvalues of the HO equation. -f''(x) + k x^2 f(x) = 2 m E f(x)
    Lambda = 2 m E  
    E = (n + .5) \omega 
    \alpha = m^2 \omega^2 
    \omega = \sqrt{alpha/m^2}
    Parameters
    ----------
    n : float
        principle quantum number. float integers
    omega : float
        oscillator frequency.

    Returns
    -------
    float
        oscillator energy 2mE. 
    '''
    return 2*mass*(.5 + n)*np.sqrt(alpha/mass**2)

def V(x,alpha):
    '''
    1-d harmonic Oscillator potential

    Parameters
    ----------
    x : float or nd array
        position.
    alpha : float
        oscillator length parameter.

    Returns
    -------
    float or ndarray
        value of potential evaluated at x.

    '''
    return alpha*x**2

def construct_H(V,grid,mass,alpha):
    '''
    Uses 2nd order finite difference scheme to construct a discretized differential H operator
    Note: mass is fixed to 1 for this demo.

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.
        
    alpha : TYPE
        oscillator parameter alpha used in V(x) = alpha*x**2.

    Returns
    -------
    H : ndarray
        2d numpy array
    '''
    dim = len(grid)
    off_diag = np.zeros(dim)
    off_diag[1] = 1
    H = -1*(-2*np.identity(dim) + scipy.linalg.toeplitz(off_diag))/(mass*h**2) + np.diag(V(grid,alpha))
    return H

def get_original_solution(H,grid,h):
    '''
    Parameters
    ----------
    H : 2d ndarray
        Hamiltonian Matrix.
    grid : ndarray
        Discretized 1d domain.
    h : float
        grid spacing.

    Returns
    -------
    evals : ndarray
        returns nd array of eigenvalues of H. 
    evects : ndarray
        returns ndarray of eigenvectors of H.
    Eigenvalues and eigenvectors are ordered in ascending order. 
    '''
    evals,evects = np.linalg.eigh(H)
    evects = evects.T
    for i,evect in enumerate(evects):
        #norm = np.sqrt(1/sci.integrate.simpson(evect*evect,grid))
        norm = 1/(np.linalg.norm(np.dot(evect,evect)))
        evects[i] = evects[i]*norm
    return evals,evects

def matrix_component_1(psi1,psi2):
    #Derivative bit
    dim = len(psi1)
    off_diag = np.zeros(dim)
    off_diag[1] = 1
    H = -1*(-2*np.identity(dim) + scipy.linalg.toeplitz(off_diag))/(mass*h**2)
    
    return psi1@H@psi2

def matrix_component_2(psi1,psi2,grid):
    #The potential bit, without alpha
    return np.sum(psi1*psi2*grid**2)

def matrix_component_3(psi1,psi2):
    return psi1@psi2

def f1(a,alpha):
    return a@(m1+alpha*m2)@a

def f2(a):
    return a@m3@a

def energy(a,alpha):
    return a@(m1+alpha*m2)@a/(a@m3@a)

def energy_grad(a,alpha):
    f1Eval = f1(a,alpha)
    f2Eval = f2(a)
    return 2*(m1+alpha*m2)@a/f2Eval - 2*f1Eval/f2Eval**2 * m3@a

def gradient_descent(a0,alpha,params={'maxIters':1000,'eta':0.1}):
    enegVals = np.zeros(params['maxIters'])
    a = a0.copy()
    for i in range(params['maxIters']):
        enegVals[i] = energy(a,alpha)
        grad = energy_grad(a,alpha)
        a -= params['eta']*grad
        
    return a, enegVals

#First define global variables
h = 10**(-2) ### grid spacing for domain (Warning around 10**(-3) it starts to get slow).
### HO global parameters 
n = 0 # principle quantum number to solve in HO
mass = 1.0 # mass for the HO system
# define the domain boundaries
x_a = -10 # left boundary 
x_b = 10 # right boundary 
x_array = np.arange(x_a,x_b+h,h)
m = len(x_array) 
print('Number of grid points: ',m)

# Select alpha values to use to solve SE exactly.
alpha_vals = [.5,5,10,15]  #Here, we choose 3 values of alpha to solve exactly. This results in 3 basis functions
# initialize solution arrays. T is the matrix that will hold wavefunction solutions. 
# T has the form T[i][j], i = alpha, j = solution components
T = np.zeros((len(alpha_vals),m)) 
# T_evals holds the eigenvalues for each evect in T. 
T_evals = np.zeros(len(alpha_vals))

for i,alpha_sample in enumerate(alpha_vals):
    H = construct_H(V,x_array,mass,alpha_sample) # construct the Hamiltonian matrix for given alpha_sample.
    evals, evects = get_original_solution(H,x_array,h) # solve the system for evals and evects.
    T[i] = evects[n]/np.linalg.norm(evects[n]) # assign the nth evect to solution array T
    T_evals[i] = evals[n] # assign the nth eigenvalue to the eigenvalue array T_eval.
    print(f'alpha = {alpha_sample}, lambda = {evals[n]}')
    print(f'alpha = {alpha_sample}, exact lambda = {get_eigenvalue(n,mass,alpha_sample)}')

#%%
m1 = np.zeros((len(alpha_vals),len(alpha_vals)))
m2 = np.zeros(m1.shape)
m3 = np.zeros(m1.shape)
for i in range(len(alpha_vals)):
    for j in range(i+1):
        m1[i,j] = matrix_component_1(T[i],T[j])
        m1[j,i] = m1[i,j]
        
        m2[i,j] = matrix_component_2(T[i],T[j],x_array)
        m2[j,i] = m2[i,j]
        
        m3[i,j] = matrix_component_3(T[i],T[j])
        m3[j,i] = m3[i,j]
        
     
a0 = np.array([1.33380471, -1.50572164, 2.22572807, -1.06054276])
# a0 = np.array([1.,0,0,0])
alphaTest = 0.3

opt = minimize(energy,a0,args=(alphaTest,))
a = opt.x

print('My energy: %.3f'%opt.fun)
print('Exact eigenvalue: %.3f'%get_eigenvalue(0,1,alphaTest))
# t0 = time.time()
# a, enegVsIter = gradient_descent(a0,alphaTest)
# t1 = time.time()
# print('RBM time: %.3e s'%(t1-t0))

# fig, ax = plt.subplots()
# ax.plot(enegVsIter)

# H = construct_H(V,x_array,mass,alphaTest) # construct the Hamiltonian matrix for given alpha_sample.
# t0 = time.time()
# evals, evects = get_original_solution(H,x_array,h) # solve the system for evals and evects.
# t1 = time.time()
# print('Finite difference time: %.3e s'%(t1-t0))
# ax.axhline(evals[0])
# ax.axhline(get_eigenvalue(0,1,alphaTest),ls='--',color='red')

newWavefunction = a@T
print(newWavefunction@newWavefunction)
# newWavefunction /= newWavefunction@newWavefunction

fig, ax = plt.subplots()
ax.plot(x_array,newWavefunction)

psi_real = make_wavefunction(0,alphaTest)
psiArr = psi_real(x_array)
psiArr /= np.sqrt(psiArr@psiArr)
print(psiArr@psiArr)

ax.plot(x_array,psiArr)

# print(energy_grad(np.array([1,0,0,0]),0.2))
# ## Checks to make sure numerical solutions are accurate:
# # Make plots comparing the numerical wavefunction to the exact wavefunction
# fig, ax = plt.subplots(1,2,figsize=(10,3))
# for i in range(len(alpha_vals)):
#     wf = make_wavefunction(n,alpha_vals[i])
#     wf_vals = wf(x_array)/np.linalg.norm(wf(x_array))
#     diff = np.abs(wf_vals - T[i])
#     ax[0].plot(x_array,np.abs(T[i]),label=r'$\alpha$ = '+str(alpha_vals[i]))
#     ax[1].semilogy(x_array,diff,label=r'$\alpha$='+str(alpha_vals[i]))
# ax[0].set_title('Numerical solutions')
# ax[1].set_title('Comparison of exact and numerical sols')
# ax[0].legend()
# ax[1].legend()
# ax[1].set_ylabel(r'$||\phi-\phi^{*}||$')
# fig.tight_layout()
# plt.show()

# # print error of eigenvalue 
# for i in range(len(alpha_vals)):
#     diff = np.abs(get_eigenvalue(n,mass,alpha_vals[i]) - T_evals[i])
#     print(f'|lambda - lambda_exact| = {diff}')

