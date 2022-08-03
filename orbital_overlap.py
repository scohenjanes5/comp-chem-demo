#copied from https://joshuagoings.com/assets/integrals.pdf
import numpy as np
from scipy.special import factorial2, hyp1f1
class BasisFunction(object):
    ''' 
    A class that contains all our basis function data
    Attributes:
    origin: array/list containing the coordinates of the Gaussian origin
    shell: tuple of angular momentum
    exps: list of primitive Gaussian exponents
    coefs: list of primitive Gaussian coefficients
    norm: list of normalization factors for Gaussian primitives
    '''
    def __init__(self,origin=[0.0,0.0,0.0],shell=(0,0,0),exps=[],coefs=[]):
        self.origin = np.asarray(origin)
        self.shell = shell
        self.exps = exps
        self.coefs = coefs
        self.norm = None
        self.normalize()
    def normalize(self):
        ''' 
        Routine to normalize the basis functions, in case they
        do not integrate to unity.
        '''
        l,m,n = self.shell
        L = l+m+n
        pi_power = np.power(np.pi, 1.5)
        fact2_pdct = factorial2(2*l - 1)*factorial2(2*m - 1)*factorial2(2*n - 1)
        # self.norm is a list of length equal to number primitives
        # normalize primitives first (PGBFs)
        self.norm = np.sqrt(np.power(2, 2*L + 1.5)*
            np.power(self.exps, L+1.5)/ fact2_pdct / pi_power)
        # now normalize the contracted basis functions (CGBFs)
        # Eq. 1.44 of Valeev integral whitepaper
        prefactor = pi_power * fact2_pdct / np.power(2.0, L)
        N = 0.0
        num_exps = len(self.exps)
        for i in range(num_exps):
            for j in range(num_exps):
                N += self.norm[i]*self.norm[j]*self.coefs[i]*self.coefs[j]/\
                    np.power(self.exps[i] + self.exps[j], L+1.5)
        N *= prefactor
        N = np.power(N, -0.5)
        for i in range(num_exps):
            self.coefs[i] *= N

def E(i,j,t,Qx,a,b):
    ''' 
    Recursive definition of Hermite Gaussian coefficients.
    Returns a float. (1d component of overlap of 2 primitives)
    a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
    i,j: orbital angular momentum number on Gaussian 'a' and 'b'
    t: number nodes in Hermite (depends on type of integral,
    e.g. always zero for overlap integrals)
    Qx: distance between origins of Gaussian 'a' and 'b'
    '''
    p = a + b
    q = a*b/p
    
    if (t < 0) or (t > (i + j)):
        # out of bounds for t
        return 0.0
    elif i == j == t == 0:
        # base case
        return np.exp(-q*Qx*Qx) # K_AB
    elif j == 0:
        # decrement index i
        Idown = E(i-1,j,t,Qx,a,b)
        Idown_tup = E(i-1,j,t+1,Qx,a,b)
        Idown_tdown = E(i-1,j,t-1,Qx,a,b)
        
        return (1/(2*p)) * Idown_tdown - \
        (q*Qx/a) * Idown + \
        (t+1) * Idown_tup
    else:
        # decrement index j
        Jdown = E(i,j-1,t,Qx,a,b)
        Jdown_tup = E(i,j-1,t+1,Qx,a,b)
        Jdown_tdown = E(i,j-1,t-1,Qx,a,b)
        
        return (1/(2*p)) * Jdown_tdown + \
        (q*Qx/b) * Jdown + \
        (t+1) * Jdown_tup

def overlap(a,lmn1,A,b,lmn2,B):
    ''' 
    Evaluates overlap integral between two Gaussians primitives
    Returns a float.
    a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
    lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
    for Gaussian 'a'
    lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
    A: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
    B: list containing origin of Gaussian 'b'
    '''
    # l1,m1,n1 = lmn1 # shell angular momentum on Gaussian 'a'
    # l2,m2,n2 = lmn2 # shell angular momentum on Gaussian 'b'
    Sproduct = 1
    for i in range(3):
        Sproduct *= E(lmn1[i], lmn2[i], 0, A[0]-B[0], a, b)

    # S1 = E(l1,l2,0,A[0]-B[0],a,b) # X
    # S2 = E(m1,m2,0,A[1]-B[1],a,b) # Y
    # S3 = E(n1,n2,0,A[2]-B[2],a,b) # Z
    return Sproduct*np.power(np.pi/(a+b),1.5)

def S(a,b):
    '''
    Evaluates overlap between two contracted Gaussians (multiple primitives)
    Returns float.
    Arguments:
    a: contracted Gaussian 'a', BasisFunction object
    b: contracted Gaussian 'b', BasisFunction object
    '''
    s = 0.0
    for i, ca in enumerate(a.coefs):
        for j, cb in enumerate(b.coefs):
            s += a.norm[i]*b.norm[j]*ca*cb*\
            overlap(a.exps[i],a.shell,a.origin,
            b.exps[j],b.shell,b.origin)
    return s

# KE Part

def kinetic(a,lmn1,A,b,lmn2,B):
    '''
    Evaluates kinetic energy integral between two Gaussian primitives
    Returns a float.
    KE part can be written as modified overlap integrals
    a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
    lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
    for Gaussian 'a'
    lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
    A: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
    B: list containing origin of Gaussian 'b'
    '''
    l1,m1,n1 = lmn1
    l2,m2,n2 = lmn2
    term0 = b*(2*(l2+m2+n2)+3)*overlap(a,lmn1,A,b,lmn2,B)

    term1 = -2*np.power(b,2)*\
        (overlap(a,lmn1,A,b,(l2+2,m2,n2),B) +
        overlap(a,lmn1,A,b,(l2,m2+2,n2),B) +
        overlap(a,lmn1,A,b,(l2,m2,n2+2),B))

    term2 = -0.5*(l2*(l2-1)*overlap(a,lmn1,A,b,(l2-2,m2,n2),B) +
        m2*(m2-1)*overlap(a,lmn1,A,b,(l2,m2-2,n2),B) +
        n2*(n2-1)*overlap(a,lmn1,A,b,(l2,m2,n2-2),B))

    return term0+term1+term2
    
def T(a,b):
    '''
    Evaluates kinetic energy between two contracted Gaussian orbitals
    Returns float.
    (all combinations of primitives)
    Arguments:
    a: contracted Gaussian 'a', BasisFunction object
    b: contracted Gaussian 'b', BasisFunction object
    '''
    t = 0.0
    for i, ca in enumerate(a.coefs):
        for j, cb in enumerate(b.coefs):
            t += a.norm[i]*b.norm[j]*ca*cb*\
                kinetic(a.exps[i],a.shell,a.origin,
                b.exps[j],b.shell,b.origin)
    return t

def gaussian_product_center(a,A,b,B):
    return (a*A+b*B)/(a+b)

def boys(n,T):
    return hyp1f1(n+0.5, n+1.5, -T) / (2.0*n + 1.0) #equation 27 where 1F1 is hyp1f1 here

##Coulombic part
def R(t,u,v,n,p,PCx,PCy,PCz,RPC):
    ''' 
    Returns the Coulomb auxiliary Hermite integrals
    Returns a float.
    Arguments:
    t,u,v: order of Coulomb Hermite derivative in x,y,z
    (see defs in Helgaker and Taylor)
    n: order of Boys function
    PCx,y,z: Cartesian vector distance between Gaussian
    composite center P and nuclear center C
    RPC: Distance between P and C
    '''
    T = p*RPC*RPC
    val = 0.0
    if t == u == v == 0:
        val += np.power(-2*p,n)*boys(n,T) #equation 25
    elif t == u == 0: #i.e. z
        if v > 1:
            val += (v-1)*R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC)
        val += PCz*R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)
    elif t == 0: #i.e. y
        if u > 1:
            val += (u-1)*R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC)
        val += PCy*R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)
    else: #i.e. x
        if t > 1:
            val += (t-1)*R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC)
        val += PCx*R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)
    return val

def nuclear_attraction(a,lmn1,A,b,lmn2,B,C):
    ''' 
    Evaluates nuclear attraction integral between two Gaussian primitives
    Returns a float.
    a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
    b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
    lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
    for Gaussian 'a'
    lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
    A: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
    B: list containing origin of Gaussian 'b'
    C: list containing origin of nuclear center 'C'
    '''
    l1,m1,n1 = lmn1
    l2,m2,n2 = lmn2
    p = a + b
    P = gaussian_product_center(a,A,b,B) # Gaussian composite center
    RPC = np.linalg.norm(P-C) #dist between centers
    
    val = 0.0
    for t in range(l1+l2+1): #every
        for u in range(m1+m2+1): #combo
            for v in range(n1+n2+1): #represented in equation 29
                val += E(l1,l2,t,A[0]-B[0],a,b) * \
                    E(m1,m2,u,A[1]-B[1],a,b) * \
                    E(n1,n2,v,A[2]-B[2],a,b) * \
                    R(t,u,v,0,p,P[0]-C[0],P[1]-C[1],P[2]-C[2],RPC)
    val *= 2 * np.pi / p #constant in front of summation.
    return val

def V(a,b,C):
    '''
    Evaluates overlap between two contracted Gaussians
    Returns float.
    Arguments:
    a: contracted Gaussian 'a', BasisFunction object
    b: contracted Gaussian 'b', BasisFunction object
    C: center of nucleus
    '''
    v = 0.0
    for i, ca in enumerate(a.coefs): #every combo of primitives that make up an orbital.
        for j, cb in enumerate(b.coefs):
            v += a.norm[i]*b.norm[j]*ca*cb*\
                nuclear_attraction(a.exps[i],a.shell,a.origin,
                b.exps[j],b.shell,b.origin,C)
    return v


myOrigin = [1.0, 2.0, 3.0]
orig2 = [1.5, 2.0, 3.0]
myShell = (0,0,0) # p-orbitals would be (1,0,0) or (0,1,0) or (0,0,1), etc.
myExps = [3.42525091, 0.62391373, 0.16885540]
myCoefs = [0.15432897, 0.53532814, 0.44463454]
a = BasisFunction(origin=myOrigin,shell=myShell,exps=myExps,coefs=myCoefs)
b = BasisFunction(origin=orig2,shell=myShell,exps=myExps,coefs=myCoefs)
# print(S(a,a)) #1 if normallized properly
# print(T(a,a)) #0.760032 hartree is the electronic kinetic energy of the hydrogen atom described by the STO-3G basis set (compared to the exact value of 0.5 hartree).
# print(V(a,a,myOrigin))

print(hyp1f1(1,2,3))