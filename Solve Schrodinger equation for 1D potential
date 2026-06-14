### solving schrodinger equation for a wave function with potential ###
from pylab import *

###function
def make_A():


    A[0,0] = 2*1j-2-2*V[0]
    A[0,1] = 1
    A[1,0] = 1
    A[N-1,N-1] =2*1j-2-2*V[N-1]
    A[N-2,N-1] = 1
    A[N-1,N-2] = 1
    for i in range(1,N-1):
        A[i,i-1] = 1
        A[i-1,i] = 1
        A[i,i] = 2*1j-2-2*V[i]
        A[i,i+1] = 1
        A[i+1,i] = 1

###
def make_B():

    B[0,0] = 2*1j+2+2*V[0]
    B[0,1] = -1
    B[1,0] = -1
    B[N-1,N-1] =2*1j+2+2*V[N-1]
    B[N-2,N-1] = -1
    B[N-1,N-2] = -1
    for i in range(1,N-1):
        B[i,i-1] = -1
        B[i-1,i] = -1
        B[i,i] = 2*1j+2+2*V[i]
        B[i,i+1] = -1
        B[i+1,i] = -1
        

###
def initial_psi():

    for i in range(N):
        psi[i] = exp((-((i*dx-i0*dx)**2)/(4*s**2))+1j*k*i*dx)

###
def make_V(): 
    for i in range(N):
        if (int(N/2) < i < int(3*N/4)):
            V[i] = 0.1

###variables
N = 500
T = 1000

k = 300
s = 0.05

dx = 0.002
dt = 2*dx*dx

i0 = int(N/4)

A = zeros((N,N),complex) 
B = zeros((N,N),complex)

###main
### initial condition ###
psi = zeros(N,complex)
V = zeros(N)

### solve the problem ###
make_V()
make_A()
make_B()
initial_psi()

for t in range(T):

    b = B.dot(psi)
    
    psi = linalg.solve(A,b)

    P=[]
    for i in range(N):
        P.append(psi[i]*conj(psi[i]))
    
    clf()
    plot(P)
    pause(0.01)
