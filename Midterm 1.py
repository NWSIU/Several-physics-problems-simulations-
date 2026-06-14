### midterm 1 (1404/09/12)
from pylab import*
from numpy import*
from matplotlib import*
from scipy import*

###functions
### creating A and B matrices with Crank-Nicholson method ###
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
    norm = sqrt(sum(abs(psi)**2) * dx)  ### normlization
    psi[:] = psi / norm
### defining a time dependent (harmonic oscillator) potential ###
def make_V(t_index=0, time_dep=True):  
    m = 1
    omega = 1.5
    x0 = i0*dx+velocity*(t_index*dt) ###time-dependent term (centre of potential)
    for i in range(N):
        V[i] = 0.5*m*(omega**2)*((i*dx-x0)**2)
###
### adding hamiltonian for calculating the mean energy
def H_psi(psi_in):
    Hpsi_local = zeros(N,complex)

    for i in range(1, N-1):
        lapx = (psi_in[i+1] - 2.0*psi_in[i] + psi_in[i-1]) / (dx*dx)
        Hpsi_local[i] = -lapx + V[i] * psi_in[i]
    # boundaries
    Hpsi_local[0] = Hpsi_local[1]
    Hpsi_local[N-1] = Hpsi_local[N-2]
    return Hpsi_local
    
###variables
### position gridding
N = 500
x_min = 0.0
x_max = 4.0
x = linspace(x_min,x_max,N)
dx = x[1]-x[0]

### time gridding
T = 1000
dt = 0.0005
time = arange(0,T*dt,dt)

velocity = 1         ### speed of the centre of potential
k = 30               # speed of the wavefunction
s = 0.05
i0 = int(N/4)

A = zeros((N,N),complex) 
B = zeros((N,N),complex)

psi = zeros(N,complex)
psi0 = zeros(N,complex)
V = zeros(N,float)
Hpsi = zeros(N,complex)

energy_vs_time = zeros(T,float)
overlap_vs_time = zeros(T,float)

### main
initial_psi()
H_psi(psi)
psi0[:] = psi.copy()          # store initial state
overlap_vs_time[0] = 1.0      # correct starting probability (background knowledge)

for n in range(T):
    make_V(t_index=n, time_dep=True)
    make_A()
    make_B()

    b = B.dot(psi)
    
    psi = linalg.solve(A,b)
    psi /= sqrt(sum(abs(psi)**2) * dx) ### normalization

    Hpsi = H_psi(psi)
    E_mean_t = real(sum(conjugate(psi) * Hpsi) * dx)
    energy_vs_time[n] = E_mean_t

    overlap = sum(conjugate(psi0) * psi) * dx
    overlap_vs_time[n] = abs(overlap)**2
   
    clf()
    plot(x, V / max(V) * max(abs(psi)**2), 'r', label="V(x) (scaled)") 
    plot(x, abs(psi)**2, 'b', label="|psi|^2")
    title(f"P(x,t) and V(x) at step {n}")
    xlabel("x")
    ylabel("amplitude of ")
    pause(0.001)

print("Final <E> =", energy_vs_time[-1])
print("Final overlap probability =", overlap_vs_time[-1])

### potential well
subplot(2, 2, 1)
plot(x, V, label='V(x)')
title("potential well V(x)")
xlabel("x")
ylabel("V(x)")
grid(True)
legend()

### mean energy
subplot(2, 1, 2)
time_array = arange(T) * dt
plot(time_array, energy_vs_time)
title("mean energy ⟨E(t)⟩")
xlabel("t")
ylabel("⟨E⟩")
grid(True)

tight_layout()
show()

### overlap with respect to time
figure(figsize=(6,4))
plot(time_array, overlap_vs_time)
title("overlap probability P_overlap(t)")
xlabel("t")
ylabel("P_overlap")
ylim(-0.02, 1.20)
grid(True)
show()
