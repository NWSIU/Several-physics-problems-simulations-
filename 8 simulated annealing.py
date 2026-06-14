from pylab import *
from numpy import *

### functions
def make_graph():
    k[:] = 0

    for x in range(L):
        for y in range(L):

            i = x + y*L
            hi[i] = h

            j = ((x+1)%L) + y*L
            li = k[i]
            lj = k[j]

            v[i][li] = j
            v[j][lj] = i

            k[i]+=1
            k[j]+=1

            Jbond = random_J()
            Jij[i][li] = Jbond
            Jij[j][lj] = Jbond

            j = x + ((y+1)%L)*L
            li = k[i]
            lj = k[j]

            v[i][li] = j
            v[j][lj] = i

            k[i]+=1
            k[j]+=1

            Jbond = random_J()
            Jij[i][li] = Jbond
            Jij[j][lj] = Jbond

###
def initial_MC():

    for i in range(N):
        s[i] = +1
        if (uniform()<0.5): 
            s[i] = -1

###
def update(i):
    global M,E

    de = 2*hi[i]*s[i]
    for li in range (k[i]):
        j = v[i][li]
        de = de + 2*Jij[i][li]*s[i]*s[j]

    if (uniform()<exp(-beta*de)):
        E += de
        M -= 2*s[i]
        s[i] = -s[i]

###
def compute_initial():
    e = 0
    m = 0
    for i in range(N):
        m +=s[i]
        e = e-hi[i]*s[i]
        for li in range(k[i]):
            j = v[i][li]
            if (j<i) :
                e = e - Jij[i][li]*s[i]*s[j]
    
    return e,m 

###
def random_J():

    if uniform() < 0.5:
        return -J
    else:
        return +J

###
def monte_carlo(T, teq):
    global beta, E, M

    beta = 1.0 / T

    for tmc in range(teq):

        O = arange(N)
        shuffle(O)

        for i in O:
            update(i)
###
def show_spins():

    lattice = reshape(s, (L, L))

    imshow(lattice, cmap='coolwarm')
    title("Final Spin Configuration")
    colorbar()
    show()

### variables
L = 50
N = L*L
K = 4

h = 0
J = 1
beta = 0.5

teq = 200

hi = zeros(N,float)
Jij = zeros((N,K),float)

s = zeros(N, int)

k = zeros(N,int)
v = zeros((N,K),int)

T0 = 3.0
Tmin = 0.05
deltaT = 0.05

### main
seed()

make_graph()
initial_MC()
E, M = compute_initial()

T = T0

energy = []
temp = []

while T>Tmin:

    print("Temp=", T)

    monte_carlo(T, teq)

    energy.append(E/N)
    temp.append(T)

    T = T - deltaT

# print("Annealing finished")
print("ground state energy per spin =", E / N)
figure()
plot(temp, energy, '-o')
xlabel("temp")
ylabel("energy/spin")
title("simulated annealing")
gca().invert_xaxis()
show()

show_spins()
