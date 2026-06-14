from pylab import *
from numpy import *
from matplotlib import *

### functions
def make_graph():
    for x in range(L):
        for y in range(L):
            i = x + y*L 

            if x+1<L:
                j = (x+1)+ y*L
                if uniform(0,1)<p:
                    li = k[i]
                    lj = k[j]
                    v[i][li] = j
                    v[j][lj] = i
                    k[i] += 1
                    k[j] += 1
            
            if y+1<L:
                j = x +(y+1)*L
                if uniform(0,1)<p:
                    li = k[i]
                    lj = k[j]
                    v[i][li] = j
                    v[j][lj] = i
                    k[i] += 1
                    k[j] += 1

###
def components():
    t = 0 
    while len(Q)>0:
        i = choice(Q)
        make_cluster(i,t)

        t += 1

###
def make_cluster(i,t):
    P = zeros(N,int)
    Pnew = zeros(N,int)

    n = 1
    P[0] = i
    comp[i] = t
    Q.remove(i)

    while n > 0:
        nnew = 0
        for l in range(n):
            ii = P[l]
            for ll in range(k[ii]):
                j = v[ii][ll]
                if comp[j] < 0:
                    comp[j] = t
                    Q.remove(j)
                    Pnew[nnew] = j
                    nnew += 1
        P[:] = Pnew[:]
        n = nnew

###
def percolating_cluster():
    mmax = 0
    cmax = 0
    nc = max(comp)+1

    for c in range(nc):
        mc = comp.count(c)
        if mc > mmax:
            mmax = mc
            cmax = c
    
    return cmax

###
def mass_R(seed,cmax):
    x0 = seed % L
    y0 = seed // L

    Rlist = []
    Mlist = []

    Rmax = L//2

    for R in range(1,Rmax):
        m = 0
        for i in range(N):
            if comp[i] == cmax:
                x = i % L
                y = i // L
                r = sqrt((x-x0)**2 + (y-y0)**2)
                if r <= R:
                    m += 1
        if m > 0:
            Rlist.append(R)
            Mlist.append(m)

    return array(Rlist),array(Mlist)

### variables           
L = 50
N = L*L          
p = 0.59

k = zeros(N,int)
v = zeros((N,4),int)

Q = [i for i in range(N)]
comp = [-1 for i in range(N)]

### main
make_graph()

components()

cmax = percolating_cluster()

cluster_sites = [i for i in range(N) if comp[i] == cmax]
seed = choice(cluster_sites)

R,M_R = mass_R(seed,cmax)

logR = log(R)
logM = log(M_R)

plot(logR,logM,'o', label='data')

coef = polyfit(logR,logM,1)
Df = coef[0]

plot(logR,coef[0]*logR +coef[1],'-', label='fit')

xlabel('log R')
ylabel('log M(R)')
legend()
title('Fractal dimension D = %.3f' % Df)

print("Fractal dimension D =", Df)
show()
