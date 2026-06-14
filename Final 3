from numpy import *
from pylab import *

### functions
def make_graph():
    k[:] = 0
    v[:] = 0

    for x in range(L):
        for y in range(L):

            i = x+y*L

            if (x+1<L):
                j = (x+1)+y*L
                if (uniform()<p):
                    li = k[i]
                    lj = k[j]
                    v[i][li] = j
                    v[j][lj] = i
                    k[i] += 1
                    k[j] += 1

            if (y+1<L):
                j = x+(y+1)*L
                if (uniform()<p):
                    li = k[i]
                    lj = k[j]
                    v[i][li] = j
                    v[j][lj] = i
                    k[i] += 1
                    k[j] += 1

###
def components():

    Q[:] = [i for i in range(N)]
    for i in range(N): comp[i] = -1

    t = 0
    while(len(Q)>0):
        i = choice(Q)
        make_cluster(i,t)
        t += 1

###
def make_cluster(i,t):

    P = [i]
    comp[i] = t
    Q.remove(i)

    while len(P)>0:

        Pnew = []

        for ii in P:
            for ll in range(k[ii]):
                j = v[ii][ll]

                if(comp[j]<0):
                    comp[j] = t
                    Q.remove(j)
                    Pnew.append(j)

        P = Pnew

###
def display():

    A = zeros((L,L))
    mmax = 0
    nc = max(comp)+1

    for c in range (nc):
        mc = comp.count(c)
        if (mc>mmax):
            mmax = mc
            cmax = c

    for x in range(L):
        for y in range(L):
            i = x+y*L
            if(comp[i]==cmax):
                A[x,y] = -1

    imshow(A,origin='lower')
    title("largest cluster")
    show()

###
def make_G():

    nc = max(comp)+1
    sizes = [comp.count(c) for c in range(nc)]
    cmax = argmax(sizes)

    G = []
    for i in range(N):
        if comp[i]==cmax:
            G.append(i)

    return G

###
def make_G_neighbors():

    G = make_G()
    NG = len(G)

    map_index = {}
    for i in range(NG):
        map_index[G[i]] = i

    G_neigh = [[] for _ in range(NG)]

    for site in G:
        i = map_index[site]

        for ll in range(k[site]):
            j = v[site][ll]
            if j in map_index:
                G_neigh[i].append(map_index[j])

    return G_neigh, NG

###
def initial_MC():
    for i in range(NG):
        s[i] = +1
        if uniform()<0.5:
            s[i] = -1

###
def delta_energy(i):

    de = 0
    for j in G_neighbors[i]:
        de += s[j]

    return 2*J*s[i]*de

###
def update():

    for _ in range(NG):

        i = randint(NG)
        dE = delta_energy(i)

        if dE <= 0 or uniform()<exp(-dE/T):
            s[i] = -s[i]

###
def run_MC():

    initial_MC()

    for _ in range(teq):
        update()

    Mlist = []

    for _ in range(tsample):
        update()
        Mlist.append(abs(sum(s)/NG))

    Mlist = array(Mlist)

    Mmean = mean(Mlist)
    Merror  = sqrt(var(Mlist)/len(Mlist))

    return Mmean, Merror

### variables
L = 40
N = L*L
p = 0.5
J = 1.0

k = zeros(N,int)
v = zeros((N,4),int)

Q = [i for i in range(N)]
comp = [-1 for i in range(N)]

teq = 1000
tsample = 2000

s = zeros(N)

### main
seed()

make_graph()
components()

display()

G_neighbors, NG = make_G_neighbors()

Tlist = linspace(0.8,4.0,15)

Mvals = []
Merror  = []

for T in Tlist:
    M,error = run_MC()

    print(f"T = {T:.2f}, <M> = {M:.3f} ± {error:.4f}")

    Mvals.append(M)
    Merror.append(error)

errorbar(Tlist,Mvals,yerr=Merror,fmt='o-')
xlabel("temperature")
ylabel("mean magnetization")
show()
