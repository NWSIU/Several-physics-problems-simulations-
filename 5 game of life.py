from pylab import*
from matplotlib import*

### functions
def count_neighbors(i,j):
    cn = (lattice[i-1,j]+lattice[i-1,j-1]
    +lattice[i-1,j+1]+lattice[i+1,j]
    +lattice[i+1,j-1]+lattice[i+1,j+1]
    +lattice[i,j-1]+lattice[i,j+1])
    return cn

### variables
L = 50
t = 0
T = 170

### main
lattice = zeros((L,L))
### initial condition 
lattice[5,5] = 1
lattice[5,6] = 1
lattice[5,7] = 1
lattice[4,7] = 1
lattice[3,6] = 1

lattice[13,14] = 1
lattice[14,15] = 1
lattice[15,13] = 1
lattice[15,14] = 1
lattice[15,15] = 1

while(t<T):
    new = lattice.copy()
    t += 1
    for i in range (L-1):
        for j in range (L-1):

            a = count_neighbors(i,j)
            
            if a == 3 :
                new[i,j] = 1
            elif a<2 or a>=4 :
                new[i,j] = 0
    lattice = new

    clf()
    imshow(lattice)
    pause(0.1)
