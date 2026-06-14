from numpy import *
from pylab import *

### functions
def initialize():

    for i in range(N):
        placed = False
        while not placed:
            x[i] = uniform(0,L)
            y[i] = uniform(0,L)

            placed = True
            for j in range(i):
                dx = x[i]-x[j]
                dy = y[i]-y[j]
                r = sqrt(dx*dx + dy*dy)
                if r < 1.2*sigma:
                    placed = False

    vx[:] = uniform(-1,1,N)
    vy[:] = uniform(-1,1,N)

###
def compute_force():

    fx[:] = 0
    fy[:] = 0

    for i in range(N):
        for j in range(i):

            dx = x[i]-x[j]
            dy = y[i]-y[j]
            r2 = dx*dx + dy*dy

            if r2 < rcut*rcut and r2 > 1e-12:

                r = sqrt(r2)

                sr6 = (sigma/r)**6
                sr12 = sr6*sr6

                F = 24*epsilon*(2*sr12 - sr6)/r2

                fxij = F*dx
                fyij = F*dy

                fx[i] += fxij
                fy[i] += fyij
                fx[j] -= fxij
                fy[j] -= fyij

###
def reflect_walls():
    for i in range(N):

        if x[i] < 0:
            x[i] = -x[i]
            vx[i] *= -1

        if x[i] > L:
            x[i] = 2*L - x[i]
            vx[i] *= -1

        if y[i] < 0:
            y[i] = -y[i]
            vy[i] *= -1

        if y[i] > L:
            y[i] = 2*L - y[i]
            vy[i] *= -1

###
def energy():
    KE = 0
    PE = 0
    KE = 0.5*m*sum(vx*vx + vy*vy)

    for i in range(N):
        for j in range(i):
            dx = x[i]-x[j]
            dy = y[i]-y[j]
            r = sqrt(dx*dx + dy*dy)

            if r < rcut:
                sr6 = (sigma/r)**6
                sr12 = sr6*sr6
                PE += 4*epsilon*(sr12 - sr6)

    return KE, PE

###
def verlet_step():

    vx[:] += 0.5*dt*fx/m
    vy[:] += 0.5*dt*fy/m

    x[:] += dt*vx
    y[:] += dt*vy

    reflect_walls()

    compute_force()

    vx[:] += 0.5*dt*fx/m
    vy[:] += 0.5*dt*fy/m

### variables
L = 1.0
N = 20
sigma = 0.05
epsilon = 1.0
m = 1.0

dt = 0.0005
T = 1.5
rcut = 3*sigma

x = zeros(N)
y = zeros(N)
vx = zeros(N)
vy = zeros(N)

fx = zeros(N)
fy = zeros(N)

### main
initialize()
compute_force()

t = 0
times = []
energies = []

figure(1)

while t < T:

    verlet_step()

    KE, PE = energy()
    times.append(t)
    energies.append(KE+PE)

    clf()
    xlim(0,L)
    ylim(0,L)
    plot(x,y,'o',markersize=10)
    pause(0.01)

    t += dt

figure(2)
plot(times, energies)
xlabel("Time")
ylabel("Total Energy")
title("Energy Conservation")
show()
