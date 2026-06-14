from pylab import *

### functions
def initial():

    for i in range(N):
        check = 0
        while(check == 0):
            x[i] = d/2 +(L-d)*uniform(0,1)
            y[i] = d/2 +(L-d)*uniform(0,1)
            check = 1

            for j in range(i):
                dx = x[i] - x[j]
                dy = y[i] - y[j]
                dr = sqrt(dx*dx + dy*dy)
                if (dr < d):
                    check = 0

    for i in range(N):
        vx[i] = uniform(-1,1)
        vy[i] = uniform(-1,1)

###
def compute_force(px,py):

    fx[:] = 0
    fy[:] = 0

    for i in range(N):
        for j in range(i):

            dx = px[i] - px[j]
            dy = py[i] - py[j]

            r2 = dx*dx + dy*dy
            r = sqrt(r2)

            if r < 3*sigma and r > 1e-6:

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
def wall_reflect():

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

### runge-kuttah
def rk2_step():

    compute_force(x,y)

    ax1 = fx/m
    ay1 = fy/m

    x_mid = x + 0.5*dt*vx
    y_mid = y + 0.5*dt*vy

    vx_mid = vx + 0.5*dt*ax1
    vy_mid = vy + 0.5*dt*ay1

    compute_force(x_mid,y_mid)

    ax2 = fx/m
    ay2 = fy/m

    x[:] += dt*vx_mid
    y[:] += dt*vy_mid

    vx[:] += dt*ax2
    vy[:] += dt*ay2

###
def lj_potential(r):
    sr6 = (sigma/r)**6
    sr12 = sr6*sr6
    return 4*epsilon*(sr12 - sr6)

###
def plot_lj():

    r = linspace(0.8*sigma, 3*sigma, 500)
    V = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

    rmin = 2**(1/6)*sigma

    figure(1)
    plot(r,V)
    axvline(rmin, linestyle='--')
    axhline(0)
    xlabel("r")
    ylabel("V(r)")
    title("Lennard-Jones Potential")
    grid()

### variables
L = 1
N = 20
d = 0.1
sigma = 0.05
epsilon = 1.0
m = 1.0

dt = 0.001
T = 1

x = zeros(N)
y = zeros(N)
vx = zeros(N)
vy = zeros(N)

fx = zeros(N)
fy = zeros(N)

### main
initial()
plot_lj()

t = 0
while t < T:

    rk2_step()
    wall_reflect()

    t += dt

    figure(2)        # ← simulation figure
    clf()
    xlim(0,L)
    ylim(0,L)
    plot(x,y,'o',markersize=10)
    pause(0.01)

show()
