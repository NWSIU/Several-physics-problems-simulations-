from pylab import*
from matplotlib import*

### functions
def ax(x_m,y_m,x_M,y_M,M):
    r = sqrt((x_M-x_m)**2+(y_M-y_m)**2)
    return (M*(x_M-x_m))/(2*m*(r**3))
###
def ay(x_m,y_m,x_M,y_M,M):
    r = sqrt((x_M-x_m)**2+(y_M-y_m)**2)
    return M*(y_M-y_m)/(2*m*(r**3))
###
def aX(x_m,y_m,x_M,y_M,m):
    r = sqrt((x_M-x_m)**2+(y_M-y_m)**2)
    return -m*(x_M-x_m)/(2*M*(r**3))
###
def aY(x_m,y_m,x_M,y_M,m):
    r = sqrt((x_M-x_m)**2+(y_M-y_m)**2)
    return -m*(y_M-y_m)/(2*M*(r**3))
### variables

dt = 0.01
T = 30000

m = 400
M = 500

x = []
y = []
X = []
Y = []

### main
# initial condition
n = 0

x_m = 7
y_m = 7
b = y_m
vx_m = -0.1
vy_m = 0

x_M = 0
y_M = 0
vx_M = 0.2
vy_M = 0

# theta = uniform(0,2*pi)
# vx_m = -cos(theta)

plot(x_m , y_m ,'o')
plot(x_M , y_M ,'o')

while(n<T):

    t = n*dt

    x1_m = x_m + vx_m*dt/2
    y1_m = y_m + vy_m*dt/2
    x1_M = x_M + vx_M*dt/2
    y1_M = y_M + vy_M*dt/2
    vx1_m = vx_m + ax(x_m,y_m,x_M,y_M,M)*dt/2
    vy1_m = vy_m + ay(x_m,y_m,x_M,y_M,M)*dt/2
    vx1_M = vx_M + aX(x_m,y_m,x_M,y_M,m)*dt/2
    vy1_M = vy_M + aY(x_m,y_m,x_M,y_M,m)*dt/2

    x2_m = x_m + vx1_m*dt
    y2_m = y_m + vy1_m*dt
    x2_M = x_M + vx1_M*dt
    y2_M = y_M + vy1_M*dt
    vx2_m = vx_m + ax(x1_m,y1_m,x1_M,y1_M,M)*dt
    vy2_m = vy_m + ay(x1_m,y1_m,x1_M,y1_M,M)*dt
    vx2_M = vx_M + aX(x1_m,y1_m,x1_M,y1_M,m)*dt
    vy2_M = vy_M + aY(x1_m,y1_m,x1_M,y1_M,m)*dt


    x_m = x2_m
    y_m = y2_m
    vx_m = vx2_m
    vy_m = vy2_m
    x_M = x2_M
    y_M = y2_M
    vx_M = vx2_M
    vy_M = vy2_M


    x.append(x_m)
    y.append(y_m)
    X.append(x_M)
    Y.append(y_M)

    n += 1

plot(x_m , y_m ,'o')
plot(x_M , y_M ,'o')
plot(x,y,color='blue')
plot(X,Y,color='red')

xlabel("x")
ylabel("y")
title("Orbital Trajectory")
show()
