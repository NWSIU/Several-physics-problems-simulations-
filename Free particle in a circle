from pylab import *
from numpy import *
from matplotlib import *

# description of the problem:
# free particle in a circular box


### variables
dt = 0.1
T = 10000
R = 5
v = 1
n = 0

### main
x = uniform(-R/1.5 , R/1.5)
y = uniform(-R/1.5 , R/1.5)

theta = uniform(0 , 2*pi)
vx = v*cos(theta)
vy = v*sin(theta)


# Circle boundary (got help from chat GPT)
theta_circ = linspace(0, 2*np.pi, 100)
a = R*cos(theta_circ)
b = R*sin(theta_circ)
plot(a, b, 'k--')

X = []
Y = []

while(n<T):
    x = x + vx*dt
    y = y + vy*dt

    X.append(x)
    Y.append(y)
    
    if x**2+y**2 > R**2:

        r = sqrt(x**2 + y**2)
        r_hat_x = x / r
        r_hat_y = y / r

        v_dot_r = vx * r_hat_x + vy*r_hat_y
        
        vx = vx - 2*v_dot_r*r_hat_x
        vy = vy - 2*v_dot_r*r_hat_y
        
    n += 1

plot(X,Y)
plot('equal')
show()
