from pylab import *
from matplotlib import *

### variables
t = 100
p = 0.5
# X = []
M = 50000


### main
avgtime = 0
avgtimeabs = 0
for s in range(M):
###
    x = 0
    x_pos = 0 
    x_neg = 0
    for n in range(t):
        if (uniform(0,1)<p):
            x+=1
        else:
            x-=1
        # X.append(x)

        if (x>0):
            x_pos += 1
        if (x<0):
            x_neg += 1
    
    avgtime += x_pos - x_neg   
    avgtimeabs += abs(x_pos - x_neg)

avgtime = avgtime/M  
avgtimeabs = avgtimeabs/M
# plot(X)
# show()
print(avgtime,avgtimeabs)
