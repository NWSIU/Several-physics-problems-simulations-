import random 
import matplotlib.pyplot as plt
import math
### variables
n = 10
k = 3

### main
U = 0
while U == 0:

    S = []
    for i in range(1,n+1):
        for j in range(k):
            S.append(i)
    
    random.shuffle(S)
    edges = []

    U = 1
    for i in range(0, len(S), 2):
        
        u = S[i]
        v = S[i + 1]
        (u,v) = (min(u,v),max(u,v))
        
        if u == v:
            U = 0
            break

        if (u,v) in edges:
            U = 0
            break

        edges.append((u,v))

print(edges)

### visualizing (got help from chat gbt) ###
# node positions on a circle
positions = {}

for i in range(1, n + 1):
    angle = 2 * math.pi * (i - 1) / n
    x = math.cos(angle)
    y = math.sin(angle)
    positions[i] = (x, y)

# draw edges
plt.figure()

for (u, v) in edges:
    x_values = [positions[u][0], positions[v][0]]
    y_values = [positions[u][1], positions[v][1]]
    plt.plot(x_values, y_values)

# draw nodes
for node in positions:
    x, y = positions[node]
    plt.scatter(x, y)
    plt.text(x, y, str(node))

plt.axis("equal")
plt.axis("off")
plt.show()
