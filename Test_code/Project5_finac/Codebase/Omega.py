import matplotlib.pyplot as plt
import numpy as np
V = []
Mon = []
data = open('Omega.txt', 'r')
for line in data:
    x1 = float(line.split()[0])
    V.append(x1)
data = open('Mon_vis.txt', 'r')
for line in data:
    x1 = float(line.split()[0])
    Mon.append(x1)
mcs = np.linspace(1, len(V), len(V))
plt.plot(Mon, V)
plt.title("The standard deviation as a function of MCc")
plt.xlabel("M")
plt.ylabel(r"$\log(\omega_M)$")
plt.show()
