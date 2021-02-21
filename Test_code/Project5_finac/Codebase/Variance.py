import matplotlib.pyplot as plt
import numpy as np
V =[]
data = open('V_vis.txt','r')
for line in data:
    x1 = float(line.split()[0])
    V.append(x1)
mcs = np.linspace(1,len(V),len(V))
plt.plot(mcs,V)
plt.title("The standard deviation as a function of MCc")
plt.xlabel("Monte Carlo cycle [#]")
plt.ylabel("Average distance from m_0")
plt.show()
