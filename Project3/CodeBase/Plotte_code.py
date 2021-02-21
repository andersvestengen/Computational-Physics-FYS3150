import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

data = open('values_Euler.txt', 'r')
xE_value = []
yE_value = []
for line in data:
    l1, l2, l3 = line.split(' ')
    yE_value.append(float(l2))
    xE_value.append(float(l1))
ax = plt.axes(projection='3d')
zE_value = np.zeros(len(xE_value))
ax.plot3D(xE_value, yE_value, zE_value)

plt.show()
data.close()
