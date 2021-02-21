import matplotlib.pyplot as plt
import numpy as np
data = open('E.txt', 'r')
counter=[]
Energy= []

for line in data:
    x1 = float(line.split()[0])
    Energy.append(x1)
num_bins = len(Energy)
n, bins, patches = plt.hist(Energy, 300, density=True, facecolor='b')
plt.title("P(E) after the steady state situation is reached.")
plt.xlabel("E/J")
plt.ylabel("P(E)")
plt.show()
