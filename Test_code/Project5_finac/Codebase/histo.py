import matplotlib.pyplot as plt
import numpy as np
data = open('Mon_vis.txt', 'r')
Mon = []
for line in data:
    x1 = float(line.split()[0])
    Mon.append(x1)
num_bins = len(Mon)
Mon = np.asarray(Mon)
plt.axvspan(0, 0.2*np.max(Mon), facecolor='saddlebrown', alpha=0.5, label="Peasants (from Borgen)",) # X coordinate
plt.axvspan(np.max(Mon)*0.2, 0.95*np.max(Mon), facecolor='silver', alpha=0.5, label="Middle class (from Toten)") # X coordinate
plt.axvspan(np.max(Mon)*0.95 , np.max(Mon), facecolor='gold', alpha=0.5, label="Upper class (from Stavanger)") # X coordinate
n, bins, patches = plt.hist(Mon, 200, facecolor='b')
def func(m):
    return np.max(np.asarray(n))*np.exp(-1/100 * m)
#plt.axhspan(0, np.max(np.asarray(n)), facecolor='0.2', alpha=0.5) #Y coordinate
m = np.linspace(0,np.max(Mon),1000)
plt.title("True beauty of capitalism")
plt.xlabel("Money $$$")
plt.ylabel("Peepz :D")
plt.plot(m,func(m),color='k',label=r"$\omega_M$")
plt.legend()
plt.show()
