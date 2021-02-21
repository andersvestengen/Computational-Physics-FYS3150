import matplotlib.pyplot as plt
import numpy as np
import collections
names = ["alpha_05.txt", "alpha_1.txt", "alpha_15.txt", "alpha_2.txt"]
lbl = [0.5, 1, 1.5, 2]
prosent = [25, 50, 90]


def func(m):
    return beta * np.exp(-beta * m)


for i in range(len(names)):
    data = open(names[i], 'r')
    Mon = []
    for line in data:
        x1 = float(line.split()[0])
        Mon.append(x1)
    Mon = np.asarray(Mon)
    bins = np.linspace(np.min(Mon), np.max(Mon), 200)
    binplace = np.digitize(Mon, bins)
    count = np.bincount(binplace)[:-1]
    beta = 1
    plt.plot(bins[1:], count[1:] / (len(Mon)), '--',
             label=r"$\alpha$={}".format(lbl[i]))
    # plt.plot(Mon[indx],func(Mon[indx]),label=r"$\alpha$={}".format(lbl[i]))

plt.xlabel("m")
plt.ylabel("P(m)")
plt.legend()
plt.xlim([0.001, 1000])
plt.ylim([1e-3, 1])
plt.yscale('log')
plt.xscale('log')
plt.show()
