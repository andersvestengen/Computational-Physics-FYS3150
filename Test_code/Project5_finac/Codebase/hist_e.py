import matplotlib.pyplot as plt
import numpy as np


def func(m):
    return np.max(np.asarray(n)) * np.exp(-1 / 100 * m)


names = ["Mon_vis.txt", "Mon_vis2.txt", "Mon_vis3.txt", "Mon_vis4.txt",
         "Mon_vis5.txt", "Mon_vis6.txt", "Mon_vis7.txt", "Mon_vis8.txt",
         "Mon_vis9.txt", "Mon_vis10.txt"]
alpha = [0, 1]
gamma = [0, 1, 2, 3, 4]
k = 0
j = 0
for i in range(len(names)):
    if i == 5:
        k += 1
        j = 0

    data = open(names[i], 'r')
    Mon = []
    for line in data:
        x1 = float(line.split()[0])
        Mon.append(x1)
    num_bins = len(Mon)
    Mon = np.asarray(Mon)
    plt.axvspan(0, 0.2 * np.max(Mon), facecolor='saddlebrown',
                alpha=0.5, label="Peasants (from Borgen)",)  # X coordinate
    plt.axvspan(np.max(Mon) * 0.2, 0.90 * np.max(Mon), facecolor='silver',
                alpha=0.5, label="Middle class (from Toten)")  # X coordinate
    plt.axvspan(np.max(Mon) * 0.90, np.max(Mon), facecolor='gold',
                alpha=0.5, label="Upper class (from Stavanger)")  # X coordinate
    n, bins, patches = plt.hist(Mon, 200, facecolor='b')
    m = np.linspace(0, np.max(Mon), 1000)

    plt.plot(m, func(m), label=r"$\omega_M$")
    # plt.axhspan(0, np.max(np.asarray(n)), facecolor='0.2', alpha=0.5) #Y coordinate
    plt.title(r"$\alpha$ = {} $\gamma$ = {}".format(alpha[k], gamma[j]))
    plt.xlabel("M")
    plt.ylabel("Amount of financal agents")
    plt.legend()
    plt.show()

    j += 1
