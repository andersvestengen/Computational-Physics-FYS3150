import matplotlib.pyplot as plt
import numpy as np
names = ["Savings025.txt", "Savings05.txt", "Savings09.txt"]
prosent = [25, 50, 90]
for i in range(len(names)):
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
    # plt.axhspan(0, np.max(np.asarray(n)), facecolor='0.2', alpha=0.5) #Y coordinate
    plt.title("Savings ={} %".format(prosent[i]))
    plt.xlabel("M")
    plt.ylabel("Amount of financal agents")
    plt.legend()
    plt.show()
