import numpy as np
import matplotlib.pyplot as plt
from scipy import special, stats
from scipy.interpolate import UnivariateSpline
from scipy.misc import derivative


sigma = 0.4
E = 50
r = 0.04
D = 0.12


data = open("u.txt", "r")

V = []
S = []
t = []

for val in data.readline().split():
    S.append(float(val))


for line in data:
    t.append(float(line.split()[0]))
    rest_list = []

    for element in line.split()[1:]:
        rest_list.append(float(element))

    V.append(rest_list)

V = np.asarray(V)
S = np.asarray(S)
t = np.asarray(t)
T = np.max(t)


for i in range(len(t)):
    plt.plot(S, V[i, :], label="V(S,t={:.1f})".format(T - t[i]))

plt.legend()
plt.xlabel("Price of underlying asset")
plt.ylabel("Value of option")
plt.title("Numerical solution to Black-Scholes")
plt.savefig("Results/Num_sol.jpeg")
plt.show()


def d1(S_t, ti,sigma):


    return 1 / (sigma * np.sqrt(ti)) * (np.log(S_t / E) + (r + sigma**2 / 2) * (ti))


def d2(S_t, ti,sigma):
    return d1(S_t, ti,sigma) - sigma * np.sqrt(ti)


def N(d):
    return stats.norm.cdf(d)  # (special.erf(d/np.sqrt(2)) + 1)/np.sqrt(2)


def Vana(S_t, ti, sigma, r):
    return N(d1(S_t, ti,sigma)) * S_t - N(d2(S_t, ti,sigma)) * E * np.exp(-r * (ti))


for i in range(1, len(t)):
    plt.plot(S[1:-1], Vana(S[1:-1], t[i], sigma, r),
             label="V(S,t={:.1f})".format(T - t[i]))

plt.legend()
plt.xlabel("Price of underlying asset")
plt.ylabel("Value of option")
plt.title("Analytical solution to Black-Scholes")
plt.savefig("Results/Av_sol.jpeg")
plt.show()


for i in range(1, len(t)):
    plt.plot(S, (np.abs(V[i, :] - Vana(S,t[i], sigma, r))),
             label="|V_dif|,t={:.1f}".format(T - t[i]))
plt.legend()
plt.xlabel("Price of underlying asset")
plt.ylabel("Value difference for option ")
plt.title("Value difference for analytical and numerical calculation of option")
plt.savefig("Results/Diff_ana_num.jpeg")
plt.show()
