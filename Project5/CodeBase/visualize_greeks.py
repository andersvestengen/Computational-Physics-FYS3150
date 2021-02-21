import numpy as np
import matplotlib.pyplot as plt
from scipy import special, stats
from scipy.interpolate import UnivariateSpline
from scipy.misc import derivative

"""
This code is used to calculate analytical values for thr greeks and vizuallize
the numerical and analytical values for the greeks.
"""
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


def d1(S_t, ti,sigma):


    return 1 / (sigma * np.sqrt(ti)) * (np.log(S_t / E) + (r + sigma**2 / 2) * (ti))


def d2(S_t, ti,sigma):
    return d1(S_t, ti,sigma) - sigma * np.sqrt(ti)


def N(d):
    return stats.norm.cdf(d)  # (special.erf(d/np.sqrt(2)) + 1)/np.sqrt(2)

# GREEEEKS


def delta(V, S):
    V_spl = UnivariateSpline(S, V, s=0, k=4)
    V_spl_1d = V_spl.derivative(n=1)
    return V_spl_1d(S)


def gamma(V, S):
    dvds = np.diff(V, 1) / np.diff(S, 1)
    ddvdds = np.diff(dvds, 1) / np.diff(0.5 * (S[:-1] + S[1:]), 1)

    return ddvdds


sdata = open("greeks_s.txt", "r")

VS = []
sig = []
tsigma = []

for val in sdata.readline().split():
    sig.append(float(val))

for line in sdata:
    tsigma.append(float(line.split()[0]))
    rest_list = []

    for element in line.split()[1:]:
        rest_list.append(float(element))

    VS.append(rest_list)

VS = np.asarray(VS)
sig = np.asarray(sig)
tsigma = np.asarray(tsigma)


def vega(VS, sig):

    firstdiv = np.diff(VS) / np.diff(sig)
    return firstdiv

rdata = open("greeks_r.txt", "r")
VR = []
rho = []
trho = []

for val in rdata.readline().split():
    rho.append(float(val))

for line in rdata:
    trho.append(float(line.split()[0]))
    rest_list = []

    for element in line.split()[1:]:
        rest_list.append(float(element))

    VR.append(rest_list)

VR = np.asarray(VR)
rho = np.asarray(rho)
trho = np.asarray(trho)




def rhos(VR, rho):
    firstdiv = np.diff(VR) / np.diff(rho)
    return firstdiv


for i in range(1, len(t)):
    plt.plot(S[1:-1], delta(V[i, 1:-1], S[1:-1]),
             label=r"$\Delta$ for t = {:.1f}".format(T - t[i]))

plt.legend()
plt.xlabel("Price of underlying asset")
plt.ylabel(r"$\Delta$ ")
plt.title(r"Greek $\Delta$ as function of stock price")
plt.savefig("Results/delta.jpeg")
plt.show()



for i in range(1, len(t)):
    plt.plot(S[1:-3], gamma(V[i, 1:-1], S[1:-1]),
             label=r"$\gamma$ for t = {:.1f}".format(T - t[i]))
plt.legend()
plt.xlabel("Price of underlying asset")
plt.ylabel(r"$\gamma$ ")
plt.title(r"Greek $\gamma$ as function of stock price")
plt.savefig("Results/gamma.jpeg")
plt.show()

# Rho and Vega
for i in range(len(tsigma)):
    plt.plot(sig[:-1], vega(VS[i], sig),
             label=r"$\nu$($\sigma$), t={:.1f}".format(T - tsigma[i]))
plt.legend()
plt.xlabel("Volatility")
plt.ylabel(r"vega ")
plt.title(r"Greek vega as function of volatility")
plt.savefig("Results/vega.jpeg")
plt.show()

for i in range(len(trho)):
    plt.plot(rho[:-1], rhos(VR[i], rho),
             label=r"$\rho$(r), t={:.1f}".format(T - trho[i]))

plt.legend()
plt.xlabel("Risk free interest rate")
plt.ylabel(r"$\rho$ ")
plt.title(r"Greek $\rho$ as function of risk free interest rate")
plt.savefig("Results/rho.jpeg")
plt.show()


# FOOOR TAUUU

vtid = np.zeros(len(t))
for i in range(len(t)):
    vtid[i] = V[i, -1]


def theta(vtid, t):
    dvdt = -np.diff(vtid) / np.diff(t)
    return dvdt


plt.plot(t[:-1], theta(vtid, t), label=r"$\Theta$ as function of time")
plt.legend()
plt.xlabel(r"Time $\tau$")
plt.ylabel(r"$\tau$ ")
plt.title(r"Greek $\Theta$ as function of $\tau$")
plt.savefig("Results/tau.jpeg")
plt.show()

# Derivation of anaytical expression

St = S[-2]

def n(x):
    return np.exp(-x**2/2)/(np.sqrt(2*np.pi))

def altdelta(St,tau):
    return np.exp(-D*tau)*N(d1(St,tau,sigma))

def altgamma(S,tau):
    return np.exp(-D*tau)*n(d1(S,tau,sigma))/(np.sqrt(tau)*S*sigma)

def altvega(S,tau,sig):
    return E*np.exp(-r*tau)*n(d2(S,tau,sig))*np.sqrt(tau)

def alttheta(S,tau):
    return -np.exp(-D*tau)*S*n(d1(S,tau,sigma))*sigma/(np.sqrt(tau)*2) - \
    E*r*np.exp(-r*tau)*N(d2(S,tau,sigma)) + D*S*np.exp(-D*tau)*N(d1(S,tau,sigma))


def altrho(r,tau):
    return E*tau*np.exp(-r*tau)*N(d1(St,tau,sigma))

for i in range(len(trho)-1):
    plt.plot(S, altdelta(S, trho[i]),
             label=r"$\Delta$ for t = {:.1f}".format(T - trho[i]))
plt.legend()
plt.xlabel("Price of underlying asset")
plt.ylabel(r"$\Delta_{analytical}$ ")
plt.title(r"Greek $\Delta_{analytical}$ as function of stock price")
plt.savefig("Results/delta_ana.jpeg")
plt.show()


for i in range(len(trho)-1):
    plt.plot(S, altgamma(S, trho[i]),
             label=r"$\gamma$ for t = {:.1f}".format(T - trho[i]))
plt.legend()
plt.xlabel("Price of underlying asset")
plt.ylabel(r"$\gamma_{analytical}$ ")
plt.title(r"Greek $\gamma_{analytical}$ as function of stock price")
plt.savefig("Results/gamma_ana.jpeg")
plt.show()

plt.plot(t[1:], alttheta(St, t[1:]),
         label=r"$\Theta$($\tau$)")
plt.legend()
plt.xlabel(r"Time $\tau$ ")
plt.ylabel(r"$\Theta_{analytical}$ ")
plt.title(r"Greek $\Theta_{analytical}$ as function of time")
plt.savefig("Results/theta_ana.jpeg")
plt.show()

for i in range(len(tsigma)-1):
    plt.plot(sig, altvega(St, tsigma[i],sig),
             label=r"$\nu$ for t = {:.1f}".format(T - tsigma[i]))
plt.legend()
plt.xlabel(r"Volatility $\sigma$")
plt.ylabel(r"$\nu_{analytical}$ ")
plt.title(r"Greek $\nu_{analytical}$ as function of volatility")
plt.savefig("Results/vega_ana.jpeg")
plt.show()


for i in range(len(trho)-1):
    plt.plot(rho, altrho(rho,trho[i]),
             label=r"$\rho$ for t = {:.1f}".format(T - trho[i]))
plt.legend()
plt.xlabel("Risk free interest rate")
plt.ylabel(r"$\rho_{analytical}$ ")
plt.title(r"Greek $\rho_{analytical}$ as function of risk free interest rate")
plt.savefig("Results/rho_ana.jpeg")
plt.show()
