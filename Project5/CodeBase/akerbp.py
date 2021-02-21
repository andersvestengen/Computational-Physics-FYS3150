import numpy as np
import matplotlib.pyplot as plt
from scipy import special, stats
from scipy.interpolate import UnivariateSpline

data = open("Akerbp.txt", "r")

V1 = []
S1 = []

for line in data:
    S1.append(float(line.split()[0]))
    V1.append(float(line.split()[1]))


data = open("Aker_numeric.txt", "r")

V2 = []
S2 = []
t = []

for val in data.readline().split():
    S2.append(float(val))

for line in data:
    t.append(float(line.split()[0]))
    rest_list = []

    for element in line.split()[1:]:
        rest_list.append(float(element))

    V2.append(rest_list)
plt.axvspan(S1[-1]-0.1, S1[-1]+0.1,edgecolor='k',facecolor='None',alpha=1, label="Date:15.12, Stock price ={}".format(S1[-1]))
plt.axvspan(S1[-2]-0.1, S1[-2]+0.1,edgecolor='g',facecolor='None',alpha=1, label="Date:14.12, Stock price ={}".format(S1[-2]))
plt.axvspan(S1[-3]-0.1, S1[-3]+0.1,edgecolor='r',facecolor='None',alpha=1, label="Date:11.12, Stock price ={}".format(S1[-3]))
plt.axvspan(S1[-4]-0.1, S1[-4]+0.1,edgecolor='b',facecolor='None',alpha=1, label="Date:10.12, Stock price ={}".format(S1[-4]))
plt.legend()
S1 = np.asarray(S1)
ac_tim = [len(t)-1,len(t)-2,len(t)-3,len(t)-4]
plt.scatter(S1,np.asarray(V1),label='Real options values')
V_num=[]
S_num = []
for j in range(4):
    indx = int(np.where(np.abs(S1[j]-S2)==min(np.abs(S1[j]-S2)))[0][0])
    S_num.append(S2[indx])
    V_num.append(V2[int(ac_tim[j])][indx])
plt.scatter(S_num,V_num,label='Numerical option values')
plt.legend()
plt.xlabel('Stock-prizes (S)' ,fontsize=12)
plt.ylabel('Option-values (V)' ,fontsize=12)
plt.title('Numerical vs Analytical option-values for Aker BP',fontsize=12)
plt.show()
