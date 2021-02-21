import matplotlib.pyplot as plt
import numpy as np

data = open('values.txt', 'r')
relerr_value = []
n_value = []


for line in data:
    relerr_value.append(float(line.split(' ')[1]))
    n_value.append(float(line.split(' ')[0]))




rel_value = np.asarray(relerr_value)
n_value = np.asarray(n_value)
hvals = np.log10(1/n_value)



plt.plot(hvals,relerr_value, label="$log10(\epsilon)$")

plt.xlabel('log10(h)')
plt.ylabel('$log10(\epsilon)$')
plt.legend()
plt.savefig('1_d.jpg')
plt.show()
