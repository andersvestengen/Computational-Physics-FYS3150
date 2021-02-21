import matplotlib.pyplot as plt
import numpy as np

data = open('valn1.txt', 'r')
y_value = [0]
x_value = [0]

data1 = open('valn2.txt', 'r')
y1_value = [0]
x1_value = [0]

data2 = open('valn3.txt', 'r')
y2_value = [0]
x2_value = [0]
ex_vals = [0]

for line in data:
    l1, l2, l3 = line.split(' ')
    y_value.append(float(l1))
    x_value.append(float(l2))

for line in data1:
    l1, l2, l3 = line.split(' ')
    y1_value.append(float(l1))
    x1_value.append(float(l2))

for line in data2:
    l1, l2, l3 = line.split(' ')
    y2_value.append(float(l1))
    x2_value.append(float(l2))
    ex_vals.append(float(l3))

# Hardkoder inn 0 og 1 for x
x_value.append(1)
y_value.append(0)
y1_value.append(0)
x1_value.append(1)
y2_value.append(0)
x2_value.append(1)
ex_vals.append(0)

plt.plot(x2_value, ex_vals, "k-",label="Exact solution")
plt.plot(x_value,y_value,"r--", label="n=10")
plt.plot(x1_value, y1_value,"b--", label="n=100" )
plt.plot(x2_value, y2_value,"y--",label="n=1000" )
plt.xlabel('x=(0,1)')
plt.ylabel('u(x)')
plt.legend()
plt.savefig("1_c.jpeg")
plt.show()

data.close()
data1.close()
data2.close()
