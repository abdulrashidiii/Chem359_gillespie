import numpy as np
import matplotlib.pyplot as plt

time, a, b, c = np.loadtxt('data', skiprows=1).T

plt.plot(time, a, label='A', c='black')
plt.plot(time, b, label='B', c='red')
plt.plot(time, c, label='C', c='blue')
#plt.plot(time, d, label='D', c='green')
#plt.plot(time, e, label='E', c='yellow')
plt.legend(loc='upper left')
plt.xlabel('Time')
plt.ylabel('Amounts')
plt.title('Gillespie Algorithm - Stochastic Simulation of Chemical Reactions')
plt.grid()
plt.show()
