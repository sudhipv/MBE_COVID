

"""## SEIRD Model

*Based on the form presented by Jha, P.K., Cao, L., and Oden, T.
"Bayesian-based predictions of COVID-19 evolution in Texas usingmultispecies mixture-theoretic continuum models"
Computational Mechanics (2020) 66:1055–1068 https://doi.org/10.1007/s00466-020-01889-z*
"
"""
# CREDITS : Brandon Robinson, Carleton University.

# from google.colab import drive
# drive.mount('/content/drive')

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

dt = 1/24
t = np.arange(0, 90, dt)

# Preallocate x for 5 realizations
x = np.zeros((len(t),5))

# i0 = 0.1     # Initial proportion of infected population
total = 1

i0 = 0.1*total
R = 0  ## Transmition Rate
e = R * i0
r = 0
d = 0
s = total - e - i0 - r - d


print('Susceptible', s)
print('exposed', e)
print('Infected', i0)
print('Recovered', r)
print('Deceased', d)

x[0,:] = np.array([s,e,i0,r,d])
# x[0,:] = np.array([1-i0,0,i0,0,0])

# Model parameters from Oden
gamma_e = 1/6
gamma_r = 1/24
gamma_d = 1/160
gamma_i = 1/7
beta_e = 3.78E-4*1000
beta_i = 3.78E-4*1000
A = 0.5


for kk in range(1,len(t)):

  x[kk,0] = x[kk-1,0] + dt*(- (1-A)*beta_e*x[kk-1,0]*x[kk-1,1] - (1-A)*beta_i*x[kk-1,0]*x[kk-1,2])
  x[kk,1] = x[kk-1,1] + dt*((1-A)* beta_e*x[kk-1,0]*x[kk-1,1] + (1-A)*beta_i*x[kk-1,0]*x[kk-1,2] - (gamma_i + gamma_e)*x[kk-1,1])
  x[kk,2] = x[kk-1,2] + dt*(gamma_i*x[kk-1,1] - (gamma_r + gamma_d)*x[kk-1,2])
  x[kk,3] = x[kk-1,3] + dt*(gamma_e*x[kk-1,1] + gamma_r*x[kk-1,2])
  x[kk,4] = x[kk-1,4] + dt*(gamma_d*x[kk-1,2])


tt = np.zeros((len(t),0))
print(np.shape(tt))
tt = x[:,0] +  x[:,1] + x[:,2] + x[:,3] + x[:,4]


# plt.plot(t, x[:,0], label='Susceptible')
plt.plot(t, x[:,1], label='Exposed')
plt.plot(t, x[:,2], label='Infectious')
plt.plot(t, x[:,3], label='Recoverd')
plt.plot(t, x[:,4], label='Dead')
# plt.plot(t, tt, label='Total')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('Proportion of the population')
plt.grid()
plt.show()

np.save('x.npy',x)
