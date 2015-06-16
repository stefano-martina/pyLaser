import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#constants
G = 2      #gain coefficient
k = 100    #photon loss rate
f = 100     #decadiment rate for spontaneous emission
p = 5000     #pump force
#pMin = 1000
#pMax = 10000

def makeMilonni(G, k, f, p):
    def ode(state, t=0.):
        n = state[0]  #photons in the system
        N = state[1]  #excited atoms

        nd = G*n*N - k*n
        Nd = -G*n*N - f*N + p

        return [nd, Nd]    
    return ode


fig = plt.figure(figsize=(13, 10));
plt.subplots_adjust(wspace=0.3, hspace=0.3)

plt.subplot(221)
plt.title('Milonni model, evolution of $n$ from $80$')
plt.xlabel('$t$')
plt.ylabel('$n$')

state0 = [80,80]
t = np.linspace(0.0, 0.2, 300)
state = scipy.integrate.odeint(makeMilonni(G, k, f, p), state0, t)
plt.plot(t,state[:,0])

plt.subplot(222)
plt.title('Milonni model, evolution of $(n,N)$ from $(80,80)$')
plt.xlabel('$n$')
plt.ylabel('$N$')

plt.plot(state[:,0], state[:,1])

plt.subplot(223)
plt.title('Vector field')
plt.xlabel('$n$')
plt.ylabel('$N$')


status = np.meshgrid( np.linspace(0,100,20),np.linspace(0,100,20) )
endStatus = makeMilonni(G, k, f, p)(status)

quiver(status[0], status[1], endStatus[0], endStatus[1])

plt.show()

