import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from pylab import *

#constants
G = 2      #gain coefficient
k = 100    #photon loss rate
f = 100     #decadiment rate for spontaneous emission
p = 10000     #pump force


def makeMilonni(G, k, f, p):
    def ode(state, t=0.):
        n = state[0]  #photons in the system
        N = state[1]  #excited atoms

        nd = G*n*N - k*n
        Nd = -G*n*N - f*N + p

        return [nd, Nd]    
    return ode


fig = plt.figure(figsize=(13, 5));
plt.subplot(121)
plt.xlabel('$t$')
plt.ylabel('$n$')

state0 = [80,80]
t = np.linspace(0.0, 0.2, 100)
state = scipy.integrate.odeint(makeMilonni(G, k, f, p), state0, t)
plt.plot(t,state[:,0])

plt.subplot(122)
plt.xlabel('$n$')
plt.ylabel('$N$')

#n = np.linspace(0., 250., 100).tolist()
#nd = list(map(makeMilonni(G, k, f, p), n))

#plt.plot(n,nd)

plt.plot(state[:,0], state[:,1])

plt.figure(2)
plt.xlabel('$n$')
plt.ylabel('$N$')


status = np.meshgrid( np.linspace(0,100,20),np.linspace(0,100,20) )
endStatus = makeMilonni(G, k, f, p)(status)
#U = endStatus[:,0]
#V = endStatus[:,1]

quiver(status[0], status[1], endStatus[0], endStatus[1])

plt.show()





'''
fig = plt.figure(figsize=(13, 5));
plt.subplots_adjust(wspace=0.3)
ax1 = plt.subplot(121,xlim=(0, 1), ylim=(0, 5))
plt.xlabel('$t$ (sec)')
plt.ylabel('$n$ (photons)')
plt.title('Haken laser model')
line1, = ax1.plot([], [], lw=2)

ax2 = plt.subplot(122,xlim=(0, 6), ylim=(-100, 100))
plt.xlabel('$n$')
plt.ylabel('$\dot{n}$')
plt.title('Haken $n$ vs $\dot{n}$')
plt.axhline(0, color='green')
line2, = ax2.plot([], [], lw=2)
N0Text = ax2.text(0.70, 0.95, '', transform=ax2.transAxes)
GText = ax2.text(0.70, 0.90, '', transform=ax2.transAxes)
kText = ax2.text(0.70, 0.85, '', transform=ax2.transAxes)
aText = ax2.text(0.70, 0.80, '', transform=ax2.transAxes)

def init1():
    line1.set_data([], [])
    return line1,

def step1(G, k, a):
    def realStep1(N0):
        n0 = 1.
        t = np.linspace(0.0, 1.0, 100)
        n = scipy.integrate.odeint(makeHaken(N0, G, k, a), n0, t)

        line1.set_data(t, n)
        #print(N0)
        return line1,
    return realStep1

def init2():
    line2.set_data([], [])
    N0Text.set_text('')
    GText.set_text('')
    kText.set_text('')
    aText.set_text('')
    return line2, N0Text, GText, kText, aText

def step2(G, k, a):
    def realStep2(N0):
        n = np.linspace(0., 10., 100).tolist()
        nd = list(map(makeHaken(N0, G, k, a), n))

        line2.set_data(n, nd)

        N0Text.set_text('N0 = %d' % N0)
        GText.set_text('G = %.1f' % G)
        kText.set_text('k = %.1f' % k)
        aText.set_text('$\\alpha$ = %d' % a)
        return line2, N0Text, GText, kText, aText
    return realStep2

frameInt = 70
anim1 = ani.FuncAnimation(fig, step1(G, k, a), init_func=init1, frames=100, interval=frameInt, blit=True)
anim2 = ani.FuncAnimation(fig, step2(G, k, a), init_func=init2, frames=100, interval=frameInt, blit=True)
plt.show()
'''
    
