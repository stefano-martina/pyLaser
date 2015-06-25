#!/usr/bin/python3

import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from pylab import *

#constants
G = 2.      #gain coefficient
k = 100.    #photon loss rate
a = 10.     #atoms drop to ground rate


def makeHaken(N0, G, k, a):
    return lambda n,t=0.: (G * N0 - k) * n - (a * G) * n * n

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

        N0Text.set_text('$N_0$ = %d' % N0)
        GText.set_text('$G$ = %.1f' % G)
        kText.set_text('$k$ = %.1f' % k)
        aText.set_text('$\\alpha$ = %d' % a)
        return line2, N0Text, GText, kText, aText
    return realStep2

frameInt = 70
anim1 = ani.FuncAnimation(fig, step1(G, k, a), init_func=init1, frames=100, interval=frameInt, blit=True)
anim2 = ani.FuncAnimation(fig, step2(G, k, a), init_func=init2, frames=100, interval=frameInt, blit=True)


fig2 = plt.figure(2)
plt.xlabel('$N_0$')
plt.ylabel('$n$')
plt.title('Haken bifurcation')

N0,n = np.meshgrid( np.linspace(0,100,20),np.linspace(0,7,20) )
U = 0
V = makeHaken(N0, G, k, a)(n)

quiver( N0,n,U, V)

n0 = 1
N0s = []
ns = []
t =  np.linspace(0.0, 1.0, 100)
for N0 in np.linspace(0,100,100).tolist():
    n = scipy.integrate.odeint(makeHaken(N0, G, k, a), n0, t)
    N0s.append(N0)
    ns.extend(n[-1])

plt.plot(N0s, ns)

plt.show()
    
