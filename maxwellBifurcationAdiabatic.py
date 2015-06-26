#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.animation as ani
#from pylab import *

#constants
k = 10.    #decay rate in laser cavity (beam trasmission) (>0)

lMin = -2.  #pumping energy parameter (in R)
lMax = 1.
lSteps = 0.01

graphLimit = [[-5., 5.],[-20., 20.]]
graph2YLimit = [-2., 2.]


funcE = lambda k, l: lambda E: k*l*(E-E**3)/(l*E*E+1)
funcEd = lambda k, l: lambda E: k*l*(-l*E**4 - (l+3)*E**2 + 1) / (l*E**2+1)**2
    

fig = plt.figure(figsize=(13, 8));
plt.subplots_adjust(hspace=0.3)

#up plot
ax1 = plt.subplot(211, xlim=(graphLimit[0][0], graphLimit[0][1]), ylim=(graphLimit[1][0],graphLimit[1][1]))
fig.canvas.set_window_title('Stability of Maxwell-Bloch')
plt.title('Stability study of Maxwell-Bloch equations with adiabatic elimination')
plt.xlabel('$E$')
plt.ylabel('$\dot{E}$')
plt.axhline(0, color='black')

line, = ax1.plot([],[], label='$\dot{E} = \kappa\lambda(E-E^3)/(\lambda E^2+1)$')
lined, = ax1.plot([],[], label='$d/dE [\ \kappa\lambda(E-E^3)/(\lambda E^2+1)\ ]$')

ax1.legend(handles=[line, lined])

stable, = ax1.plot([],[],'go')
unstable, = ax1.plot([],[],'ro')

asymptote, = ax1.plot([],[],'bo')

lText = plt.figtext(0.8, 0.74, '')

#down plot
ax2 = plt.subplot(212, ylim=(graph2YLimit[0], graph2YLimit[1]))
plt.axvline(0, color='black')

#plt.title('Stability study of Maxwell-Bloch equations with adiabatic elimination')
plt.xlabel('$\lambda$')
plt.ylabel('$E$')

l1 = np.linspace(lMin, -1, 1000, endpoint=False)
l2 = np.linspace(-1, 0, 1000, endpoint=False)
l3 = np.linspace(0, lMax, 1000)

fpA = lambda l: -1
fpB = lambda l: -1/math.sqrt(-l)
fpC = lambda l: 0
fpD = lambda l: 1/math.sqrt(-l)
fpE = lambda l: 1

line2, = ax2.plot([],[], 'black')

ax2.plot(l1, list(map(fpA, l1)), 'g')
ax2.plot(l2, list(map(fpA, l2)), 'r')
ax2.plot(l3, list(map(fpA, l3)), 'g')

ax2.plot(l1, list(map(fpB, l1)), 'r')
ax2.plot(l2, list(map(fpB, l2)), 'g')

ax2.plot(l1, list(map(fpC, l1)), 'g')
ax2.plot(l2, list(map(fpC, l2)), 'g')
ax2.plot(l3, list(map(fpC, l3)), 'r')

ax2.plot(l1, list(map(fpD, l1)), 'r')
ax2.plot(l2, list(map(fpD, l2)), 'g')

ax2.plot(l1, list(map(fpE, l1)), 'g')
ax2.plot(l2, list(map(fpE, l2)), 'r')
ax2.plot(l3, list(map(fpE, l3)), 'g')


pause = False

def onClick(event):
    global pause
    pause ^= True
    
fig.canvas.mpl_connect('button_press_event', onClick)

def init():
    line.set_data([], [])
    lined.set_data([], [])
    line2.set_data([], [])
    lText.set_text('')
    stable.set_data([], [])
    unstable.set_data([], [])
    asymptote.set_data([], [])
    return line, lined, line2, lText, stable, unstable, asymptote

def makeGenerator(lMin, lMax, lSteps):
    def generator():
        l = lMin
        while l < lMax:
            if not pause:
                l = l + lSteps
            yield l
    return generator

def step(k, lMin, lSteps):
    def realStep(l):
        E = np.linspace(graphLimit[0][0],graphLimit[0][1], 1000)
        fE = funcE(k, l)(E)
        line.set_data(E,fE)
        fEd = funcEd(k, l)(E)
        lined.set_data(E,fEd)
        lText.set_text('$\lambda$ = %.2f' % l)

        stableX = []
        stableY = []
        unstableX = []
        unstableY = []

        for root in [-1, 0, 1]:
            if funcEd(k, l)(root) < 0:
                stableX.append(root)
                stableY.append(0)
            else:
                unstableX.append(root)
                unstableY.append(0)

        stable.set_data(stableX, stableY)
        unstable.set_data(unstableX, unstableY)

        if l < 0:
            asymptote.set_data([-1/math.sqrt(-l), 1/math.sqrt(-l)], [0.,0.])
        
        line2.set_data([l,l],[graph2YLimit[0], graph2YLimit[1]])
        
        return line, lined, line2, lText, stable, unstable, asymptote
    return realStep


anim = ani.FuncAnimation(fig, step(k, lMin, lSteps), frames=makeGenerator(lMin, lMax, lSteps), init_func=init, interval=20, blit=False, repeat=True)

plt.show()

