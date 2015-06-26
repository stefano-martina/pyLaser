#!/usr/bin/python

#import numpy as np
#import matplotlib.pyplot as plt
import matplotlib.animation as ani
from pylab import *

#constants
k = 10    #decay rate in laser cavity (beam trasmission) (>0)
g1 = 1000   #decay rates of atomic polarization (>0)
g2 = 1000   #decay rates for population inversion (>0)

lMin = -2  #pumping energy parameter (in R)
lMax = 1
lSteps = 0.01

graphLimit = [[-5, 5],[-20, 20]]


maxwellAdiabaticEl = lambda k, g1, g2, l: lambda E: k*l*(E-E**3)/(l*E*E+1)
    

fig = plt.figure(figsize=(13, 7));
fig.canvas.set_window_title('Stability of Maxwell-Bloch')
plt.title('Stability study of Maxwell-Bloch equations with adiabatic elimination')
plt.xlabel('$E$')
plt.ylabel('$\dot{E}$')
plt.axhline(0, color='green')

line, = plt.plot([],[])

xlim(graphLimit[0])
ylim(graphLimit[1])

plt.figtext(0.7, 0.80, '$\dot{E} = \kappa\lambda(E-E^3)/(\lambda E^2+1)$')
lText = plt.figtext(0.7, 0.75, '')


pause = False

def onClick(event):
    global pause
    pause ^= True
    
fig.canvas.mpl_connect('button_press_event', onClick)

def init():
    line.set_data([], [])
    lText.set_text('')
    return line, lText

def makeGenerator(lMin, lMax, lSteps):
    def generator():
        l = lMin
        while l < lMax:
            if not pause:
                l = l + lSteps
            yield l
    return generator

def step(k, g1, g2, lMin, lSteps):
    def realStep(l):
        E = np.linspace(graphLimit[0][0],graphLimit[0][1], 1000)
        dE = maxwellAdiabaticEl(k, g1, g2, l)(E)
        line.set_data(E,dE)
        lText.set_text('$\lambda$ = %.2f' % l)
        return line, lText
    return realStep


anim = ani.FuncAnimation(fig, step(k, g1, g2, lMin, lSteps), frames=makeGenerator(lMin, lMax, lSteps), init_func=init, interval=20, blit=False, repeat=True)

plt.show()

