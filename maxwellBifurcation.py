#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
#import sys
import scipy.integrate
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as ani
#from pylab import *

#constants
kMin = 1.
kMax = 1000.
kStep = 0.1
k = 5.    #decay rate in laser cavity (beam trasmission) (>0)

g1Min = 1.
g1Max = 1000.
g1Step = 10.
g1 = 1.   #decay rates of atomic polarization (>0)

g2Min = 1.
g2Max = 1000.
g2Step = 10.
g2 = 1.   #decay rates for population inversion (>0)

lMin = 11.  #pumping energy parameter (in R)
lMax = 100.
lStep = 0.01
lInt = 10.

l = lMin

graphLimit = [[0.45, 0.45, 0.45],[0.55, 0.55, 0.55]]
viewLimit = [[-3, -3, -7],[3, 3, 7]]
#graphLimit = [[-0.5, -0.5, -0.5],[0.5, 0.5, 0.5]]
#viewLimit = [[-2, -2, -2],[2, 2, 2]]
#graphLimit = [[-2, -2, -2],[2, 2, 2]]
#viewLimit = [[-10, -10, -10],[10, 10, 10]]

gridNum = 4

tMin = 0.1
tMax = 100.
tStep = 0.1
t = 100.  #integration time
integrationSteps = 10000

startPoints = []
for x in np.linspace(graphLimit[0][0], graphLimit[1][0], gridNum):
    for y in np.linspace(graphLimit[0][1], graphLimit[1][1], gridNum):
        for z in np.linspace(graphLimit[0][2], graphLimit[1][2], gridNum):
            startPoints.append([x,y,z])

startPoints = [[0.5,0.5,0.5]]


#E = S[0]
#P = S[1]
#D = S[2]

#full system
#Ed = k(P-E)
#Pd = g1(ED-P)
#Dd = g2(l+1-D-lEP)
maxwell = lambda k, g1, g2, l: lambda S, t:[k*(S[1]-S[0]),g1*(S[0]*S[2]-S[1]), g2*(l+1-S[2]-l*S[0]*S[1])]

#jacobian
#  -k       P       0
#  g1D    -g1P     g1E
# -g2lP   -g2lE     g2
maxwellJac = lambda k, g1, g2, l: lambda S, t:[[-k, S[1], 0], [g1*S[2], -g1*S[1], g1*S[0]], [-g2*l*S[1], -g2*l*S[0], g2]]

#adiabatic elimination system
#Ed = kE((l+1)/(lE^2+1) -1)
maxwellAdiabaticEl = lambda k, l: lambda E, t: k*l*(E-E**3)/(l*E*E+1)
adiabaticP = lambda l: lambda E: E*(l+1)/(l*E*E+1)
adiabaticD = lambda l: lambda E: (l+1)/(l*E*E+1)

fig = plt.figure(figsize=(13, 7));
ax = fig.gca(projection='3d')
fig.canvas.set_window_title('Maxwell-Bloch')
ax.set_title('Study of Maxwell-Bloch equations in E')
ax.set_xlabel('$E$')
ax.set_ylabel('$P$')
ax.set_zlabel('$D$')
ax.set_xlim(viewLimit[0][0], viewLimit[1][0])
ax.set_ylim(viewLimit[0][1], viewLimit[1][1])
ax.set_zlim(viewLimit[0][2], viewLimit[1][2])

line = []
lineA = []
for i in range(0, len(startPoints)):
    line[len(line):len(line)], = [plt.plot([],[],[])]
    lineA[len(lineA):len(lineA)], = [plt.plot([],[],[])]

plt.figtext(0.7, 0.80, '$\dot{E} = \kappa(P-E)$')
plt.figtext(0.7, 0.75, '$\dot{P} = \gamma_1(ED-P)$')
plt.figtext(0.7, 0.70, '$\dot{D} = \gamma_2(\lambda+1-D-\lambda EP)$')
kText = plt.figtext(0.7, 0.65, '')
g1Text = plt.figtext(0.7, 0.60, '')
g2Text = plt.figtext(0.7, 0.55, '')
lText = plt.figtext(0.7, 0.50, '')
tText = plt.figtext(0.7, 0.45, '')

pause = True
reverse = False
adiabatic = False

def onClick(event):
#    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata))

    global pause
    global reverse
    global adiabatic
    global k
    global g1
    global g2
    global l
    global t
    if event.key == ' ':
        pause ^= True
    
    if event.key == 'r':
        reverse ^= True

    if event.key == 'a':
        adiabatic ^= True

    if event.key == '9':
        if k < kMax:
            k = min(k + kStep, kMax)

    if event.key == '8':
        if k > kMin:
            k = max(k - kStep, kMin)

    if event.key == 'o':
        if g1 < g1Max:
            g1 = min(g1 + g1Step, g1Max)

    if event.key == 'i':
        if g1 > g1Min:
            g1 = max(g1 - g1Step, g1Min)

    if event.key == 'l':
        if g2 < g2Max:
            g2 = min(g2 + g2Step, g2Max)

    if event.key == 'k':
        if g2 > g2Min:
            g2 = max(g2 - g2Step, g2Min)

    if event.key == '.':
        if l < lMax:
            l = min(l + lStep, lMax)

    if event.key == ',':
        if l > lMin:
            l = max(l - lStep, lMin)

    if event.key == 'm':
        if t < tMax:
            t = min(t + tStep, tMax)

    if event.key == 'n':
        if t > tMin:
            t = max(t - tStep, tMin)

    if event.key == 'q':
        exit()

#fig.canvas.mpl_connect('button_press_event', onClick)
fig.canvas.mpl_connect('key_press_event', onClick)

def init():
    for i in range(0, len(startPoints)):
        line[i].set_data([], [])
        line[i].set_3d_properties([])
        lineA[i].set_data([], [])
        lineA[i].set_3d_properties([])
    kText.set_text('')
    g1Text.set_text('')
    g2Text.set_text('')
    lText.set_text('')
    tText.set_text('')
    return line, kText, g1Text, g2Text, lText

def makeGenerator(lMin, lMax, lStep):
    def generator():
        global l
        if not reverse:
            l = lMin
        else:
            l = lMax
            
        while l <= lMax+lStep and l >= lMin-lStep:
            if not pause:
                if not reverse:
                    l = l + lStep
                else:
                    l = l - lStep
            yield l
    return generator

def step(l):
    global adiabatic
    global k
    global g1
    global g2
    global t
    ts = np.linspace(0.0, t, integrationSteps)
    i = 0
    for sp in startPoints:
        state = scipy.integrate.odeint(maxwell(k, g1, g2, l), sp, ts, Dfun=maxwellJac(k, g1, g2, l))#, mxstep=1000)
        
        line[i].set_data(state[:,0],state[:,1])
        line[i].set_3d_properties(state[:,2])

        if adiabatic:
            state = scipy.integrate.odeint(maxwellAdiabaticEl(k, l), sp[0], ts)#, mxstep=1000)
            Es = state[:,0]

            Ps = list(map(adiabaticP(l), Es))
            Ds = list(map(adiabaticD(l), Es))
            lineA[i].set_data(Es,Ps)
            lineA[i].set_3d_properties(Ds)
        else:
            lineA[i].set_data([], [])
            lineA[i].set_3d_properties([])

        
        i = i + 1
            
    kText.set_text('$\kappa$ = %.2f' % k)
    g1Text.set_text('$\gamma_1$ = %.2f' % g1)
    g2Text.set_text('$\gamma_2$ = %.2f' % g2)
    lText.set_text('$\lambda$ = %.2f' % l)
    tText.set_text('$t$ = %.2f' % t)
    return line, kText, g1Text, g2Text, lText, tText


anim = ani.FuncAnimation(fig, step, frames=makeGenerator(lMin, lMax, lStep), init_func=init, blit=False, repeat=True) #, interval=lInt

plt.show()

