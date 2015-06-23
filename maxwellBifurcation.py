#import numpy as np
#import matplotlib.pyplot as plt
#import sys
import scipy.integrate
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as ani
from pylab import *

#constants
k = 10    #decay rate in laser cavity (beam trasmission) (>0)
g1 = 10   #decay rates of atomic polarization (>0)
g2 = 10   #decay rates for population inversion (>0)

lMin = -2  #pumping energy parameter (in R)
lMax = 1
lStep = 0.01
lInt = 10

graphLimit = [[-0.5, -0.5, -0.5],[0.5, 0.5, 0.5]]
viewLimit = [[-1, -1, -1],[1, 1, 2]]
gridNum = 4

tMax = 0.5
integrationSteps = 1000

startPoints = []
for x in np.linspace(graphLimit[0][0], graphLimit[1][0], gridNum):
    for y in np.linspace(graphLimit[0][1], graphLimit[1][1], gridNum):
        for z in np.linspace(graphLimit[0][2], graphLimit[1][2], gridNum):
            startPoints.append([x,y,z])

#Ed = k(P-E)
#Pd = g1(ED-P)
#Dd = g2(l+1-D-lEP)

#E = S[0]
#P = S[1]
#D = S[2]
maxwell = lambda k, g1, g2, l: lambda S, t:[k*(S[1]-S[0]),g1*(S[0]*S[2]-S[1]), g2*(l+1-S[2]-l*S[0]*S[1])]

fig = plt.figure(figsize=(13, 7));
ax = fig.gca(projection='3d')
fig.canvas.set_window_title('Stability of Maxwell-Bloch')
ax.set_title('Stability study of Maxwell-Bloch equations in E')
ax.set_xlabel('$E$')
ax.set_ylabel('$P$')
ax.set_zlabel('$D$')
ax.set_xlim(viewLimit[0][0], viewLimit[1][0])
ax.set_ylim(viewLimit[0][1], viewLimit[1][1])
ax.set_zlim(viewLimit[0][2], viewLimit[1][2])

line = []
for i in range(0, len(startPoints)):
    line[len(line):len(line)], = [plt.plot([],[],[])]

plt.figtext(0.7, 0.80, '$\dot{E} = \kappa(P-E)$')
plt.figtext(0.7, 0.75, '$\dot{P} = \gamma_1(ED-P)$')
plt.figtext(0.7, 0.70, '$\dot{D} = \gamma_2(\lambda+1-D-\lambda EP)$')
lText = plt.figtext(0.7, 0.65, '')

pause = False
reverse = False
l = lMin

def onClick(event):
#    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata))

    global pause
    global reverse
    global l
    if event.key == ' ':
        pause ^= True
    
    if event.key == 'r':
        reverse ^= True

    if event.key == '.':
        if l < lMax:
            l = l + lStep

    if event.key == ',':
        if l > lMin:
            l = l - lStep

    if event.key == 'q':
        exit()

#fig.canvas.mpl_connect('button_press_event', onClick)
fig.canvas.mpl_connect('key_press_event', onClick)

def init():
    for i in range(0, len(startPoints)):
        line[i].set_data([], [])
        line[i].set_3d_properties([])
    lText.set_text('')
    return line, lText

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

def step(k, g1, g2, lMin, lStep):
    def realStep(l):
        t = np.linspace(0.0, tMax, integrationSteps)
        i = 0
        for sp in startPoints:
            state = scipy.integrate.odeint(maxwell(k, g1, g2, l), sp, t)
            
            line[i].set_data(state[:,0],state[:,1])
            line[i].set_3d_properties(state[:,2])
            i = i + 1
            
        lText.set_text('$\lambda$ = %.2f' % l)
        return line, lText
    return realStep


anim = ani.FuncAnimation(fig, step(k, g1, g2, lMin, lStep), frames=makeGenerator(lMin, lMax, lStep), init_func=init, interval=lInt, blit=False, repeat=True)

plt.show()

