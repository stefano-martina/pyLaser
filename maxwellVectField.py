#!/usr/bin/python

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

#constants
k = 5    #decay rate in laser cavity (beam trasmission) (>0)
g1 = 1   #decay rates of atomic polarization (>0)
g2 = 1   #decay rates for population inversion (>0)
l = 11    #pumping energy parameter (in R)

graphLimitVar = 3
graphLimit = [[-graphLimitVar, -graphLimitVar, -graphLimitVar],[graphLimitVar, graphLimitVar, graphLimitVar]]
#graphLimit = [[-0.1, -2, -2],[0.1, 2, 2]]

quiverLen = 0.3
meshPoints = [10, 10, 10]

def makeMaxwell(k, g1, g2, l):
    def ode(state, t=0.):
        E = state[0]  # electric field
        P = state[1]  # mean polarization
        D = state[2]  # population inversion

        nE = k*(P - E)
        nP = g1*(E*D - P)
        nD = g2*(l+1 - D - l*E*P)
        
        return [nE, nP, nD]    
    return ode


fig = plt.figure(figsize=(13, 7))
ax = fig.gca(projection='3d')
fig.canvas.set_window_title('Maxwell-Bloch')
ax.set_title('Vector field induced by Maxwell-Bloch equations')
ax.set_xlabel('$E$')
ax.set_ylabel('$P$')
ax.set_zlabel('$D$')
plt.plot([-1,0,1],[-1,0,1],[1,0,1], 'ro')
#plt.plot([0],[0],[0], 'ro')

xMesh = np.linspace(graphLimit[0][0], graphLimit[1][0],meshPoints[0])
yMesh = np.linspace(graphLimit[0][1], graphLimit[1][1],meshPoints[1])
zMesh = np.linspace(graphLimit[0][2], graphLimit[1][2],meshPoints[2])
status = np.meshgrid( xMesh, yMesh, zMesh)
endStatus = makeMaxwell(k, g1, g2, l)(status)
q = ax.quiver(status[0], status[1], status[2], endStatus[0], endStatus[1], endStatus[2], length=quiverLen)

plt.show()

