import scipy.integrate
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from pylab import *

#constants
k = 100    #decay rate in laser cavity (beam trasmission) (>0)
g1 = 1000   #decay rates of atomic polarization (>0)
g2 = 1000   #decay rates for population inversion (>0)
#l = 100    #pumping energy parameter (in R)

lMin = -3  #pumping energy parameter (in R)
lMax = 1000
lSteps = 100

#startPoint = [80,80]
graphLimitVar = 2
graphLimit = [[-graphLimitVar, -graphLimitVar, -graphLimitVar],[graphLimitVar, graphLimitVar, graphLimitVar]]
#graphLimit = [[-0.1, -2, -2],[0.1, 2, 2]]

quiverLen = 0.1
meshPoints = [10, 10, 10]

#tMax = 0.5
#integrationSteps = 1000

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
#ax3 = plt.subplot(111, xlim=(graphLimit[0][0], graphLimit[1][0]), ylim=(graphLimit[0][1], graphLimit[1][1]), zlim=(graphLimit[0][2],graphLimit[0][0]graphLimit[0][0]graphLimit[0][0]graphLimit[0][0] graphLimit[1][2]))
#plt.title('Vector field')
ax.set_xlabel('$E$')
ax.set_ylabel('$P$')
ax.set_zlabel('$D$')
plot([-1,0,1],[-1,0,1],[1,0,1], 'ro')
#plot([0],[0],[0], 'ro')

xMesh = np.linspace(graphLimit[0][0], graphLimit[1][0],meshPoints[0])
yMesh = np.linspace(graphLimit[0][1], graphLimit[1][1],meshPoints[1])
zMesh = np.linspace(graphLimit[0][2], graphLimit[1][2],meshPoints[2])
status = np.meshgrid( xMesh, yMesh, zMesh)
endStatus = makeMaxwell(k, g1, g2, lMin)(status)
q = ax.quiver(status[0], status[1], status[2], endStatus[0], endStatus[1], endStatus[2], length=quiverLen)



#diagr. biforc

stabPoints0 = [[],[]]
instabPoints0 = [[],[]]
stabPoints1 = [[],[]]
stabPointsm1 = [[],[]]
instabPoints1 = [[],[]]
instabPointsm1 = [[],[]]

fig2 = plt.figure(figsize=(13, 7));

for l in np.linspace(-10, 10, 1000):
    if l > -k/(k-1) :
        instabPoints0[0].append(l);
        instabPoints0[1].append(0);
    else:
        stabPoints0[0].append(l);
        stabPoints0[1].append(0);

    if k > (l*(l+1))/(1-l):
        instabPoints1[0].append(l);    
        instabPoints1[1].append(1);
     
        instabPointsm1[0].append(l);
        instabPointsm1[1].append(-1);    
    else:
        stabPoints1[0].append(l);    
        stabPoints1[1].append(1);
     
        stabPointsm1[0].append(l);
        stabPointsm1[1].append(-1);    

plot(stabPoints0[0], stabPoints0[1], color='green')
plot(stabPoints1[0], stabPoints1[1], color='green')
plot(stabPointsm1[0], stabPointsm1[1], color='green')
plot(instabPoints0[0], instabPoints0[1], color='red')
plot(instabPoints1[0], instabPoints1[1], color='red')
plot(instabPointsm1[0], instabPointsm1[1], color='red')

ylim([-1.5,1.5])

# fig.canvas.set_window_title('Milonni model')
# plt.subplots_adjust(wspace=0.3, hspace=0.3)

# ax1 = plt.subplot(241, xlim=(0, tMax), ylim=(0, graphLimit[0]))
# plt.title('Evolution of $n$ from $'+str(startPoint[0])+'$')
# plt.xlabel('$t$')
# plt.ylabel('$n$')

# line1, = ax1.plot([],[])

# ax1b = plt.subplot(242, xlim=(0, tMax), ylim=(0, graphLimit[1]))
# plt.title('Evolution of $N$ from $'+str(startPoint[1])+'$')
# plt.xlabel('$t$')
# plt.ylabel('$N$')

# line1b, = ax1b.plot([],[])

# ax2 = plt.subplot(222, xlim=(0, graphLimit[0]), ylim=(0, graphLimit[1]))
# plt.title('Evolution of $(n,N)$ from $('+str(startPoint[0])+','+str(startPoint[1])+')$')
# plt.xlabel('$n$')
# plt.ylabel('$N$')

# cumulativeEndStatus = [[],[]];
# line2, = ax2.plot([],[])
# line2b, = ax2.plot([],[], color='red')

# ax3 = plt.subplot(223, xlim=(-150, graphLimit[0]), ylim=(-150, graphLimit[1]))
# #plt.title('Vector field')
# plt.xlabel('$n$')
# plt.ylabel('$N$')

# status = np.meshgrid( np.linspace(-150,graphLimit[0],20),np.linspace(-150,graphLimit[1],20) )
# endStatus = makeMilonni(G, k, f, pMin)(status)
# q = ax3.quiver(status[0], status[1], endStatus[0], endStatus[1])

# GText = plt.figtext(0.6, 0.40, '')
# kText = plt.figtext(0.6, 0.35, '')
# fText = plt.figtext(0.6, 0.30, '')
# pText = plt.figtext(0.6, 0.25, '')

# plt.figtext(0.68, 0.31, '{', fontsize=40) #size='xx-large')
# plt.figtext(0.7, 0.35, '$\dot{n} = GnN - kn$')
# plt.figtext(0.7, 0.30, '$\dot{N} = -GnN - fN + p$')

pause = False

def onClick(event):
    global pause
    pause ^= True
    
fig.canvas.mpl_connect('button_press_event', onClick)

def init():
    return

    # line1.set_data([], [])
    # line1b.set_data([], [])
    # line2.set_data([], [])
    # line2b.set_data([], [])

    # GText.set_text('')
    # kText.set_text('')
    # fText.set_text('')
    # pText.set_text('')

    # return line1, line1b, line2, line2b, GText, kText, fText, pText

def step(k, g1, g2, lMin, lSteps):
    def realStep(i):
        if not pause:
            endStatus = makeMaxwell(k, g1, g2, lMin+(i*lSteps))(status)
            q.set_UVC(endStatus[0], endStatus[1], endStatus[2])
            return q

            
        #     if i==0:
        #         cumulativeEndStatus[:] = [[],[]]
        #     t = np.linspace(0.0, tMax, integrationSteps)

        #     state = scipy.integrate.odeint(makeMilonni(G, k, f, pMin+(i*pSteps)), startPoint, t)
        #     line1.set_data(t,state[:,0])
        #     line1b.set_data(t,state[:,1])
        #     cumulativeEndStatus[0].append(state[:,0][-1])
        #     cumulativeEndStatus[1].append(state[:,1][-1])

        #     line2.set_data(state[:,0], state[:,1])
        #     line2b.set_data(cumulativeEndStatus[0], cumulativeEndStatus[1])

        #     endStatus = makeMilonni(G, k, f, pMin+(i*pSteps))(status)
        #     q.set_UVC(endStatus[0], endStatus[1])

        #     GText.set_text('$G$ = %.1f' % G)
        #     kText.set_text('$k$ = %.1f' % k)
        #     fText.set_text('$f$ = %.1f' % f)
        #     pText.set_text('$p$ = %d' % (pMin+(i*pSteps)))
        
        # return line1, line1b, line2, line2b, cumulativeEndStatus, q, GText, kText, fText, pText
    return realStep

#anim = ani.FuncAnimation(fig, step(k, g1, g2, lMin, lSteps), init_func=init, frames=math.ceil((lMax-lMin)/lSteps), interval=10, blit=False) #, fargs=(q, status)
plt.show()

