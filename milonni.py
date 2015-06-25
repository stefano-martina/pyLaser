#!/usr/bin/python3

import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from pylab import *

#constants
G = 2      #gain coefficient
k = 100    #photon loss rate
f = 100     #decadiment rate for spontaneous emission

pMin = 0  #pump force
pMax = 12000
pSteps = 100

startPoint = [80,80]
graphLimit = [90, 90]

tMax = 0.5
integrationSteps = 1000

def makeMilonni(G, k, f, p):
    def ode(state, t=0.):
        n = state[0]  #photons in the system
        N = state[1]  #excited atoms

        nd = G*n*N - k*n
        Nd = -G*n*N - f*N + p

        return [nd, Nd]    
    return ode


fig = plt.figure(figsize=(13, 7));
#fig = pyplot.gcf()
fig.canvas.set_window_title('Milonni model')
plt.subplots_adjust(wspace=0.3, hspace=0.3)

ax1 = plt.subplot(241, xlim=(0, tMax), ylim=(0, graphLimit[0]))
plt.title('Evolution of $n$ from $'+str(startPoint[0])+'$')
plt.xlabel('$t$')
plt.ylabel('$n$')

line1, = ax1.plot([],[])

ax1b = plt.subplot(242, xlim=(0, tMax), ylim=(0, graphLimit[1]))
plt.title('Evolution of $N$ from $'+str(startPoint[1])+'$')
plt.xlabel('$t$')
plt.ylabel('$N$')

line1b, = ax1b.plot([],[])

ax2 = plt.subplot(222, xlim=(0, graphLimit[0]), ylim=(0, graphLimit[1]))
plt.title('Evolution of $(n,N)$ from $('+str(startPoint[0])+','+str(startPoint[1])+')$')
plt.xlabel('$n$')
plt.ylabel('$N$')

cumulativeEndStatus = [[],[]];
line2, = ax2.plot([],[])
line2b, = ax2.plot([],[], color='red')

ax3 = plt.subplot(223, xlim=(0, graphLimit[0]), ylim=(0, graphLimit[1]))
#plt.title('Vector field')
plt.xlabel('$n$')
plt.ylabel('$N$')

status = np.meshgrid( np.linspace(0,graphLimit[0],20),np.linspace(0,graphLimit[1],20) )
endStatus = makeMilonni(G, k, f, pMin)(status)
q = ax3.quiver(status[0], status[1], endStatus[0], endStatus[1])

GText = plt.figtext(0.6, 0.40, '')
kText = plt.figtext(0.6, 0.35, '')
fText = plt.figtext(0.6, 0.30, '')
pText = plt.figtext(0.6, 0.25, '')

plt.figtext(0.68, 0.31, '{', fontsize=40) #size='xx-large')
plt.figtext(0.7, 0.35, '$\dot{n} = GnN - kn$')
plt.figtext(0.7, 0.30, '$\dot{N} = -GnN - fN + p$')

plt.draw()

pause = False

def onClick(event):
    global pause
    pause ^= True
    
fig.canvas.mpl_connect('button_press_event', onClick)

def init():
    line1.set_data([], [])
    line1b.set_data([], [])
    line2.set_data([], [])
    line2b.set_data([], [])

    GText.set_text('')
    kText.set_text('')
    fText.set_text('')
    pText.set_text('')

    return line1, line1b, line2, line2b, GText, kText, fText, pText

def step(G, k, f, pMin, pSteps):
    def realStep(i):#, q, status):
        if not pause:
            if i==0:
                cumulativeEndStatus[:] = [[],[]]
            t = np.linspace(0.0, tMax, integrationSteps)

            state = scipy.integrate.odeint(makeMilonni(G, k, f, pMin+(i*pSteps)), startPoint, t)
            line1.set_data(t,state[:,0])
            line1b.set_data(t,state[:,1])
            cumulativeEndStatus[0].append(state[:,0][-1])
            cumulativeEndStatus[1].append(state[:,1][-1])

            line2.set_data(state[:,0], state[:,1])
            line2b.set_data(cumulativeEndStatus[0], cumulativeEndStatus[1])

            endStatus = makeMilonni(G, k, f, pMin+(i*pSteps))(status)
            q.set_UVC(endStatus[0], endStatus[1])

            GText.set_text('$G$ = %.1f' % G)
            kText.set_text('$k$ = %.1f' % k)
            fText.set_text('$f$ = %.1f' % f)
            pText.set_text('$p$ = %d' % (pMin+(i*pSteps)))
        
        return line1, line1b, line2, line2b, cumulativeEndStatus, q, GText, kText, fText, pText
    return realStep


#biforcazioni

#two fix points stab and instab for n
stabPointsAn = [[],[]]
instabPointsAn = [[],[]]
stabPointsBn = [[],[]]
instabPointsBn = [[],[]]

#two fix points stab and instab for N
stabPointsAN = [[],[]]
instabPointsAN = [[],[]]
stabPointsBN = [[],[]]
instabPointsBN = [[],[]]


for p in np.linspace(pMin, pMax, 1000):
    if p > (k*f)/G:
        instabPointsAn[0].append(0);
        instabPointsAn[1].append(p);

        instabPointsAN[0].append(p);
        instabPointsAN[1].append(p/f);
    else:
        stabPointsAn[0].append(0);
        stabPointsAn[1].append(p);

        stabPointsAN[0].append(p);
        stabPointsAN[1].append(p/f);

    if p > 0 and p < (f*k)/G:
        instabPointsBn[0].append((p/k)-(f/G));
        instabPointsBn[1].append(p);

        instabPointsBN[0].append(p);
        instabPointsBN[1].append(k/G);
    else:
        stabPointsBn[0].append((p/k)-(f/G));
        stabPointsBn[1].append(p);

        stabPointsBN[0].append(p);
        stabPointsBN[1].append(k/G);

#fixpoints in n
fig2 = plt.figure(figsize=(13, 7));
fig2.canvas.set_window_title('Biforcazioni')
plt.subplots_adjust(wspace=0.3, hspace=0.3)

ax = plt.subplot(211)
plt.title('fixed points in $n$')
plt.xlabel('$n$')
plt.ylabel('$p$')

plot(stabPointsAn[0], stabPointsAn[1], color='green')
plot(stabPointsBn[0], stabPointsBn[1], color='green')
plot(instabPointsAn[0], instabPointsAn[1], color='red')
plot(instabPointsBn[0], instabPointsBn[1], color='red')

ax.annotate('A', xy=(stabPointsAn[0][0], stabPointsAn[1][0]), xytext=(stabPointsAn[0][0] + 3, stabPointsAn[1][0]), arrowprops=dict(facecolor='black', shrink=0.05))

ax.annotate('B', xy=(instabPointsBn[0][0], instabPointsBn[1][0]), xytext=(instabPointsBn[0][0] + 1, instabPointsBn[1][0] + 1500), arrowprops=dict(facecolor='black', shrink=0.05))

#fixpoint in N
ax = plt.subplot(212)
plt.title('fixed points in $N$')
plt.xlabel('$p$')
plt.ylabel('$N$')

plot(stabPointsAN[0], stabPointsAN[1], color='green')
plot(stabPointsBN[0], stabPointsBN[1], color='green')
plot(instabPointsAN[0], instabPointsAN[1], color='red')
plot(instabPointsBN[0], instabPointsBN[1], color='red')

ax.annotate('A', xy=(stabPointsAN[0][0], stabPointsAN[1][0]), xytext=(stabPointsAN[0][0], stabPointsAN[1][0] + 10), arrowprops=dict(facecolor='black', shrink=0.05))

ax.annotate('B', xy=(instabPointsBN[0][0], instabPointsBN[1][0]), xytext=(instabPointsBN[0][0], instabPointsBN[1][0] + 10), arrowprops=dict(facecolor='black', shrink=0.05))

plt.draw()

anim = ani.FuncAnimation(fig, step(G, k, f, pMin, pSteps), init_func=init, frames=math.ceil((pMax-pMin)/pSteps), interval=10, blit=False) #, fargs=(q, status)

plt.show()




