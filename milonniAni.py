import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from pylab import *

#constants
G = 2      #gain coefficient
k = 100    #photon loss rate
f = 100     #decadiment rate for spontaneous emission
#p = 5000     #pump force
pMin = 1000
pMax = 10000
pSteps = 100

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

ax1 = plt.subplot(221, xlim=(0, 0.2), ylim=(0, 90))
plt.title('Milonni model, evolution of $n$ from $80$')
plt.xlabel('$t$')
plt.ylabel('$n$')

line1, = ax1.plot([],[])

ax2 = plt.subplot(222, xlim=(0, 90), ylim=(0, 90))
plt.title('Milonni model, evolution of $(n,N)$ from $(80,80)$')
plt.xlabel('$n$')
plt.ylabel('$N$')

cumulativeEndStatus = [[],[]];
line2, = ax2.plot([],[])
line2b, = ax2.plot([],[], color='red')

ax3 = plt.subplot(223, xlim=(0, 100), ylim=(0, 100))
plt.title('Vector field')
plt.xlabel('$n$')
plt.ylabel('$N$')

status = np.meshgrid( np.linspace(0,100,20),np.linspace(0,100,20) )
endStatus = makeMilonni(G, k, f, pMin)(status)
q = ax3.quiver(status[0], status[1], endStatus[0], endStatus[1])

ax4 = plt.subplot(224)

GText = ax4.text(0.2, 0.95, '', transform=ax4.transAxes)
kText = ax4.text(0.2, 0.90, '', transform=ax4.transAxes)
fText = ax4.text(0.2, 0.85, '', transform=ax4.transAxes)
pText = ax4.text(0.2, 0.80, '', transform=ax4.transAxes)


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line2b.set_data([], [])
    GText.set_text('')
    kText.set_text('')
    fText.set_text('')
    pText.set_text('')
    return line1, line2, line2b, GText, kText, fText, pText

def step(G, k, f, pMin, pSteps):
    def realStep(i):#, q, status):
        if i==0:
            cumulativeEndStatus[:] = [[],[]]
    
        state0 = [80,80]
        t = np.linspace(0.0, 0.2, 300)
        state = scipy.integrate.odeint(makeMilonni(G, k, f, pMin+(i*pSteps)), state0, t)
        line1.set_data(t,state[:,0])
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

        return line1, line2, line2b, cumulativeEndStatus, q, GText, kText, fText, pText
    return realStep

anim = ani.FuncAnimation(fig, step(G, k, f, pMin, pSteps), init_func=init, frames=math.ceil((pMax-pMin)/pSteps), interval=10, blit=False) #, fargs=(q, status)
plt.show()

