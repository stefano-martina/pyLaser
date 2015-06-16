#from scipy.integrate import odeint
#from numpy import arange
#from pylab import plot,xlabel,ylabel,title,legend,show
import scipy.integrate
import numpy as np
#import pylab as plt
import matplotlib.pyplot as plt
import matplotlib.animation as ani

#constants
G = 2.      #gain coefficient
k = 100.    #photon loss rate
a = 10.     #atoms drop to ground rate


def makeHaken(N0, G, k, a):
    return lambda n,t=0.: (G * N0 - k) * n - (a * G) * n * n

'''
def haken(n, t0=0.):
    #constants
    G = 2.      #gain coefficient
    N0 = 100.  #initial exited atoms
    k = 100.    #photon loss rate
    a = 10.     #atoms drop to ground rate

    #calculate
    nd = (G * N0 - k) * n - (a * G) * n * n

    return nd
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

#n = np.linspace(0., 10., 100).tolist()
#zero = list(map(lambda x:0., n))
#zero = [0] * 100
#ax2.plot(n,zero)


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

#anim1.save('test1-1.mp4', fps=30)#, extra_args=['-vcodec', 'libx264'])
#anim2.save('test1-2.mp4', fps=30)

'''
plt.figure(1, figsize=(13, 5))
plt.subplots_adjust(wspace=0.3)
plt.subplot(121)
plt.plot(t, n)
plt.xlabel('$t$ (sec)')
plt.ylabel('$n$ (photons)')
plt.title('Haken laser model')
#plt.legend(('$n$ (photons)', '$\dot{n}$ (photons/sec)'))
#plt.show()


n = np.linspace(0., 10., 100).tolist()
nd = list(map(makeHaken(100), n))
zero = list(map(lambda x:0, n))

#plt.figure(2)
plt.subplot(122)
plt.plot(n, nd)
plt.plot(n, zero)
plt.xlabel('$n$')
plt.ylabel('$\dot{n}$')
plt.title('Haken $n$ vs $\dot{n}$')
plt.show()
'''
    
