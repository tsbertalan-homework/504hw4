#from matplotlib.pyplot import figure, plot, xlabel, ylabel,\
#title, xlim, ylim, legend, savefig, show, xscale, yscale
import numpy as np
import matplotlib.pyplot as plt

#for i in range(ICs):
#    [xold, yold] = np.random.rand(2)

def do_IC(k1, kn1, k2, kn2, k3, x0, y0, dt=0.01, numsteps=256):

    def xprime(x, y):
        return k1 * (1 - x - y) - kn1 * x - k3 * x * y

    def yprime(x, y):
        return k2 * (1 - x - y) - kn2 * y ** 2 - k3 * x * y

    def pprime(x, y):
        return k3 * x * y

    xold = x0
    yold = y0
    tlist = []
    xlist = []
    ylist = []
    plist = []
    pold = 0
    told = 0
    tlist.append(told)
    xlist.append(xold)
    ylist.append(yold)
    plist.append(pold)
    for j in range(numsteps):
        x = xold + dt * xprime(xold, yold)
        y = yold + dt * yprime(xold, yold)
        p = pold + dt * pprime(xold, yold)
        t = told + dt
        xlist.append(x)
        xold = x
        ylist.append(y)
        yold = y
        plist.append(p)
        pold = p
        tlist.append(t)
        told = t
    return (tlist, xlist, ylist, plist)

def make_ICs(resolution=32):
    ICs = []
    for x0 in np.arange(0, 1. + 1./resolution, 1./resolution):
        ICs.append((x0, 0))
        ICs.append((x0, 1))
        ICs.append((0, x0))
        ICs.append((1, x0))
    return ICs

cases = [] # (k1,  kn1, k2,  kn2, k3, title)
#cases.append((1.0, 0.0, 1.0, 0.0, 10.0, 'Irreversible'))
#cases.append((1.0, 1.0, 1.0, 1.0, 10.0, 'Fully Reversible'))
#cases.append((0.1, 1.0, 0.1, 1.0, 10.0, 'Slow'))

k1change=.5
for k1 in np.arange(0.0, 1.0, k1change):
    cases.append((k1, 0.0, 1.0, 0.0, 10.0, r'Irreversible, $k_1=%.2f$, $k_{-1}=%.2f$' % (k1, 0)))
for k1 in np.arange(0.0, 1.0, k1change):
    cases.append((k1, 01, 1.0, 0.0, 10.0, r'Reversible, $k_1=%.2f$, $k_{-1}=%.2f$' % (k1, .5)))


fig2 = plt.figure(2, figsize=(18, 12))
#fig2 = plt.figure(2, figsize=(24, 16))
numcases = len(cases)
axes = []
for i in range(1,numcases+1):
    axes.append(fig2.add_subplot(2, numcases/2, i))
for (axisnum, (k1, kn1, k2, kn2, k3, titlestring)) in enumerate(cases):
    print "Axis %i of %i: %s" % (axisnum, numcases, titlestring)
    ICs = make_ICs(resolution=64)
    for (x0, y0) in ICs:
        (tlist, xlist, ylist, plist) = do_IC(k1, kn1, k2, kn2, k3, x0, y0)
        axes[axisnum].plot(xlist, ylist, color='k')

    axes[axisnum].set_title(titlestring)
    if axisnum==0:
        axes[axisnum].set_ylabel(r'Adsorbed $O_1$')
    if axisnum==numcases/2:
        axes[axisnum].set_xlabel(r'Adsorbed $CO_1$')

plt.suptitle(r'504 HW4 - Tom Bertalan' + '\n' + \
r'$k_2=%.2f$, $k_{-2}=%.2f$, $k_3=%.2f$' % (k2, kn2, k3))
#savefig('p2-%i_cases.png' % (numcases))
#fig2.tight_layout()
plt.show()


