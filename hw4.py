#from matplotlib.pyplot import figure, plot, xlabel, ylabel,\
#title, xlim, ylim, legend, savefig, show, xscale, yscale
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

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

# Problem 1 ###################################################################
#cases = [] # (k1,  kn1, k2,  kn2, k3, title)
#cases.append((1.1, 0.0, 2.1, 0.0, 10.0, 'Irreversible'))
#cases.append((1.1, 1.1, 2.1, 2.1, 10.0, 'Fully Reversible'))
#cases.append((0.1, 1.1, 0.1, 1.1, 10.0, 'Slow Reversible'))
#cases.append((1.1, 0.1, 2.1, 0.1, 10.0, 'Fast Reversible'))
#
#fig1 = plt.figure(2, figsize=(18, 12))
#numcases = len(cases)
#axes = []
#for i in range(1,numcases+1):
#    axes.append(fig1.add_subplot(1, numcases, i))
#for (axisnum, (k1, kn1, k2, kn2, k3, titlestring)) in enumerate(cases):
#    print "Axis %i of %i: %s" % (axisnum, numcases, titlestring)
#    ICs = make_ICs(resolution=16)
#    for (x0, y0) in ICs:
#        (tlist, xlist, ylist, plist) = do_IC(k1, kn1, k2, kn2, k3, x0, y0)
#        axes[axisnum].plot(xlist, ylist, color='k')
#
#    axes[axisnum].set_title(titlestring)
#    if axisnum==0:
#        axes[axisnum].set_y(r'Adsorbed $O_1$')
#        axes[axisnum].set_xlabel(r'Adsorbed $CO_1$')
#plt.suptitle(r'504 HW4 - Tom Bertalan' + '\n' + \
#r'$k_2=%.2f$, $k_{-2}=%.2f$, $k_3=%.2f$' % (k2, kn2, k3))
##savefig('p1-%i_cases.png' % (numcases))
##fig2.tight_layout()
#plt.show()

# Problem 2 ###################################################################
#cases = []
#k1change=.25
#kn1=0
#for k1 in np.arange(0.0, 2, k1change):
#    cases.append((k1, kn1, 1.0, 0.0, 10.0, r'Irreversible, $k_1=%.2f$, $k_{-1}=%.2f$' % (k1, kn1)))
#kn1=.05
#for k1 in np.arange(0.0, 2, k1change):
#    cases.append((k1, kn1, 1.0, 0.0, 10.0, r'Reversible, $k_1=%.2f$, $k_{-1}=%.2f$' % (k1, kn1)))
#
#fig2 = plt.figure(2, figsize=(18, 12))
##fig2 = plt.figure(2, figsize=(24, 16))
#numcases = len(cases)
#axes = []
#for i in range(1,numcases+1):
#    axes.append(fig2.add_subplot(2, numcases/2, i))
#for (axisnum, (k1, kn1, k2, kn2, k3, titlestring)) in enumerate(cases):
#    print "Axis %i of %i: %s" % (axisnum, numcases, titlestring)
#    ICs = make_ICs(resolution=64)
#    for (x0, y0) in ICs:
#        (tlist, xlist, ylist, plist) = do_IC(k1, kn1, k2, kn2, k3, x0, y0)
#        axes[axisnum].plot(xlist, ylist, color='k')
#
#    axes[axisnum].set_title(titlestring)
#    if axisnum==0:
#        axes[axisnum].set_ylabel(r'Adsorbed $O_1$')
#    if axisnum==numcases/2:
#        axes[axisnum].set_xlabel(r'Adsorbed $CO_1$')
#plt.suptitle(r'504 HW4 - Tom Bertalan' + '\n' + \
#r'$k_2=%.2f$, $k_{-2}=%.2f$, $k_3=%.2f$' % (k2, kn2, k3))
##savefig('p2-%i_cases.png' % (numcases))
##fig2.tight_layout()
#plt.show()


# Problem 2 more
k2 = 1
k3 = 10
fk1 = lambda y: k2*k3*(1-y)/k3/y
fx = lambda y: fk1(y)*(1-y)/(fk1(y)+k3*y)
fr = lambda x,y: k3*x*y
l1 = list(np.arange(0.001,        1,   .001))
l2 = list(np.arange(0.0001,    .001,  .0001))
l3 = list(np.arange(0.00001,  .0001, .00001))
l4 = list(np.arange(0.000001,.00001,.000001))
ylist = l1
ylist.extend(l2)  # We need extra resolution at low y
ylist.extend(l3)  # to properly capture the behavior.
ylist.extend(l4)  # Only a couple dozen extra points.
ylist.sort()
xlist = []
k1list= []
rlist = []

kn1 = 0.05
fk1n = lambda y: (kn1+k3*y)*(sqrt(1+4*k2*(1-y)/k3/y)-1)/2
fxn = lambda y: fk1n(y)*(1-y)/(fk1n(y)-kn1+k3*y)
frn = fr
ylistn = []
xlistn = []
k1listn= []
rlistn = []
for y in ylist:
    x = fx(y)
    xlist.append(x)
    k1list.append(fk1(y))
    rlist.append(fr(x,y))

    ylistn.append(y)
    xn = fxn(y)
    xlistn.append(xn)
    k1listn.append(fk1n(y))
    rlistn.append(frn(xn,y))

# Make a plot comparing reversible and reversible
fig2b = plt.figure(22, figsize=(12, 8))

ax2O = fig2b.add_subplot(3,1,1)
ax2O.set_xlim((0,1))
ax2O.plot(k1list, ylist, 'k--')
ax2O.plot(k1listn, ylistn, 'k-')
ax2O.set_ylabel(r'Adsorbed $O_1$')

ax2CO = fig2b.add_subplot(3,1,2)
ax2CO.set_xlim((0,1))
ax2CO.plot(k1list,xlist, 'k--')
ax2CO.plot(k1listn,xlistn, 'k-')
ax2CO.set_ylabel(r'Adsorbed $CO_1$')

ax2r = fig2b.add_subplot(3,1,3)
ax2r.set_xlim((0,1))
ax2r.set_ylim((0,.5))
ax2r.plot(k1list, rlist, 'k--')
ax2r.plot(k1listn, rlistn, 'k-')
ax2r.set_ylabel(r'Rate of Product Formation, ${d[CO_2]}/{dt}$')
ax2r.set_xlabel(r'Partial Pressure of Unadsorbed $CO_2$, $k_1$')

ax2r.legend([r'$k_{-1}=0$, Irreversible', r'$k_{-1}=%.2f$, Reversible'%kn1])
plt.suptitle('Reversible vs. Irreversible Adsorption and Reaction\n'\
             r'$k_2=%.2f$, $k_3=%.2f$, $k_{-2}=0$' % (k2, k3))
plt.savefig('hw4_2b.pdf')
#plt.show()

# Make a plot showing irreversible only, since it has
# qualitatively interesting large-k1 behavior that won't
# fit on the previous plot
fig2c = plt.figure(23, figsize=(12, 8))

ax2O = fig2c.add_subplot(3,1,1)
ax2O.set_xlim((0,100))
ax2O.plot(k1list, ylist, 'k--')
ax2O.set_ylabel(r'Adsorbed $O_1$')

ax2CO = fig2c.add_subplot(3,1,2)
ax2CO.set_xlim((0,100))
ax2CO.plot(k1list,xlist, 'k--')
ax2CO.set_ylabel(r'Adsorbed $CO_1$')

ax2r = fig2c.add_subplot(3,1,3)
ax2r.set_xlim((0,100))
ax2r.plot(k1list, rlist, 'k--')
ax2r.set_ylabel(r'Rate of Product Formation, ${d[CO_2]}/{dt}$')
ax2r.set_xlabel(r'Partial Pressure of Unadsorbed $CO_2$, $k_1$')

ax2r.legend([r'$k_{-1}=0$, Irreversible'])
plt.suptitle('Irreversible Adsorption and Reaction\n'\
             r'$k_2=%.2f$, $k_3=%.2f$, $k_{-2}=0$' % (k2, k3))
plt.savefig('hw4_2c.pdf')
#plt.show()


# Problem 3 ###################################################################
k1=.5
kn1=.05
k2=1
kn2=0
k3=10
# roots of the cubic polynomial in y that results from setting
#    dx/dt = 0 = dy/dt
coeffs = [
          -k2*k3**2,
          -2*k2*kn1*k3 + k2*k3**2 - k1*k3**2,
          -k1**2*k3 - k1*k3*kn1 + 2*k2*k3*kn1 - k2*kn1**2,
          k2*kn1**2
         ]

print sum(coeffs)
yroots = np.roots(coeffs)
print "The roots are y =", yroots
fx = lambda y: k1*(1-y)/(k1+kn1+k3*y)
xroots = map(fx, yroots)
print "Similarly,    x =", xroots

fig3 = plt.figure(3, figsize=(12, 8))
ax3 = fig3.add_subplot(1,1,1)
ICs_ss = zip(xroots, yroots)
ICs_ss.append((0, 1))

for (x0,y0) in ICs_ss:
    x0 = x0 + np.random.rand()*.25
    y0 = y0 + np.random.rand()*.25
    (tlist, xlist, ylist, plist) = do_IC(k1, kn1, k2, kn2, k3, x0, y0)
    ax3.plot(xlist, ylist, 'k')
    ax3.scatter([x0],[y0], color='k')
ax3.set_title('Small random perturbations from steady-states')
ax3.set_xlabel(r'x, adsorbed $CO_1$')
ax3.set_ylabel(r'y, adsorbed $O_1$')
ax3.set_xbound(lower=-.08)
ax3.set_ybound(lower=0)
plt.savefig('hw4_3.pdf')

fig3b = plt.figure(32, figsize=(12, 8))
ax3b = fig3b.add_subplot(1,1,1)
ICs = make_ICs()

for (x0,y0) in ICs:
    (tlist, xlist, ylist, plist) = do_IC(k1, kn1, k2, kn2, k3, x0, y0)
    ax3b.plot(xlist, ylist, 'k')
    ax3b.set_xlim((0,1))
    ax3b.set_ylim((0,1))
ax3b.set_xlim((0, 1))
ax3b.set_ylim((0, 1))
ax3b.set_title('Phase plot for many ICs.')
ax3b.set_xlabel(r'x, adsorbed $CO_1$')
ax3b.set_ylabel(r'y, adsorbed $O_1$')
plt.savefig('hw4_3b.pdf')

# Problem 3, Linear stability analysis ########################################
f1x = lambda x,y: -k1 -kn1 - k3*y
f2x = lambda x,y: -2*k2+2*k2*x+2*k2*y-k3*y
f1y = lambda x,y: -k1-k3*x
f2y = lambda x,y: -2*k2+2*k2*x+2*k2*y-k3*x

for (x,y) in ICs_ss:
    jacobian = np.array([[f1x(x,y), f1y(x,y)],
                          [f2x(x,y), f2y(x,y)]])
    print ''
    print 'At (', x, ',', y, '), the Jacobian is:'
    print jacobian
    [a,b,c,d] = jacobian.flatten().tolist()
    coeffs = [1, -d-a, a*d-c*b]
    eigv = np.roots(coeffs)
    print 'Eigenvalues are:', eigv


