##############################################################################
################ Code online at github.com/tsbertalan/504hw4 #################
##############################################################################
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

def do_IC(k1, kn1, k2, kn2, k3, x0, y0, dt=0.01, numsteps=256):
    def xprime(x, y):
        return k1 * (1 - x - y) - kn1 * x - k3 * x * y
    def yprime(x, y):
        return k2 * (1 - x - y)**2 - kn2 * y ** 2 - k3 * x * y
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

# Problem 2 ###################################################################
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
# coefficients of said polynomial:

coeffs = [
     k2*k3**2,
     k1*k3**2 - 2*k2*k3**2 + 2*k2*k3*kn1,
     k1**2*k3 - k1*k3**2 + k1*k3*kn1 + k2*k3**2 - 4*k2*k3*kn1 + k2*kn1**2,
     -k3*k1**2 - k1*k3*kn1 +2*k2*k3*kn1 - 2*k2*kn1**2,
     k2*kn1**2
     ]

yroots = np.roots(coeffs)
fx = lambda y: k1*(1-y)/(k1+kn1+k3*y)
xroots = map(fx, yroots)

fig3 = plt.figure(3, figsize=(12, 8))
ax3 = fig3.add_subplot(1,1,1)
ax3.set_xlim((0, 1))
ax3.set_ylim((0, 1.2))
ICs_ss = zip(xroots, yroots)
print '#### Problem 3, steady-states ####'
print "Steady-states (x, y) are:"
for (x, y) in ICs_ss:
    print '(%.3f, %.3f)' % (x, y)
repeats = 12
for (x0,y0) in ICs_ss * repeats:
    x0 = x0 + np.random.rand()*1/8. - 1/16.
    y0 = y0 + np.random.rand()*1/8. - 1/16.
    (tlist, xlist, ylist, plist) = do_IC(k1, kn1, k2, kn2, k3, x0, y0)
    ax3.plot(xlist, ylist, 'k')
for (x0,y0) in ICs_ss:
    ax3.scatter([x0], [y0], color='k')
ax3.set_title('Small random perturbations from steady-states')
ax3.set_xlabel(r'x, adsorbed $CO_1$')
ax3.set_ylabel(r'y, adsorbed $O_1$')
plt.savefig('hw4_3.pdf')

fig3b = plt.figure(32, figsize=(12, 8))
ax3b = fig3b.add_subplot(1,1,1)
ICs = make_ICs(resolution=128)

for (x0,y0) in ICs:
    (tlist, xlist, ylist, plist) = do_IC(k1, kn1, k2, kn2, k3, x0, y0)
    ax3b.plot(xlist, ylist, 'k')
    ax3b.set_xlim((0,1))
    ax3b.set_ylim((0,1))
for (x0,y0) in ICs_ss:
    ax3b.scatter([x0], [y0], color='k')
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
print ''
print '#### Problem 3, Linear stability analysis ####'
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


