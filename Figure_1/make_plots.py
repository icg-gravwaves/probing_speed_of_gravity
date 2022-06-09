import matplotlib.pyplot as plt
import math

import numpy as np
import sys
import os
from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import corner
import matplotlib.lines as mlines
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import glob

plt.style.use(['seaborn-talk','paper'])

### FIGURE 1

plotrange = [[-1.0,   3.0],[-1.5,   6.5]]
params = {'labels': [r"$c_B$", r"$c_M$"]}

lwidth_cont=1.2
lwidth_lines=0.8

fig = plt.figure()
ax = fig.add_subplot(111)

x = np.linspace(-5,5,1000)
# the function
y2 = (1/2.)*(1+np.tanh(3*x))
y3 = (1/2.)*(1+np.tanh(6*x))
y = (1/2.)*(1+np.tanh(x))

alpha=2 # Note that we changed nomenclature during writing the paper. This
        # alpha is actually the sigma parameter.
beta=2
gamma=4

y4=1/4.*(1+np.tanh(gamma*x))*(1+np.tanh(beta*x))+1/4.*(1-np.tanh(gamma*x))*(1+np.tanh(alpha*x))

plt.plot(x,y, label=r"$\sigma=1$")
plt.plot(x,y2, linestyle='dashed', label=r"$\sigma=3$")
plt.plot(x,y3, linestyle='dotted', label=r"$\sigma=6$")
plt.legend(loc=4, frameon=True)

ax.set_xticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$10^{-2} f_\star$',r'$10^{-1} f_\star$',r'$f_\star$',r'$10 f_\star$',r'$10^2 f_\star$'])

ax.set_yticks([0,1])
ax.set_yticklabels([r'$c_{\rm gw}(0)$',r'c'])

cmap = '#555555'
ax.plot([-3.,3.], [0.,0.],color=cmap,linestyle=':',linewidth=lwidth_lines)
ax.plot([0.,0.], [-10.,10.],color=cmap,linestyle=':',linewidth=lwidth_lines)
ax.plot([1.,1.], [-10.,10.],color=cmap,linestyle=':',linewidth=lwidth_lines)
ax.plot([2.,2.], [-10.,10.],color=cmap,linestyle=':',linewidth=lwidth_lines)
ax.plot([-1.,-1.], [-10.,10.],color=cmap,linestyle=':',linewidth=lwidth_lines)
ax.plot([-2.,-2.], [-10.,10.],color=cmap,linestyle=':',linewidth=lwidth_lines)

ax.plot([-3.,3.], [1.,1.],color=cmap,linestyle=':',linewidth=lwidth_lines)
ax.plot([-3.,3.], [-1.,-1.],color=cmap,linestyle=':',linewidth=lwidth_lines)

plt.xlim(-2.2,2.2)
plt.ylim(-0.2,1.2)

#ax=axes[1,1]

#ax.set_xscale('log')

#ax.set_xlabel('frequency', size=14)
#ax.set_ylabel(r'$c_{\rm gw}$', size=20)
#ax.tick_params(labelsize=16)

#ax.set_title(r'$c_{\rm gw}$ vs. frequency')

#line1 = mlines.Line2D([], [], color=cmapT4[0], label=r"$\alpha=1$")
#line2 = mlines.Line2D([], [], color=cmapT4[1], linestyle='dashed', label=r"$\alpha=3$")
#line3 = mlines.Line2D([], [], color=cmapT4[3], linestyle='dotted', label=r"$\alpha=6$")

#plt.show()

plt.savefig("figure1.pdf")
plt.savefig("figure1.png")

#### APPENDIX FIGURE

cmap1 = '#555555'

plotrange = [[-1.0,   3.0],[-1.5,   6.5]]
params = {'labels': [r"$c_B$", r"$c_M$"]}

lwidth_cont=1.2
lwidth_lines=0.8

fig = plt.figure()
ax = fig.add_subplot(111)

x = np.linspace(-5,5,1000)
# the function
y2 = (1/2.)*(1+np.tanh(3*x))
y3 = (1/2.)*(1+np.tanh(6*x))
y = (1/2.)*(1+np.tanh(x))

alpha=1
beta=1
gamma=6

y4=1/4.*(1+np.tanh(gamma*x))*(1+np.tanh(beta*x))+1/4.*(1-np.tanh(gamma*x))*(1+np.tanh(alpha*x))


alpha=3
beta=1
gamma=6

y5=1/4.*(1+np.tanh(gamma*x))*(1+np.tanh(beta*x))+1/4.*(1-np.tanh(gamma*x))*(1+np.tanh(alpha*x))

alpha=1
beta=3
gamma=6

y6=1/4.*(1+np.tanh(gamma*x))*(1+np.tanh(beta*x))+(1/4.)*(1-np.tanh(gamma*x))*(1+np.tanh(alpha*x))

plt.plot(x,y4, label=r"$\sigma=1, \beta=1$")
plt.plot(x,y5, linestyle='dashed', label=r"$\sigma=3, \beta = 1$")
plt.plot(x,y6, linestyle='dotted', label=r"$\sigma=1, \beta = 3$")

plt.legend(loc=4, frameon=True)

ax.set_xticks([-2,-1,0,1,2])
ax.set_xticklabels([r'$10^{-2} f_\star$',r'$10^{-1} f_\star$',r'$f_\star$',r'$10 f_\star$',r'$10^2 f_\star$'])

ax.set_yticks([0,1])
ax.set_yticklabels([r'$c_{\rm gw}(0)$',r'c'])

ax.plot([-3.,3.], [0.,0.],color=cmap1,linestyle=':',linewidth=lwidth_lines)
ax.plot([0.,0.], [-10.,10.],color=cmap1,linestyle=':',linewidth=lwidth_lines)
ax.plot([1.,1.], [-10.,10.],color=cmap1,linestyle=':',linewidth=lwidth_lines)
ax.plot([2.,2.], [-10.,10.],color=cmap1,linestyle=':',linewidth=lwidth_lines)
ax.plot([-1.,-1.], [-10.,10.],color=cmap1,linestyle=':',linewidth=lwidth_lines)
ax.plot([-2.,-2.], [-10.,10.],color=cmap1,linestyle=':',linewidth=lwidth_lines)

ax.plot([-3.,3.], [1.,1.],color=cmap1,linestyle=':',linewidth=lwidth_lines)
ax.plot([-3.,3.], [-1.,-1.],color=cmap1,linestyle=':',linewidth=lwidth_lines)

plt.xlim(-2.2,2.2)
plt.ylim(-0.2,1.2)

plt.savefig("figureA1.pdf")
plt.savefig("figureA1.png")
