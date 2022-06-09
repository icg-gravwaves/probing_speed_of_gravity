import matplotlib
import numpy
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy import ndimage

plt.style.use(['seaborn-talk','paper'])

fstar_values = numpy.logspace(-5,1,25)
alpha_values = numpy.linspace(0.2,10,50)

zoom_level = 10

fstar_valuesZ = numpy.logspace(-5,1,25*zoom_level)
alpha_valuesZ = numpy.linspace(0.2,10,50*zoom_level)

data_LIGO = numpy.loadtxt('LIGO_output.txt', delimiter=',')
data_LIGO = data_LIGO[data_LIGO[:,2] < 0]
data_LIGO = numpy.array(sorted(data_LIGO, key = lambda x: x[1]))
data_LIGO = numpy.array(sorted(data_LIGO, key = lambda x: x[0]))

data_LISA = numpy.loadtxt('LISA_output.txt', delimiter=',')
data_LISA = data_LISA[data_LISA[:,2] < 0]
data_LISA = numpy.array(sorted(data_LISA, key = lambda x: x[1]))
data_LISA = numpy.array(sorted(data_LISA, key = lambda x: x[0]))


f_stars = data_LISA[:,0]
alphas = data_LISA[:,1]
exclusions = data_LISA[:,3]
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=[17,6])


ax1 = axes[0]
x_1, y_1 = numpy.meshgrid(fstar_valuesZ, alpha_valuesZ)
sc_plot = ax1.scatter(f_stars, alphas, s=15, marker='o', c=exclusions, norm=matplotlib.colors.LogNorm(vmax=1, vmin=1E-20))
log_grid_excs = numpy.log10(numpy.resize(exclusions,[25,50]).T)
result = ndimage.zoom(log_grid_excs, zoom_level, order=3, prefilter=False)
result = 10**result
ctours = ax1.contour(fstar_valuesZ, alpha_valuesZ, result, levels=[1E-15, 1E-12, 1E-9, 1E-6, 1E-3], norm=matplotlib.colors.LogNorm(vmax=1, vmin=1E-20))
#ax1.clabel(ctours, inline=True, fontsize=20, fmt='%.0e')
ax1.semilogx()
#cb = fig.colorbar(sc_plot, ax=ax)
ax1.set_ylabel(r'$\sigma$')
ax1.set_xlabel(r'$f_{\star}$')
ax1.set_xticks([1E-5, 1E-3, 1E-1, 10])
#cb.set_label('$\delta$ $c_{GW}$ exclusion')

ax2 = axes[1]
f_stars = data_LIGO[:,0]
alphas = data_LIGO[:,1]
exclusions = data_LIGO[:,3]
sc_plot = ax2.scatter(f_stars, alphas, s=15, marker='o', c=exclusions, norm=matplotlib.colors.LogNorm(vmax=1, vmin=1E-20))
log_grid_excs = numpy.log10(numpy.resize(exclusions,[25,50]).T)
result = ndimage.zoom(log_grid_excs, zoom_level, order=3, prefilter=False)
result = 10**result
ctours = ax2.contour(fstar_valuesZ, alpha_valuesZ, result, levels=[1E-15, 1E-12, 1E-9, 1E-6, 1E-3], norm=matplotlib.colors.LogNorm(vmax=1, vmin=1E-20))
ax2.semilogx()
#cb = fig.colorbar(sc_plot, ax=ax)
#ax2.set_ylabel(r'$\alpha$')
ax2.set_xlabel(r'$f_{\star}$')
ax2.set_xticks([1E-5, 1E-3, 1E-1, 10])
ax2.tick_params('y', labelleft=False)

#cb.set_label('$\delta$ $c_{GW}$ exclusion')

f_stars = data_LIGO[:,0]
alphas = data_LIGO[:,1]

ax3 = axes[2]
exclusions1 = data_LIGO[:,3]
exclusions2 = data_LISA[:,3]
exclusions_comb = numpy.amin(numpy.array([exclusions1, exclusions2]), axis=0)
print(exclusions_comb)
sc_plot = ax3.scatter(f_stars, alphas, s=15, marker='o', c=exclusions_comb, norm=matplotlib.colors.LogNorm(vmax=1, vmin=1E-20))

log_grid_excs = numpy.log10(numpy.resize(exclusions_comb,[25,50]).T)
result = ndimage.zoom(log_grid_excs, zoom_level, order=3, prefilter=False)
result = 10**result
ctours = ax3.contour(fstar_valuesZ, alpha_valuesZ, result, levels=[1E-15, 1E-12, 1E-9, 1E-6, 1E-3], norm=matplotlib.colors.LogNorm(vmax=1, vmin=1E-20))

ax3.semilogx()
#ax3.set_ylabel(r'$\alpha$')
ax3.set_xlabel(r'$f_{\star}$')
ax3.set_xticks([1E-5, 1E-3, 1E-1, 10])
ax3.tick_params('y', labelleft=False)

#cb = ax3.cax.colorbar(sc_plot)
#cbar = grid.cbar_axes[0].colorbar(sc_plot)
#cbar.set_label('$\delta$ $c_{GW}$ exclusion')

fig.subplots_adjust(right=0.93)
cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
cbar = fig.colorbar(sc_plot, cax=cbar_ax)
cbar.set_label('$\delta$ $c_{GW}$ exclusion')

fig.savefig('figure3_sl.pdf', bbox_inches='tight')
fig.savefig('figure3_sl.png', bbox_inches='tight')

