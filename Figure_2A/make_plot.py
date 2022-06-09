import numpy as np
import pycbc
import pycbc.waveform
import pycbc.conversions
import pycbc.psd
import pycbc.filter
import lal
from scipy import interpolate
import sys
import copy
from pycbc.waveform import get_fd_waveform
from pycbc.waveform.utils import td_taper
from pycbc import conversions, cosmology


import matplotlib.pyplot as plt

plt.style.use(['seaborn-talk','paper'])

dataLIGO = np.load('figure2_dataLIGO.npz')
dataET = np.load('figure2_dataET.npz')
non_gr_n_values_lg = dataLIGO['non_gr_n_values']
outputs_lg = dataLIGO['outputs']
non_gr_n_values_et = dataET['non_gr_n_values']
outputs_et = dataET['outputs']
print(non_gr_n_values_et, outputs_et)

#plt.rc('font', family='serif')
#plt.rc('xtick', labelsize='x-small')
#plt.rc('ytick', labelsize='x-small')

fig = plt.figure()#figsize=(4, 3))
ax = fig.add_subplot(111)
ax.semilogy(non_gr_n_values_lg/2, outputs_lg, label=r'aLIGO; $f_{\mathrm{ref}}$=30Hz', color='#1f77b4')
outputs_lg_100 = outputs_lg * (30. / 100)**non_gr_n_values_lg
ax.semilogy(non_gr_n_values_lg/2, outputs_lg_100, label=r'aLIGO; $f_{\mathrm{ref}}$=100Hz', color='#ff7f0e')
outputs_lg_500 = outputs_lg * (30. / 500)**non_gr_n_values_lg
ax.semilogy(non_gr_n_values_lg/2, outputs_lg_500, label=r'aLIGO; $f_{\mathrm{ref}}$=500Hz', color='#2ca02c')

ax.semilogy(non_gr_n_values_et/2, outputs_et, label=r'ET; $f_{\mathrm{ref}}$=30Hz', color='#1f77b4', linestyle='--')
outputs_et_100 = outputs_et * (30. / 100)**non_gr_n_values_et
ax.semilogy(non_gr_n_values_et/2, outputs_et_100, label=r'ET; $f_{\mathrm{ref}}$=100Hz', color='#ff7f0e', linestyle='--')
outputs_et_500 = outputs_et * (30. / 500)**non_gr_n_values_et
ax.semilogy(non_gr_n_values_et/2, outputs_et_500, label=r'ET; $f_{\mathrm{ref}}$=500Hz', color='#2ca02c', linestyle='--')


#ax.semilogy(non_gr_n_values_et, outputs_et, label='ET')
ax.set_xlim(0,5)
ax.set_ylim(1E-24, 1E-18)
ax.set_xlabel(r'$\sigma$: power law slope parameter')
ax.set_ylabel(r'$\delta c_{\rm GW}(f_{\mathrm{ref}})$: Deviation from $c$ at $f_{\mathrm{ref}}$')
ax.legend()
ax.set_title(r'\textbf{Approximate exclusion capability}')
fig.savefig('figure2a.pdf', bbox_inches='tight')
fig.savefig('figure2a.png', bbox_inches='tight')

