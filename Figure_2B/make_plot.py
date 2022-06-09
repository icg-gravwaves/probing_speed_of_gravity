import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['seaborn-talk','paper'])

data_super = np.load('figure2_data_superluminal.npz')
data_sub = np.load('figure2_data_subluminal.npz')
non_gr_n_values_sub = data_sub['non_gr_n_values']
non_gr_deltacgw_sub = data_sub['non_gr_deltacgw']
snr_sub = data_sub['snr']

non_gr_n_values_super = data_super['non_gr_n_values']
non_gr_deltacgw_super = data_super['non_gr_deltacgw']
snr_super = data_super['snr']

for idx, freq in enumerate([30,40,50,60,80,100,150,200,250,300,500,700,1000]):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(True)
    non_gr_deltacgwT = np.log10(10**non_gr_deltacgw_sub * (30. / freq)**non_gr_n_values_sub)
    im = ax.scatter(non_gr_n_values_sub/2, 10**non_gr_deltacgwT, s=6,
                    marker='o')
    non_gr_deltacgwT = np.log10(10**non_gr_deltacgw_super * (30. / freq)**non_gr_n_values_super)
    im = ax.scatter(non_gr_n_values_super/2, 10**non_gr_deltacgwT, s=6,
                    marker='s')

    ax.semilogy()
    ax.set_xlim(0.5,5)
    ax.set_ylim(1E-20, 1E-16)
    ax.set_xlabel(r'$\sigma$: power law slope parameter')
    ax.set_ylabel(r'$\delta c_{\rm GW}(' + str(freq) + r')$: Deviation from $c$ at ' + str(freq) + 'Hz')
    ax.set_title(r'\textbf{Exclusion from GW170817 data}', fontweight='bold')
    ax.legend(
        [
            r'Subluminal $\delta c_{\rm{GW}}$',
            r'Superluminal $\delta c_{\rm{GW}}$'
        ], 
        markerscale=3
    )
    fig.savefig(f'figure2b_video_panel_{idx}.png', bbox_inches='tight')
    if freq == 30:
        fig.savefig('figure2b.pdf')
        fig.savefig('figure2b.png')

