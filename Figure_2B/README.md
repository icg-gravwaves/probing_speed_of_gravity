# Figure 2B

Here we describe how to reproduce Figure 2B

## Step 1: Run the PyCBC Inference Bayesian analysis code

We use `pycbc_inference` to make the data for this figure, but it first needs to know about our stretched/squeezed waveform as described in the manuscript.

To do this we add the following function into `pycbc/waveform/waveform.py`

```
from scipy import interpolate
MPC_S_SI = lal.PC_SI * 1E6 / lal.C_SI
from pycbc.waveform.utils import td_taper
def get_dopplered_tf2_waveform(dtype=numpy.complex128, return_hc=True,
                               is_sequence=False,
                               **params):
    cparams = copy.deepcopy(params)
    cparams['approximant'] = 'TaylorF2'
    if is_sequence:
        hp, hc = get_fd_waveform_sequence(**cparams)
    else:
        hp, hc = get_fd_waveform(**cparams)
    non_gr_fac_n = params['non_gr_fac_n']
    non_gr_fac_a_at_30 = params['non_gr_fac_a_at_30']
    non_gr_fac_a = non_gr_fac_a_at_30 * 30**non_gr_fac_n
    tsize = int(1.0 / params['delta_t'] /  params['delta_f'])
    fsize = tsize // 2 + 1
    hp.resize(fsize)
    td_wav = hp.to_timeseries()
    td_wav.roll(- int(len(td_wav) * 0.2))
    sts_shifted = numpy.arange(len(td_wav), dtype=float) * td_wav.delta_t + float(td_wav._epoch)
    sts_shifted += int(len(td_wav) * 0.2) * td_wav.delta_t
    sts_shifted_red = sts_shifted[sts_shifted<0]
    mc = conversions.mchirp_from_mass1_mass2(params['mass1'],
                                             params['mass2'])
    f_t_basic = (256. / 5. * (mc * lal.MTSUN_SI)**(5./3.) * numpy.pi**(8./3.)
                     * (-sts_shifted_red))**(-3./8.)
    new_times = (sts_shifted_red - non_gr_fac_a * 40 * MPC_S_SI
                     * f_t_basic**(-non_gr_fac_n))
    new_times = numpy.concatenate([new_times,sts_shifted[sts_shifted >= 0]])
    interfunc = interpolate.interp1d(new_times, td_wav.data,
                                     bounds_error=False,
                                     fill_value='extrapolate')
    new_data = interfunc(sts_shifted)
    td_wav.data[:] = new_data[:]
    td_wav = td_taper(td_wav, -len(td_wav) * td_wav.delta_t, -0.9 * len(td_wav) * td_wav.delta_t, side='left')
    td_wav = td_taper(td_wav, -0.1*len(td_wav) * td_wav.delta_t, 0, side='right')
    td_wav.roll(+ int(len(td_wav) * 0.2))

    hp = td_wav.to_frequencyseries()
    hp = hp.astype(dtype)
    
    if return_hc:
        hc.resize(fsize)
        td_wav = hc.to_timeseries()
        td_wav.roll(- int(len(td_wav) * 0.2))
        interfunc = interpolate.interp1d(new_times, td_wav.data, 
                                         bounds_error=False, 
                                         fill_value='extrapolate')
        new_data = interfunc(sts_shifted)
        td_wav.data[:] = new_data[:]
        td_wav = td_taper(td_wav, -len(td_wav) * td_wav.delta_t, -0.9 * len(td_wav) * td_wav.delta_t, side='left')
        td_wav = td_taper(td_wav, -0.1*len(td_wav) * td_wav.delta_t, 0, side='right')
        td_wav.roll(+ int(len(td_wav) * 0.2))

        hc = td_wav.to_frequencyseries()
        hc = hc.astype(dtype)
    else:
        hc = None

    return hp, hc
```

and declare it as a waveform by adding:

```
    if apx == 'TaylorF2':
        apx_new = 'TaylorF2_DOPPLER'
        cpu_fd[apx_new] = get_dopplered_tf2_waveform
```

at this line

https://github.com/gwastro/pycbc/blob/v2.0.3/pycbc/waveform/waveform.py#L1003

This will correspond to superluminal cgw. We simulate subluminal cgw by changing the line:

```
    new_times = (sts_shifted_red - non_gr_fac_a * 40 * MPC_S_SI
                     * f_t_basic**(-non_gr_fac_n))
```

to

```
    new_times = (sts_shifted_red + non_gr_fac_a * 40 * MPC_S_SI
                     * f_t_basic**(-non_gr_fac_n))
```

(doubtless this would need to be more elegant if being run on production LVK/LISA data!).

Then we need a configuration file, ours can be downloaded here:

https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/inference_emcee.ini

Then we need to run `pycbc_inference`. We did this using MPI according to:

```
mpiexec -n 50 pycbc_inference --config-files inference_emcee.ini --output-file results_superluminal.hdf --verbose --use-mpi
```

This produces `results_superluminal.hdf` and `results_subluminal.hdf`, which we upload to the repository here:


https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/results_superluminal.hdf


https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/results_subluminal.hdf

Note that this is slow (takes weeks on 50 cores). The waveform generator should be made faster is using this in the future.

## Step 2: Make the plot

We then make the plot in a few stages. We first extract the necessary data from the HDF files. This is done by downloaded the edited executable:

https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/pycbc_inference_plot_posterior

and then running

```
./pycbc_inference_plot_posterior --input-file results_superluminal.hdf --plot-scatter --output-file scatter_2param_mk2.pdf --z-arg snr --parameters non_gr_fac_n "log10(non_gr_fac_a_at_30)"
mv figure2_data.npz figure2_data_superluminal.npz

./pycbc_inference_plot_posterior --input-file results_subluminal.hdf --plot-scatter --output-file scatter_2param_mk2.pdf --z-arg snr --parameters non_gr_fac_n "log10(non_gr_fac_a_at_30)"

mv figure2_data.npz figure2_data_subluminal.npz
```
to generate two .npz files. These files are also available here:

https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/figure2_data_subluminal.npz

https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/figure2_data_superluminal.npz

Finally the plot is generated from this be running the script

https://icg-gravwaves.github.io/ian_harry/probing_speed_of_gravity/Figure_2B/make_plot.py

## Step 3: Also make the video

The `make_plot.py` command also makes the still frames used to make the animation. We combined this into an mp4 using ffmpeg:

```
ffmpeg -r 1 -f image2 -i figure2b_video_panel_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" -start_number 0 figure2b.mp4
```
