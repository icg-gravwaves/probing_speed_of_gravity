import numpy as np
import pycbc
import pycbc.waveform
import pycbc.conversions
import pycbc.psd
import pycbc.filter
from pycbc.cosmology import redshift
import lal
from scipy import interpolate
import sys
import copy
from pycbc.waveform import get_fd_waveform
from pycbc.waveform.utils import td_taper
import sys

filenum = int(sys.argv[1])

file_name = f"output_{filenum}.txt"
file_pointer = open(file_name, "w")

MPC_S_SI = lal.PC_SI * 1e6 / lal.C_SI


def get_delta_t(f_array, **params):
    cgw_delta0 = params["non_gr_cgwdelta_0"]
    fstar = params["non_gr_transition_freq"]
    alpha = params["non_gr_alpha"]
    delta_cgw = cgw_delta0 + 0.5 * (cgw_delta0) * (
        np.tanh(alpha * np.log(f_array / fstar)) + 1
    )
    # print(delta_cgw)
    delta_t = delta_cgw * params["distance"] * MPC_S_SI
    delta_t -= delta_t[-1]
    # print(delta_t)
    return delta_t


def get_dopplered_tf2_waveform_lisa(
    dtype=np.complex128,
    return_hc=True,
    is_sequence=False,
    do_roll_and_taper=True,
    check_limits=False,
    **params,
):
    cparams = copy.deepcopy(params)
    cparams["approximant"] = "TaylorF2"
    if is_sequence:
        hp, hc = get_fd_waveform_sequence(**cparams)
    else:
        hp, hc = get_fd_waveform(**cparams)
    tsize = int(1.0 / params["delta_t"] / params["delta_f"])
    fsize = tsize // 2 + 1
    hp.resize(fsize)
    td_wav = hp.to_timeseries()
    if do_roll_and_taper:
        td_wav.roll(-int(len(td_wav) * 0.2))
    sts_shifted = np.arange(len(td_wav), dtype=float) * td_wav.delta_t + float(
        td_wav._epoch
    )
    if do_roll_and_taper:
        sts_shifted += int(len(td_wav) * 0.2) * td_wav.delta_t
    sts_shifted_red = sts_shifted[sts_shifted < 0]
    mc = pycbc.conversions.mchirp_from_mass1_mass2(params["mass1"], params["mass2"])
    # print(max(sts_shifted_red))
    f_t_basic = (
        256.0
        / 5.0
        * (mc * lal.MTSUN_SI) ** (5.0 / 3.0)
        * np.pi ** (8.0 / 3.0)
        * (-sts_shifted_red)
    ) ** (-3.0 / 8.0)
    new_times = sts_shifted_red + get_delta_t(f_t_basic, **params)
    new_times = np.concatenate([new_times, sts_shifted[sts_shifted >= 0]])
    if check_limits:
        print(max(new_times[sts_shifted < 0]))
    interfunc = interpolate.interp1d(
        new_times, td_wav.data, bounds_error=False, fill_value="extrapolate"
    )
    new_data = interfunc(sts_shifted)
    td_wav.data[:] = new_data[:]
    if do_roll_and_taper:
        td_wav = td_taper(
            td_wav,
            -len(td_wav) * td_wav.delta_t,
            -0.9 * len(td_wav) * td_wav.delta_t,
            side="left",
        )
        td_wav = td_taper(td_wav, -0.1 * len(td_wav) * td_wav.delta_t, 0, side="right")
        td_wav.roll(+int(len(td_wav) * 0.2))

    hp = td_wav.to_frequencyseries()
    hp = hp.astype(dtype)

    if return_hc:
        raise ValueError()
        hc.resize(fsize)
        td_wav = hc.to_timeseries()
        interfunc = interpolate.interp1d(
            new_times, td_wav.data, bounds_error=False, fill_value="extrapolate"
        )
        new_data = interfunc(td_wav.sample_times)
        td_wav.data[:] = new_data[:]
        hc = td_wav.to_frequencyseries()
        hc = hc.astype(dtype)
    else:
        hc = None

    return hp, hc


fstar_values = np.logspace(-5, 1, 25)
alpha_values = np.linspace(0.2, 10, 50)

curr_points = {}

sgastar_mass = 4.154e6
z_source = redshift(1000)
print(z_source, "redshift")
obs_mass = sgastar_mass * (1 + z_source)
print(obs_mass, "is the observed mass")


hp1, _ = get_dopplered_tf2_waveform_lisa(
    mass1=obs_mass,
    mass2=obs_mass,
    spin1z=0.0,
    spin2z=0.0,
    f_lower=5e-5,
    delta_f=8.0 / (33554432),
    delta_t=5,
    distance=1000,
    non_gr_cgwdelta_0=0,
    non_gr_transition_freq=1,
    non_gr_alpha=1,
    return_hc=False,
)
psd = pycbc.psd.from_txt("LISA_PSD.txt", len(hp1), hp1.delta_f, 2e-6, is_asd_file=False)
snr = pycbc.filter.sigma(hp1, psd=psd, low_frequency_cutoff=2e-6)
dist_overlap = 1 - 1 / (2 * snr)

print("Starting. Distinguishable overlap is", dist_overlap)

for fdx, fstar in enumerate([fstar_values[filenum]]):
    for adx, alpha in enumerate(alpha_values):
        for sign in [-1, 1]:

            def gen_waveform_and_overlap(curr_cgwdelta):
                hp2, _ = get_dopplered_tf2_waveform_lisa(
                    mass1=obs_mass,
                    mass2=obs_mass,
                    spin1z=0.0,
                    spin2z=0.0,
                    f_lower=5e-5,
                    delta_f=8.0 / (33554432),
                    delta_t=5,
                    distance=1000,
                    non_gr_cgwdelta_0=curr_cgwdelta_0 * sign,
                    non_gr_transition_freq=fstar,
                    non_gr_alpha=alpha,
                    return_hc=False,
                )
                overlap, _ = pycbc.filter.match(
                    hp1, hp2, psd=psd, low_frequency_cutoff=2e-6
                )
                return overlap

            print("Considering", fstar, alpha)
            curr_cgwdelta_0 = 1e-15
            overlap1 = gen_waveform_and_overlap(curr_cgwdelta_0)

            if overlap1 < dist_overlap:
                while overlap1 < dist_overlap:
                    curr_cgwdelta_0 = curr_cgwdelta_0 / 10
                    overlap1 = gen_waveform_and_overlap(curr_cgwdelta_0)
                    if curr_cgwdelta_0 < 1e-25:
                        break
                if curr_cgwdelta_0 < 1e-25:
                    file_pointer.write(f"{fstar}, {alpha}, {sign}, {curr_cgwdelta_0}\n")
                    continue
                low, high = curr_cgwdelta_0, curr_cgwdelta_0 * 10
            else:
                while overlap1 > dist_overlap:
                    curr_cgwdelta_0 = curr_cgwdelta_0 * 10
                    overlap1 = gen_waveform_and_overlap(curr_cgwdelta_0)
                    if curr_cgwdelta_0 > 10:
                        break
                if curr_cgwdelta_0 > 10:
                    file_pointer.write(f"{fstar}, {alpha}, {sign}, {curr_cgwdelta_0}\n")
                    continue
                low, high = curr_cgwdelta_0 / 10, curr_cgwdelta_0

            for i in range(8):
                mid = np.exp((np.log(low) + np.log(high)) / 2)
                overlap1 = gen_waveform_and_overlap(mid)
                if overlap1 < dist_overlap:
                    high = mid
                else:
                    low = mid
            mid = np.exp((np.log(low) + np.log(high)) / 2)
            overlap1 = gen_waveform_and_overlap(mid)
            print(
                "Obtained",
                np.exp((np.log(low) + np.log(high)) / 2),
                overlap1,
            )
            file_pointer.write(f"{fstar}, {alpha}, {sign}, {mid}\n")
