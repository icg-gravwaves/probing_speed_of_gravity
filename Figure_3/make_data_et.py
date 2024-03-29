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

hp1, _ = get_dopplered_tf2_waveform_lisa(
    mass1=1.36,
    mass2=1.36,
    spin1z=0.0,
    spin2z=0.0,
    f_lower=5,
    delta_f=1.0 / 16384.0,
    delta_t=1.0 / 2048.0,
    distance=40,
    non_gr_cgwdelta_0=0,
    non_gr_transition_freq=1,
    non_gr_alpha=1,
    return_hc=False,
)
psdl = pycbc.psd.EinsteinTelescopeP1600143(len(hp1), hp1.delta_f, 1)
snr = pycbc.filter.sigma(hp1, psd=psdl, low_frequency_cutoff=5)
dist_overlap = 1 - 1 / (2 * snr)

print("Starting. Distinguishable overlap is", dist_overlap)

for fdx, fstar in enumerate([fstar_values[filenum]]):
    for adx, alpha in enumerate(alpha_values):
        for sign in [-1, 1]:

            def gen_waveform_and_overlap(curr_cgwdelta):
                hp2, _ = get_dopplered_tf2_waveform_lisa(
                    mass1=1.36,
                    mass2=1.36,
                    spin1z=0.0,
                    spin2z=0.0,
                    f_lower=5,
                    delta_f=1.0 / 16384.0,
                    delta_t=1.0 / 2048.0,
                    distance=40,
                    non_gr_cgwdelta_0=curr_cgwdelta_0,
                    non_gr_transition_freq=fstar,
                    non_gr_alpha=alpha,
                    return_hc=False,
                )
                overlap1, _ = pycbc.filter.match(
                    hp1, hp2, psd=psdl, low_frequency_cutoff=5
                )
                return overlap1

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
                    file_pointer.write(f"{fstar}, {alpha}, {curr_cgwdelta_0}\n")
                    continue
                low, high = curr_cgwdelta_0, curr_cgwdelta_0 * 10
            else:
                while overlap1 > dist_overlap:
                    curr_cgwdelta_0 = curr_cgwdelta_0 * 10
                    overlap1 = gen_waveform_and_overlap(curr_cgwdelta_0)
                    if curr_cgwdelta_0 > 10:
                        break
                if curr_cgwdelta_0 > 10:
                    file_pointer.write(f"{fstar}, {alpha}, {curr_cgwdelta_0}\n")
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
            file_pointer.write(f"{fstar}, {alpha}, {mid}\n")
