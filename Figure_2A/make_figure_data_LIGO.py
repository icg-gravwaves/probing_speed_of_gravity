# There is a lot of repetition in this code, and it could definitely be
# cleaned up a lot. However, this is what was used. Please let me know if this
# needs any clarification!

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

MPC_S_SI = lal.PC_SI * 1e6 / lal.C_SI


def get_delta_t_lisa(f_array, **params):
    cgw_delta0 = params["non_gr_cgwdelta_0"]
    fstar = params["non_gr_transition_freq"]
    alpha = params["non_gr_alpha"]
    delta_cgw = cgw_delta0 + 0.5 * (cgw_delta0) * (
        np.tanh(alpha * np.log(f_array / fstar)) + 1
    )
    delta_t = delta_cgw * params["distance"] * MPC_S_SI
    delta_t -= delta_t[-1]
    return delta_t


def get_delta_t_ligo(f_array, **params):
    non_gr_fac_n = params["non_gr_fac_n"]
    non_gr_fac_a_at_30 = params["non_gr_fac_a_at_30"]
    non_gr_fac_a = non_gr_fac_a_at_30 * 30**non_gr_fac_n
    return non_gr_fac_a * 40 * MPC_S_SI * f_array ** (-non_gr_fac_n)


get_delta_t = get_delta_t_ligo


def get_dopplered_tf2_waveform_lisa(
    dtype=np.complex128,
    return_hc=True,
    is_sequence=False,
    do_roll_and_taper=True,
    check_limits=False,
    **params
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


def get_newtimes_tf2_waveform_lisa(
    dtype=np.complex128,
    return_hc=True,
    is_sequence=False,
    do_roll_and_taper=True,
    check_limits=False,
    **params
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
    f_t_basic = (
        256.0
        / 5.0
        * (mc * lal.MTSUN_SI) ** (5.0 / 3.0)
        * np.pi ** (8.0 / 3.0)
        * (-sts_shifted_red)
    ) ** (-3.0 / 8.0)
    new_times = sts_shifted_red + get_delta_t(f_t_basic, **params)
    new_times = np.concatenate([new_times, sts_shifted[sts_shifted >= 0]])
    return new_times


non_gr_n_values = np.linspace(1, 10, 79)

outputs = np.zeros([79])
outputs_neg = np.zeros([79])

hp1, _ = get_dopplered_tf2_waveform_lisa(
    mass1=1.36,
    mass2=1.36,
    spin1z=0.0,
    spin2z=0.0,
    f_lower=20,
    delta_f=1.0 / 256.0,
    delta_t=1.0 / 4096.0,
    distance=40,
    non_gr_fac_n=1,
    non_gr_fac_a_at_30=0,
    return_hc=False,
)
psdl = pycbc.psd.aLIGOZeroDetHighPower(len(hp1), hp1.delta_f, 10)
snr = 32.4
dist_overlap = 1 - 1 / (2 * (snr) ** 2)

print("Starting. Distinguishable overlap is", dist_overlap)

for ndx, nval in enumerate(non_gr_n_values):
    for sign in [-1, 1]:
        if sign == 1:
            curr_output = outputs
        else:
            curr_output = outputs_neg
        print("Considering", nval, sign)
        curr_a_at_30 = 1e-15
        hp2, _ = get_dopplered_tf2_waveform_lisa(
            mass1=1.36,
            mass2=1.36,
            spin1z=0.0,
            spin2z=0.0,
            f_lower=20,
            delta_f=1.0 / 256.0,
            delta_t=1.0 / 4096.0,
            distance=40,
            non_gr_fac_n=nval,
            non_gr_fac_a_at_30=curr_a_at_30 * sign,
            return_hc=False,
        )
        overlap1, _ = pycbc.filter.match(hp1, hp2, psd=psdl, low_frequency_cutoff=10)
        if overlap1 < dist_overlap:
            while overlap1 < dist_overlap:
                curr_a_at_30 = curr_a_at_30 / 10
                hp2, _ = get_dopplered_tf2_waveform_lisa(
                    mass1=1.36,
                    mass2=1.36,
                    spin1z=0.0,
                    spin2z=0.0,
                    f_lower=20,
                    delta_f=1.0 / 256.0,
                    delta_t=1.0 / 4096.0,
                    distance=40,
                    non_gr_fac_n=nval,
                    non_gr_fac_a_at_30=curr_a_at_30 * sign,
                    return_hc=False,
                )
                overlap1, _ = pycbc.filter.match(
                    hp1, hp2, psd=psdl, low_frequency_cutoff=10
                )
                if curr_a_at_30 < 1e-25:
                    break
            if curr_a_at_30 < 1e-25:
                curr_output[ndx] = curr_a_at_30
                continue
            low, high = curr_a_at_30, curr_a_at_30 * 10
        else:
            while overlap1 > dist_overlap:
                curr_a_at_30 = curr_a_at_30 * 10
                hp2, _ = get_dopplered_tf2_waveform_lisa(
                    mass1=1.36,
                    mass2=1.36,
                    spin1z=0.0,
                    spin2z=0.0,
                    f_lower=20,
                    delta_f=1.0 / 256.0,
                    delta_t=1.0 / 4096.0,
                    distance=40,
                    non_gr_fac_n=nval,
                    non_gr_fac_a_at_30=curr_a_at_30 * sign,
                    return_hc=False,
                )
                overlap1, _ = pycbc.filter.match(
                    hp1, hp2, psd=psdl, low_frequency_cutoff=10
                )
                if curr_a_at_30 > 10:
                    break
            if curr_a_at_30 > 10:
                curr_output[ndx] = curr_a_at_30
                continue
            low, high = curr_a_at_30 / 10, curr_a_at_30

        for i in range(8):
            mid = np.exp((np.log(low) + np.log(high)) / 2)
            hp2, _ = get_dopplered_tf2_waveform_lisa(
                mass1=1.36,
                mass2=1.36,
                spin1z=0.0,
                spin2z=0.0,
                f_lower=20,
                delta_f=1.0 / 256.0,
                delta_t=1.0 / 4096.0,
                distance=40,
                non_gr_fac_n=nval,
                non_gr_fac_a_at_30=sign * mid,
                return_hc=False,
            )
            overlap1, _ = pycbc.filter.match(
                hp1, hp2, psd=psdl, low_frequency_cutoff=10
            )
            if overlap1 < dist_overlap:
                high = mid
            else:
                low = mid
        mid = np.exp((np.log(low) + np.log(high)) / 2)
        hp2, _ = get_dopplered_tf2_waveform_lisa(
            mass1=1.36,
            mass2=1.36,
            spin1z=0.0,
            spin2z=0.0,
            f_lower=20,
            delta_f=1.0 / 256.0,
            delta_t=1.0 / 4096.0,
            distance=40,
            non_gr_fac_n=nval,
            non_gr_fac_a_at_30=sign * mid,
            return_hc=False,
        )
        overlap1, _ = pycbc.filter.match(hp1, hp2, psd=psdl, low_frequency_cutoff=10)
        print(
            "Obtained",
            np.exp((np.log(low) + np.log(high)) / 2),
            overlap1,
        )
        curr_output[ndx] = mid

print(outputs)
print(outputs_neg)

np.savez(
    "figure2_dataLIGO.npz",
    non_gr_n_values=non_gr_n_values,
    outputs=outputs,
    outputs_sl=outputs_neg * sign,
)
