# Paul Grenier 2021-03-03
"""Regroup plotting functions used in the Gain Managed Nonlinear Amplification simulation"""
import numpy as np
from matplotlib import pyplot as plt
from util.GMNA_Sim import cal_amp_roi, cal_amp_phase


def plot_linear_AW_AT(pulse, AW, AT, nb_curve=4):
    """
    Function that produce linear plots of the spectral (AW) and temporal (AT) squared amplitude.
    It plots the amplitudes calculated at "nb_curve" locations from beginning to end of fiber.
    :param pulse: one of pulse instance used during the simulation, for the frequency
     and time grid definition
    :param AW: numpy float array" Set of pulse spectral amplitude at the end of each fiber segments
    :param AT: numpy float array: Set of pulse temporal amplitude at the end of each fiber
    :param nb_curve: int: Desired number of curve to be plot
    :return: None
    """
    #
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    fig, [ax0, ax1] = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ax0_right = ax0.twinx()
    ax1_right = ax1.twinx()
    nb_data = np.shape(AW)[1]
    selected_data = np.linspace(0, nb_data, num=nb_curve, dtype=int)
    selected_data[-1] = nb_data -1   # Make sure to have the last curve at output end
    for i, no in enumerate(selected_data):
        roi = cal_amp_roi(AW[:, no], ratio_min=0.1)
        ax0.plot(pulse.wl_nm[roi], np.abs(AW[roi, no]) ** 2 / np.max(np.abs(AW[roi, no]) ** 2),
                 color=colors[i])
        ax0_right.plot(pulse.wl_nm[roi], cal_amp_phase(AW[roi, no], idx_0=np.sum(roi) // 2),
                       color=colors[i], linestyle=':')
        #
        roi = cal_amp_roi(AT[:, no], ratio_min=0.1)
        ax1.plot(pulse.T_ps[roi], np.abs(AT[roi, no]) ** 2 / np.max(np.abs(AT[roi, no]) ** 2),
                 color=colors[i])
        ax1_right.plot(pulse.T_ps[roi], cal_amp_phase(AT[roi, no], idx_0=np.sum(roi) // 2),
                       color=colors[i], linestyle=':')
    ax0.set_ylabel('Power (a.u.)')
    ax0.set_xlabel('Wavelength (nm)')
    ax1.set_xlabel('Time (ps)')
    ax1_right.set_ylabel('Phase (rad)')
    plt.show()

def plot_log_AW_AT(pulse, AW, AT, nb_curve=4):
    """
    Function that produce log plots of the spectral (AW) and temporal (AT) squared amplitude.
    It plots the amplitudes calculated at "nb_curve" locations from beginning to end of fiber.
    :param pulse: one of pulse instance used during the simulation, for the frequency
     and time grid definition
    :param AW: numpy float array" Set of pulse spectral amplitude at the end of each fiber segments
    :param AT: numpy float array: Set of pulse temporal amplitude at the end of each fiber
    :param nb_curve: int: Desired number of curve to be plot
    :return: None
    """
    #
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    fig, [ax0, ax1] = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ax0_right = ax0.twinx()
    ax1_right = ax1.twinx()
    nb_data = np.shape(AW)[1]
    selected_data = np.linspace(0, nb_data, num=nb_curve, dtype=int)
    selected_data[-1] = nb_data - 1  # Make sure to have the last curve at output end
    for i, no in enumerate(selected_data):
        roi = cal_amp_roi(AW[:, no], ratio_min=0.1)
        y = np.log10(np.abs(AW[roi, no]) ** 2 / np.max(np.abs(AW[roi, no]) ** 2)) * 10
        ax0.plot(pulse.wl_nm[roi], y, color=colors[i])
        ax0_right.plot(pulse.wl_nm[roi], cal_amp_phase(AW[roi, no], idx_0=np.sum(roi) // 2),
                       color=colors[i], linestyle=':')
        #
        roi = cal_amp_roi(AT[:, no], ratio_min=0.1)
        y = np.log10(np.abs(AT[roi, no]) ** 2 / np.max(np.abs(AT[roi, no]) ** 2)) * 10
        ax1.plot(pulse.T_ps[roi], y, color=colors[i])
        ax1_right.plot(pulse.T_ps[roi], cal_amp_phase(AT[roi, no], idx_0=np.sum(roi) // 2),
                       color=colors[i], linestyle=':')
    ax0.set_ylabel('Normalized power dB')
    ax0.set_xlabel('Wavelength (nm)')
    ax1.set_xlabel('Time (ps)')
    ax1_right.set_ylabel('Phase (rad)')
    plt.show()
    return


def plot_AW_AT_vs_length(pulse, AW, AT, length, ratio_min=[0.001, 0.02]):
    """
    Function that produce color map log plots of the spectral (AW) and temporal (AT)
    amplitude vs position in the fiber. "ratio_min" are the minimum amplitude to define
    the ROI for the AW and AT plot
    :param pulse: one of pulse instance used during the simulation, for the frequency
     and time grid definition
    :param AW: numpy float array" Set of pulse spectral amplitude at the end of each fiber segments
    :param AT: numpy float array: Set of pulse temporal amplitude at the end of each fiber
    :param length: float: Fiber length
    :param ratio_min: list of 2 float: Ratio minimum for the AW and AT
    :return: None
    """
    # Spectral amplitude
    fig, [ax0, ax1] = plt.subplots(2, 1, figsize=(10, 5), constrained_layout=True)
    roi = cal_amp_roi(AW[:, -1], ratio_min=ratio_min[0])
    zW = np.log10(np.abs(AW[roi]) ** 2) * 10
    extent = [0, length, pulse.wl_nm[roi][-1], pulse.wl_nm[roi][0]]
    im = ax0.imshow(zW, aspect='auto', extent=extent, vmin=np.max(zW) - 30.0, vmax=np.max(zW), cmap='jet')
    ax0.set_ylabel('Wavelength (nm)')
    fig.colorbar(im, ax=ax0, label='Normaliser Spectral power (dB)')
    # Temporal amplitude
    roi = cal_amp_roi(AT[:, -1], ratio_min=ratio_min[1])
    zA = np.flip(np.log10(np.abs(AT[roi]) ** 2), axis=0) * 10
    extent = [0, length, pulse.T_ps[roi][0], pulse.T_ps[roi][-1]]
    im = ax1.imshow(zA, aspect='auto', extent=extent, vmin=np.max(zA) - 30.0, vmax=np.max(zA), cmap='jet')
    ax1.set_xlabel('Fiber position (m)')
    ax1.set_ylabel('Time (ps)')
    fig.colorbar(im, ax=ax1, label='Peak power dBW')
    plt.show()
    return


def plot_PeakPower_EPP_vs_length(pulse, AT, length):
    """Function that produce plots of the peak power and pulse energy
    amplitude vs position in the fiber."""
    #
    fig, ax0 = plt.subplots(1, 1, figsize=(6, 5), constrained_layout=True)
    ax0_right = ax0.twinx()
    nb_section = np.shape(AT)[1]
    x = np.linspace(0, length, num=nb_section)
    y = np.max(np.abs(AT)**2, axis=0)/1000    # Peak Power in kW
    ax0.plot(x, y, color = 'b', label='Peak power')
    #
    y = pulse.dT_mks * np.trapz(abs(AT)**2, axis=0) * 1e9   # Energy in nJ
    ax0_right.plot(x, y, color = 'r', label='Energy')
    ax0.set_ylabel('Peak power (kW)', color='b')
    ax0.set_xlabel('Fiber position (m)')
    ax0_right.set_ylabel('Energy (nJ)', color='r')
    plt.show()
    return


def finalize_my_plot(ax, x_label, y_label, xlim=None):
    """Function to set few aspects to a plot for a nicer look: inner tick, thicker line, ..."""
    #
    ax.tick_params(which='minor', direction='in', left=True, right=False, top=True, bottom=False,
                   width=1, length=3)
    ax.tick_params(which='major', direction='in', left=True, right=False, top=True, bottom=True,
                   width=2, length=5, labelsize=12)
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    if xlim is not None:
        ax.set_xlim(xlim)
    return
