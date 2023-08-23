# Paul Grenier 2021-02-26
#
# Module containing methode in support to the "Gain Managed Nonlinear Amplification" simulation
#

import matplotlib.pyplot as plt
import numpy as np
import pynlo
from pyfiberamp import helper_funcs as hf
from pyfiberamp.fibers import DoubleCladFiber
from pyfiberamp.spectroscopies import Spectroscopy
from pyfiberamp.steady_state import SteadyStateSimulation
from pyfiberamp.mode_solver import GaussianMode
from scipy import constants, interpolate

c_mks = constants.value('speed of light in vacuum')  # m/s\n
c_nmps = constants.value('speed of light in vacuum') * 1e9 / 1e12  # c in nm/ps


# TODO: method to experiment and implement
#  chirp_pulse_W() --> calculate achievable minimum duration
#
# TODO : within class fiber:
#  set_gamma_function() - Capability to have gamma to change as a function of z --> For taper fiber model
#


def cal_pulse_PowerSpectrum(pulse):
    """Function to calculate pulse's spectral power in [W per frequency bin]
    Return: A vector array of power values, same length as the pulse"""
    return np.absolute(pulse.AW) ** 2 / pulse.dF_mks * pulse.frep_mks  # [W] Average power in each bin


def cal_amp_roi(x, ratio_min=0.01):
    """
    Function to calculate region of interest (ROI), defined as the region where the absolute value of x
    is larger or equal to "ratio_min * maximum value of x"
    :param x: array of  complex
    :param ratio_min: float The minimum ratio of amplitude defining the ROI
    :returns: ROI - An boolean array same length as x
    """
    #
    roi = np.greater_equal(np.abs(x), np.max(np.abs(x)) * ratio_min)
    return roi


def cal_amp_phase(x, idx_0=None):
    """
    Function to calculate the phase of the complex vector x
    Can accept the index idx_0 where the phase is set to zero
    :param x: array of complex
    :param idx_0: Index where the phase is to be set = 0
    :returns: array of float with unwrapped phase in radian, relative to idx_0, if provided
    """
    #
    x_phase = np.unwrap(np.angle(x))
    if idx_0 is not None:
        x_phase = x_phase - x_phase[idx_0]
    return x_phase


# TODO: Anomaly when the pulse is highly distorted and the ROI is not continuous in frequency domain
def channel_roi_PowerSpectrum(pulse, amplifier):
    """Function that calculates parameters for the signal channel array based on the current pulse characteristics,
    that are stored in dictionary amplifier (roi, minmax, nb_of_channel, ...)
    Then it prepares the signal channel inputs array to be used in amplifier simulation model.
    Return: Vector arrays of laser signal power [W] in each channel"""

    # Select the ROI in the pulse spectrum
    pulse_PowerSpectrum = cal_pulse_PowerSpectrum(pulse)  # [W] Average power in each bin
    roi_f_idx = cal_amp_roi(pulse_PowerSpectrum, ratio_min=amplifier['spectrum']['roi_factor'])
    roi_f_minmax = [np.min(pulse.F_mks[roi_f_idx]), np.max(pulse.F_mks[roi_f_idx])]

    # Verify if pulse spectrum ROI is too wide
    warning = False
    if roi_f_minmax[0] < amplifier['spectrum']['f_Hd_minmax'][0]:
        warning = True
        amplifier['spectrum']['roi_warning'].append(roi_f_minmax[0] / 1e12)  # Log exceeding value
        roi_f_minmax[0] = amplifier['spectrum']['f_Hd_minmax'][0]
    if roi_f_minmax[1] > amplifier['spectrum']['f_Hd_minmax'][1]:
        warning = True
        amplifier['spectrum']['roi_warning'].append(roi_f_minmax[1] / 1e12)  # Log exceeding value
        roi_f_minmax[1] = amplifier['spectrum']['f_Hd_minmax'][1]
    if warning:
        roi_f_idx = (np.greater_equal(pulse.F_mks, roi_f_minmax[0]) &
                     np.less_equal(pulse.F_mks, roi_f_minmax[1]))
        print('*   *   *   *   *\n Pulse ROI spectrum exceeding fixed limit \n*   *   *   *   *')
    amplifier['spectrum']['roi_f_idx'] = roi_f_idx
    amplifier['spectrum']['roi_f_minmax'] = roi_f_minmax

    # Set number of spectral point/channel to be use in amplifier calculation
    nb_channel_factor = amplifier['channel']['nb_ch_minmax'][1] / (amplifier['spectrum']['f_Hd_minmax'][1] -
                                                                   amplifier['spectrum']['f_Hd_minmax'][0])
    nb_of_channel = int(np.ceil(nb_channel_factor * (roi_f_minmax[1] - roi_f_minmax[0])))
    if nb_of_channel < amplifier['channel']['nb_ch_minmax'][0]:
        nb_of_channel = amplifier['channel']['nb_ch_minmax'][0]
    if nb_of_channel > amplifier['channel']['nb_ch_minmax'][1]:
        nb_of_channel = amplifier['channel']['nb_ch_minmax'][1]
    amplifier['channel']['nb_ch'] = nb_of_channel

    # Set other arrays and values for the defined number of signal channels 
    amplifier['channel']['f_mks'] = np.linspace(roi_f_minmax[0], roi_f_minmax[1], num=nb_of_channel)
    amplifier['channel']['dF_mks'] = amplifier['channel']['f_mks'][1] - amplifier['channel']['f_mks'][0]
    amplifier['channel']['wl_mks'] = 299792458. / amplifier['channel']['f_mks']
    amplifier['channel']['dWl_mks'] = amplifier['channel']['wl_mks'] / amplifier['channel']['f_mks'] * \
                                      amplifier['channel']['dF_mks']

    # Set the channel input power for the amplifier simulation
    # f_PowerSpectrum = interpolate.interp1d(pulse.F_mks[roi_f_idx],
    #                                       pulse_PowerSpectrum[roi_f_idx], kind='cubic')
    nb_average = int(amplifier['channel']['dF_mks'] // pulse.dF_mks + 1)
    pulse_f_averaged = moving_average(pulse.F_mks, nb_average)
    pulse_PowerSpectrum_averaged = moving_average(pulse_PowerSpectrum, nb_average)
    func_PowerSpectrum = interpolate.interp1d(pulse_f_averaged, pulse_PowerSpectrum_averaged, kind='cubic')
    channel_PowerSpectrum = func_PowerSpectrum(amplifier['channel']['f_mks']) * \
                            amplifier['channel']['dF_mks'] / pulse.dF_mks  # To be in 'W per channel' unit

    return channel_PowerSpectrum


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


def my_pulse(amplifier):
    """Function to create pulse instance using parameter in the dictionary amplifier
    Number de points and Time window need to be carefully selected for each case.
    Return : A pulse instance"""
    #
    # Pulse instance and general parameters of the pulse
    pulse = pynlo.light.PulseBase.Pulse()
    pulse.set_center_wavelength_m(amplifier['signal']['wl'])
    pulse.set_time_window_s(amplifier['signal']['T_window'])
    pulse.set_NPTS(amplifier['signal']['NPTS'])
    # print('f min = {:0.2f}, f max = {:0.2f}, dF = {:0.5f} THz'.format(pulse_in.F_THz[0],
    #                                                                  pulse_in.F_THz[-1], pulse_in.dF_THz))
    # Gaussian pulse
    pulse.set_AT(np.exp(-2.77 * 0.5 * pulse.T_mks ** 2 / (amplifier['signal']['FWHM'] ** 2)))
    # Energy per pulse, rep rate and Power spectrum
    pulse.set_epp(amplifier['signal']['EPP'])
    pulse.set_frep_MHz(amplifier['signal']['frep'] / 1e6)
    return pulse


def run_amplifier_sim(amplifier):
    """Run a simulation of the amplifier power evolution only,
    i.e. without nonlinear pulse propagation calculation.
    This to have a first look at the signal and pump power evolution along the active fiber.
    Generate plot of pump and signal evolution along the active fiber"""

    # Load proper cross section data
    yb_spectroscopy = Spectroscopy.from_files(
        absorption_cross_section_file=amplifier['fiber']['CrossSection']['abs'],
        emission_cross_section_file=amplifier['fiber']['CrossSection']['emi'],
        upper_state_lifetime=amplifier['fiber']['UpStateLifetime'],
        interpolate='spline')  # alternatively: interpolate='linear'

    # Create the fiber
    yb_dc_fiber = DoubleCladFiber(spectroscopy=yb_spectroscopy,
                                  ion_number_density=amplifier['fiber']['yb_number_density'],
                                  length=amplifier['fiber']['length'],
                                  core_radius=amplifier['fiber']['core_radius'],
                                  core_na=amplifier['fiber']['core_na'],
                                  background_loss=amplifier['fiber']['background_loss'],
                                  ratio_of_core_and_cladding_diameters=amplifier['fiber']['core_cladding_ratio'])

    # yb_dc_fiber.default_signal_mode_shape_parameters['functional_form'] = 'gaussian'

    # Create the simulation instance
    simulation = SteadyStateSimulation(fiber=yb_dc_fiber)
    simulation.fiber = yb_dc_fiber
    simulation.add_forward_pump(amplifier['pump']['wl'], amplifier['pump']['initial_power_forward'])
    #
    signal_power = amplifier['signal']['EPP'] * amplifier['signal']['frep']
    # Approximate the signal BW to half of the calculated wavelentgth channel
    nb_channel = amplifier['channel']['nb_ch']
    signal_BW = amplifier['channel']['dWl_mks'][int(nb_channel / 2)] * nb_channel / 2
    simulation.add_forward_signal(wl=amplifier['signal']['wl'],
                                  input_power=signal_power,
                                  wl_bandwidth=signal_BW,
                                  mode=GaussianMode(mfd=amplifier['fiber']['mode_diameter'],
                                                    core_radius=amplifier['fiber']['core_radius']))
    # Modefield according to spec sheet
    # Run the simulation
    result = simulation.run(tol=1e-4)

    result.plot()


def run_amplifier_and_propagation_sim(pulse, amplifier, nb_section, add_noise=False, noise_factor=10):
    """
    Function that run the full simulation of amplification + propagation of the pulse
    The active fiber is divided in nb_section and the amplification + propagation are
    applied successively on each section. Option to add noise to the pulse at each section.
    :param pulse: Pulse instance at the entrance of fiber
    :param amplifier: Dictionary for the amplifier parameters
    :param nb_section: integer Number of section into which fiber is divided
    :param add_noise: boolean If noise is to be added to the pulse during propagation
    :param noise_factor: integer Number of time noise is added, this after each segment
    :returns: 4 arrays of float The evolution at each segment of the Pump power, Signal power,
    Pulse spectral amplitude, and Pulse temporal amplitude, and the last "pulse_out" of the
    PyNLO propagation method calculation
    """
    #
    # Array to store parameter evolution along the fiber
    pump_power_evo = np.zeros(nb_section + 1)
    signal_power_evo = np.zeros(nb_section + 1)
    AW_evo = np.zeros((pulse.NPTS, nb_section + 1), dtype=np.complex64)
    AT_evo = np.zeros((pulse.NPTS, nb_section + 1), dtype=np.complex64)

    # Value at fiber entrance
    channel_PowerSpectrum = channel_roi_PowerSpectrum(pulse, amplifier)
    pump_power_evo[0] = amplifier['pump']['initial_power_forward']
    signal_power_evo[0] = np.sum(channel_PowerSpectrum)
    AW_evo[:, 0] = pulse.AW
    AT_evo[:, 0] = pulse.AT

    section_length = amplifier['fiber']['length'] / nb_section

    # Create the fiber-section for the amplifier similation
    yb_spectroscopy = Spectroscopy.from_files(
        absorption_cross_section_file=amplifier['fiber']['CrossSection']['abs'],
        emission_cross_section_file=amplifier['fiber']['CrossSection']['emi'],
        upper_state_lifetime=amplifier['fiber']['UpStateLifetime'],
        interpolate='spline')  # alternatively: interpolate='linear'
    yb_dc_fiber_section = DoubleCladFiber(spectroscopy=yb_spectroscopy,
                                          ion_number_density=amplifier['fiber']['yb_number_density'],
                                          length=section_length,
                                          core_radius=amplifier['fiber']['core_radius'],
                                          core_na=amplifier['fiber']['core_na'],
                                          background_loss=amplifier['fiber']['background_loss'],
                                          ratio_of_core_and_cladding_diameters=amplifier['fiber'][
                                              'core_cladding_ratio'])
    # yb_dc_fiber_section.default_signal_mode_shape_parameters['functional_form'] = 'gaussian'

    Steps = 2
    #
    for section_no in range(1, nb_section + 1):
        # Create the simulation instance
        simulation = SteadyStateSimulation(fiber=yb_dc_fiber_section)
        simulation.fiber = yb_dc_fiber_section
        simulation.add_forward_pump(amplifier['pump']['wl'], pump_power_evo[section_no - 1])
        # simulation.add_ase(wl_start=1000e-9, wl_end=1100e-9, n_bins=50)
        # for loop in reverse order to have channel in increasing order of wavelength
        for i in range(-1, -1 * (amplifier['channel']['nb_ch'] + 1), -1):
            simulation.add_forward_signal(wl=amplifier['channel']['wl_mks'][i],
                                          input_power=channel_PowerSpectrum[i],
                                          wl_bandwidth=amplifier['channel']['dWl_mks'][i],
                                          mode=GaussianMode(mfd=amplifier['fiber']['mode_diameter'],
                                                            core_radius=amplifier['fiber']['core_radius'])
                                          )  # Modefield according to spec sheet

        # Run the simulation
        result = simulation.run(tol=1e-4)

        # Get the amplifier simulation results
        channel_slice = result.channels.get_slices()
        # Flip, to get back in increasing order of frequency
        channel_output = np.flip(result.powers[channel_slice['forward_signal'], -1])
        channel_gain = channel_output / channel_PowerSpectrum

        # Calculate the spectral gain for the pulse
        spectral_gain = np.ones_like(pulse.F_mks)
        f_gain = interpolate.interp1d(amplifier['channel']['f_mks'], channel_gain, kind='cubic')
        spectral_gain[amplifier['spectrum']['roi_f_idx']] = f_gain(pulse.F_mks[amplifier['spectrum']['roi_f_idx']])
        spectral_gain = np.log(spectral_gain) / section_length

        # Create the fiber for the nonlinear propagation with the calculated spectral gain
        nonlinear_fiber = pynlo.media.fibers.fiber.FiberInstance()
        nonlinear_fiber.generate_fiber_PG(section_length, center_wl_nm=amplifier['fiber']['wl_nm'],
                                          betas=amplifier['fiber']['betas'], gamma_W_m=amplifier['fiber']['gamma'],
                                          gain=spectral_gain)
        # Run the nonlinear propagation simulation
        if add_noise:
            for j in range(noise_factor):
                pulse.add_noise()
        evol = pynlo.interactions.FourWaveMixing.SSFM.SSFM(local_error=0.005,
                                                           USE_SIMPLE_RAMAN=amplifier['fiber']['simple_raman'],
                                                           disable_Raman=np.logical_not(amplifier['fiber']['raman']),
                                                           disable_self_steepening=np.logical_not(
                                                               amplifier['fiber']['steep']))
        y, AW, AT, pulse_out = evol.propagate(pulse_in=pulse, fiber=nonlinear_fiber, n_steps=Steps)

        # New input pulse is the output pulse from previous iteration
        pulse = pulse_out.create_cloned_pulse()
        # New pulse ROI and channel input power for amplifier model
        channel_PowerSpectrum = channel_roi_PowerSpectrum(pulse, amplifier)

        # Store this iteration results
        pump_power_evo[section_no] = result.powers[channel_slice['forward_pump']][0, :][-1]
        signal_power_evo[section_no] = np.sum(channel_PowerSpectrum)
        AW_evo[:, section_no] = AW[:, -1]
        AT_evo[:, section_no] = AT[:, -1]
    return pump_power_evo, signal_power_evo, AW_evo, AT_evo, pulse_out
