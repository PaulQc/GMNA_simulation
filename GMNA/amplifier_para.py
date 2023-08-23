# Parameter for the amplifier and propagation calculation
# Several set have been established and are compiled here.
# Simulation routine must choose one
#
from scipy import constants
c_mks = constants.value('speed of light in vacuum')  # m/s\n
c_nmps = constants.value('speed of light in vacuum') * 1e9 / 1e12  # c in nm/ps


# Parameters established for Wise's group results with Fiber Nufern PM-YDF-5-130
amplifier_Wise1 = {
    'fiber': {
        'length': 5.0,                      # [m]
        'core_radius': 5e-6/2,              # [m]
        'pump_cladding_radius': 130e-6/2,   # [m]
        'core_cladding_ratio': 0.0385,      # More precisely doping radius / pump radius
        'core_na': 0.12,
        'mode_diameter': 6.5e-6,            # Modefield according to spec sheet [m]
        'yb_number_density': 1.5e26,        # [1/m^3]
        'UpStateLifetime': 1.e-3,           # Upper state lifetime [s]
        'CrossSection':{
            'abs': './util/Nufern_abs.txt',
            'emi': './util/Nufern_emi.txt'
            },
        'background_loss': 0.03,            # [1/m]
        'wl_nm': 1025,                      # Center WL of fiber (for propagation methode) [nm]
        'betas': (0.025, 3.3e-5, -2.8e-8),  # B2, B3, B4 [ps^n/m]
        'gamma': 0.0025,                    # Gamma [1/(W m)]
        'simple_raman': False,              # See PyNLO description
        'raman': True,                      # Enable Raman effect?
        'steep': True                       # Enable self steepening?
        },
    'pump': {
        'wl': 976e-9,                       # Amplifier pump wavelength [m]
        'initial_power_forward': 1.,        # Amplifier pump power [W]
        'initial_power_backward': None      # Amplifier pump power [W]
        },
    'signal': {
        'wl': 1025e-9,                      # Amplifier pulse signal wavelength [m]
        'EPP': 0.5e-9,                      # Signal energy per pulse [J]
        'frep': 5e6,                        # Signal repetition frequency [Hz]
        'FWHM': 2.0e-12,                    # Signal pulse FWHM [s]
        'T_window': 40e-12,                 # Time window of signal pulse grid [s]
        'NPTS': 2**12                       # Number of point of signal pulse grid
        },
    'spectrum': {
        'f_Hd_minmax': (c_mks/1200e-9,
                        c_mks/980e-9),      # Hard minimum and maximum frequency values [Hz]
        'roi_factor': 0.01,                 # Minimum relative amplitude defining Region Of Interest
        'roi_warning': [],                  # Log values exceeding hard limit for the frequency value
        'roi_f_idx': [],                    # Index, in the pulse grid, of the current ROI
        'roi_f_minmax': [],                 # Minimum and maximum frequency value of the current ROI [Hz]
         },
    'channel': {
        'nb_ch_minmax': (10, 200),          # Min & Max number of signal channel for the amplifier model
        'nb_ch': None,                      # Number of simulation channel for the current pulse
        'f_mks': [],                        # Channel frequency for the current pulse [Hz]
        'dF_mks': None,                     # Channel spectral width for the current pulse [Hz]
        'wl_mks': [],                       # Channel wavelength for the current pulse [m]
        'dWl_mks': [],                      # Channel wavelength width for the current pulse [m]
        }
    }


# Parameters established for Wise's group results with Fiber Nufern PM-YDF-30-400
amplifier_Wise2 = {
    'fiber': {
        'length': 2.5,                      # [m]
        'core_radius': 30e-6/2,             # [m]
        'pump_cladding_radius': 400e-6/2,   # [m]
        'core_cladding_ratio': 0.075,       # More precisely doping radius / pump radius
        'core_na': 0.06,
        'mode_diameter': 21.e-6,            # Mode field according to spec sheet [m]
        'yb_number_density': 0.6e26,        # [1/m^3]
        'UpStateLifetime': 1.e-3,           # Upper state lifetime [s]
        'CrossSection':{
            'abs': './util/Nufern_abs.txt',
            'emi': './util/Nufern_emi.txt'
            },
        'background_loss': 0.03,            # [1/m]
        'wl_nm': 1025,                      # Center WL of fiber (for propagation methode) [nm]
        'betas': (0.019, 4.1e-5, -5.1e-8),  # B2, B3, B4 [ps^n/m]
        'gamma': 0.0003,                    # Gamma [1/(W m)]
        'simple_raman': False,              # See PyNLO description
        'raman': True,                      # Enable Raman effect?
        'steep': True                       # Enable self steepening?
        },
    'pump': {
        'wl': 976e-9,                       # Amplifier pump wavelength [m]
        'initial_power_forward': 15,        # Amplifier pump power [W]
        'initial_power_backward': None      # Amplifier pump power [W]
        },
    'signal': {
        'wl': 1025e-9,                      # Amplifier pulse signal wavelength [m]
        'EPP': 1.0e-9,                      # Signal energy per pulse [J]
        'frep': 5.4e6,                      # Signal repetition frequency [Hz]
        'FWHM': 0.6e-12,                    # Signal pulse FWHM [s]
        'T_window': 30e-12,                 # Time window of signal pulse grid [s]
        'NPTS': 2**12                       # Number of point of signal pulse grid
        },
    'spectrum': {
        'f_Hd_minmax': (c_mks/1250e-9,
                        c_mks/950e-9),      # Hard minimum and maximum frequency values [Hz]
        'roi_factor': 0.01,                 # Minimum relative amplitude defining Region Of Interest
        'roi_warning': [],                  # Log values exceeding hard limit for the frequency value
        'roi_f_idx': [],                    # Index, in the pulse grid, of the current ROI
        'roi_f_minmax': [],                 # Minimum and maximum frequency value of the current ROI [Hz]
         },
    'channel': {
        'nb_ch_minmax': (10, 400),          # Min & Max number of signal channel for the amplifier model
        'nb_ch': None,                      # Number of simulation channel for the current pulse
        'f_mks': [],                        # Channel frequency for the current pulse [Hz]
        'dF_mks': None,                     # Channel spectral width for the current pulse [Hz]
        'wl_mks': [],                       # Channel wavelength for the current pulse [m]
        'dWl_mks': [],                      # Channel wavelength width for the current pulse [m]
        }
    }

