"""Gadolinium concentration curves of blood, myocardium and infarct (Fermi model).

Tissue concentration is modelled as the convolution of an arterial input
function (AIF) with a Fermi-shaped impulse residue function (Jerosch-Herold
et al., Med Phys 1998)::

    C_tissue(t) = MBF * (C_art (*) R)(t - shift),
    R(t) = (1 + beta) / (1 + beta * exp(alpha * t)),

where ``MBF`` is the myocardial blood flow [ml/g/s], ``alpha`` [s^-1] sets the
wash-out rate and ``beta`` [-] the duration of the plateau (``R(0) = 1``), and
``shift`` [s] is the transit delay between the left ventricle and the tissue.
The AIF is the MRXCAT first-pass curve (gamma-variate, no recirculation),
rescaled to a 12 mM peak for 0.1 mmol/kg and linearly with dose.

Command line: ``t1mes_simulate_fermi --help``.
"""
import argparse

import numpy as np
from scipy.interpolate import PchipInterpolator


######################
# CONSTANTS
######################

# Arterial input function [a.u.] sampled every second, measured with a half
# dose (0.05 mmol/kg) on 6 volunteers and fitted with a gamma-variate for the
# first pass only (from MRXCAT, Wissmann et al., JCMR 2014).
AIF005 = np.array(
    [0, 0, 0, 0, 0.0009, 0.0178, 0.1343, 0.5466, 1.4530, 2.8276,
     4.3365, 5.5112, 6.0165, 5.7934, 5.0205, 3.9772, 2.9161, 1.9988,
     1.2913, 0.7916, 0.4632, 0.2598, 0.1404, 0.0733, 0.0371, 0.0182,
     0.0087, 0.0041, 0.0019, 0.0008, 0.0004, 0.0002, 0.0001, 0.0000,
     0.0000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
     ]
)
AIF_TIME_STEP = 1.0  # [s]

# Peak arterial concentration for a 0.1 mmol/kg dose [mM] (see get_concentration_fermi)
AIF_PEAK_01 = 12.0

MBF_REST = 1.0 / 60    # myocardial blood flow at rest [ml/g/s] (1.0 ml/g/min)
MBF_STRESS = 3.5 / 60  # myocardial blood flow under stress [ml/g/s] (3.5 ml/g/min)


######################
# MODEL FUNCTIONS
######################

def fermi_function(t, alpha, beta):
    """Fermi impulse residue function ``R(t) = (1 + beta) / (1 + beta * exp(alpha * t))``.

    Args:
        t (np.ndarray): time [s].
        alpha (float): wash-out rate [s^-1].
        beta (float): plateau parameter [-].

    Returns:
        np.ndarray: impulse residue function (``R(0) = 1``).
    """
    return (1 + beta) / (1 + beta * np.exp(alpha * t))


def convolve_concentration_impulse_residue_function(
    t_upsample,
    arterial_concentration,
    impulse_residue_function,
    myocardial_blood_flow,
):
    """Tissue concentration ``MBF * (C_art (*) R)(t)`` (causal discrete convolution).

    Args:
        t_upsample (np.ndarray): uniformly sampled time [s].
        arterial_concentration (np.ndarray): arterial concentration on ``t_upsample`` [mM].
        impulse_residue_function (np.ndarray): impulse residue function on ``t_upsample``.
        myocardial_blood_flow (float): myocardial blood flow [ml/g/s].

    Returns:
        np.ndarray: tissue concentration on ``t_upsample`` [mM].
    """
    dt = (t_upsample[-1] - t_upsample[0]) / (len(t_upsample) - 1)
    n = len(arterial_concentration)
    conv = np.convolve(arterial_concentration, impulse_residue_function)[:n]
    return myocardial_blood_flow * conv * dt


def shift_concentration(t, concentration, t_shift):
    """Delay a concentration curve by ``t_shift`` [s] (monotone PCHIP interpolation).

    Args:
        t (np.ndarray): time [s].
        concentration (np.ndarray): concentration on ``t`` [mM].
        t_shift (float): positive delay [s].

    Returns:
        np.ndarray: ``concentration(t - t_shift)`` sampled on ``t``.
    """
    interpolator = PchipInterpolator(t + t_shift, concentration, extrapolate=True)
    return interpolator(t)


def get_initial_concentration(T_post, T_native, relaxivity, time_unit='ms'):
    """Concentration explaining a relaxation time change (inverse of the relaxivity model)::

        C = (1/T_post - 1/T_native) / relaxivity

    Args:
        T_post (float): relaxation time after injection, in ``time_unit``.
        T_native (float): native relaxation time, in ``time_unit``.
        relaxivity (float): relaxivity of the contrast agent [s^-1 mM^-1].
        time_unit (str, optional): 'ms' (default) or 's'.

    Returns:
        float: contrast agent concentration [mM].
    """
    assert T_post < T_native, "You must have T_post < T_native"
    assert time_unit in ['ms', 's'], "time_unit should be either 'ms' or 's'"
    time_factor = 1e-3 if time_unit == 'ms' else 1
    return 1 / relaxivity * (1 / (T_post * time_factor) - 1 / (T_native * time_factor))


######################
# MAIN FUNCTIONS
######################

def get_concentration_fermi(times_post, alpha, beta, contrast_dose=0.075, is_stress=False,
                            time_step=0.2, lv_myo_shift=3):
    """Arterial and tissue gadolinium concentration with the Fermi model.

    Args:
        times_post (np.ndarray): time after injection of each output sample [s].
        alpha (float): Fermi wash-out rate [s^-1].
        beta (float): Fermi plateau parameter [-].
        contrast_dose (float, optional): injected dose [mmol/kg]. Defaults to 0.075.
        is_stress (bool, optional): use the stress (True) or rest (False) blood flow. Defaults to False.
        time_step (float): time step of the internal convolution grid [s]. Defaults to 0.2.
        lv_myo_shift (float): delay between left ventricle and tissue [s]. Defaults to 3.

    Returns:
        np.ndarray: arterial concentration at ``times_post`` [mM].
        np.ndarray: tissue concentration at ``times_post`` [mM].

    Peak arterial concentration (Gadovist, 1.0 mmol/ml), rough estimate:
        - 75 kg x 0.1 mmol/kg -> 7.5 ml of contrast agent injected;
        - stroke volume ~80-100 ml with ~80 % ejection fraction -> ~75 ml per beat;
        - normalised AIF peak ~12 % of the bolus per beat -> 0.9 ml in 75 ml of blood,
          i.e. a dilution of 0.012;
        - peak concentration ~1000 mM x 0.012 = 12 mM.
    """
    t = np.arange(0, len(AIF005)) * AIF_TIME_STEP
    myocardial_blood_flow = MBF_STRESS if is_stress else MBF_REST

    # Arterial concentration scaled to the requested dose
    aif01 = AIF005 * AIF_PEAK_01 / np.max(AIF005)
    arterial_concentration = aif01 * contrast_dose / 0.1

    # Resample on a fine grid; the AIF is 0 after its last sample (first pass only)
    t_upsample = np.arange(0, times_post[-1], time_step)
    arterial_concentration_upsample = np.interp(t_upsample, t, arterial_concentration)

    impulse_residue_function = fermi_function(t_upsample, alpha, beta)
    myocardial_concentration_upsample = convolve_concentration_impulse_residue_function(
        t_upsample,
        arterial_concentration_upsample,
        impulse_residue_function,
        myocardial_blood_flow,
    )
    myocardial_concentration_upsample = shift_concentration(
        t_upsample, myocardial_concentration_upsample, lv_myo_shift)

    # Sample at the requested times
    arterial_concentration = np.interp(times_post, t_upsample, arterial_concentration_upsample)
    myocardial_concentration = np.interp(times_post, t_upsample, myocardial_concentration_upsample)
    return arterial_concentration, myocardial_concentration


def get_concentrations(
    times_post,
    alpha_myo,
    beta_myo,
    alpha_inf,
    beta_inf,
    lv_myo_shift=3,
    lv_inf_shift=3,
    contrast_dose=0.075,
    is_stress=False,
    time_step=.2,
    plot=True,
):
    """Concentration of blood, myocardium and infarct over time.

    Args:
        times_post (np.ndarray): time after injection [s].
        alpha_myo (float): Fermi alpha of the myocardium [s^-1].
        beta_myo (float): Fermi beta of the myocardium [-].
        alpha_inf (float): Fermi alpha of the infarct [s^-1].
        beta_inf (float): Fermi beta of the infarct [-].
        lv_myo_shift (float): delay between left ventricle and myocardium [s]. Defaults to 3.
        lv_inf_shift (float): delay between left ventricle and infarct [s]. Defaults to 3.
        contrast_dose (float, optional): injected dose [mmol/kg]. Defaults to 0.075.
        is_stress (bool, optional): stress (True) or rest (False) blood flow. Defaults to False.
        time_step (float, optional): time step of the internal convolution grid [s]. Defaults to 0.2.
        plot (bool, optional): plot the curves with matplotlib. Defaults to True.

    Returns:
        np.ndarray: arterial (blood) concentration [mM].
        np.ndarray: myocardium concentration [mM].
        np.ndarray: infarct concentration [mM].
    """
    _, c_myo = get_concentration_fermi(
        times_post, alpha_myo, beta_myo, contrast_dose,
        is_stress, time_step, lv_myo_shift,
    )
    c_art, c_inf = get_concentration_fermi(
        times_post, alpha_inf, beta_inf, contrast_dose,
        is_stress, time_step, lv_inf_shift,
    )

    if plot:
        import matplotlib.pyplot as plt
        plt.plot(times_post / 60, c_art)
        plt.plot(times_post / 60, c_myo)
        plt.plot(times_post / 60, c_inf)
        plt.legend(['Blood', 'Myocardium', 'Infarction'])
        plt.xlabel('Time post-injection (minutes)', fontsize='large')
        plt.ylabel('Gadolinium concentration (mM)', fontsize='large')
        plt.show()
    return c_art, c_myo, c_inf


def main():
    parser = argparse.ArgumentParser(
        description="Compute gadolinium concentrations of blood, myocardium and infarct with the Fermi model.")
    parser.add_argument('-t_start', '--time_start', type=float, default=0.,
                        help='Start time after injection [s]. Default: 0')
    parser.add_argument('-d', '--duration', type=float, default=60. * 5,
                        help='Duration of the sampled window [s]. Default: 300 (5 minutes)')
    parser.add_argument('-s', '--nb_sampled', type=int, default=300,
                        help='Number of time points in the window. Default: 300')
    parser.add_argument('--alpha_myo', type=float, default=.01, help='Fermi alpha of the myocardium [s^-1]. Default: 0.01')
    parser.add_argument('--beta_myo', type=float, default=.01, help='Fermi beta of the myocardium. Default: 0.01')
    parser.add_argument('--alpha_inf', type=float, default=.01, help='Fermi alpha of the infarct [s^-1]. Default: 0.01')
    parser.add_argument('--beta_inf', type=float, default=.005, help='Fermi beta of the infarct. Default: 0.005')
    parser.add_argument('--lv_myo_shift', type=float, default=3, help='Delay between LV and myocardium [s]. Default: 3')
    parser.add_argument('--lv_inf_shift', type=float, default=3, help='Delay between LV and infarct [s]. Default: 3')
    parser.add_argument('--contrast_dose', type=float, default=0.1, help='Contrast dose [mmol/kg]. Default: 0.1')
    parser.add_argument('--is_stress', action=argparse.BooleanOptionalAction, default=True,
                        help='Stress (--is_stress) or rest (--no-is_stress) myocardial blood flow')
    parser.add_argument('--time_step', type=float, default=0.2,
                        help='Time step of the internal convolution grid [s]. Default: 0.2')
    parser.add_argument('--plot', action=argparse.BooleanOptionalAction, default=True,
                        help='Plot the concentration curves')
    parser.add_argument('-o', '--output_file', type=str, default='concentration_fermi_default.npz',
                        help='Output file (.npz). Default: concentration_fermi_default.npz')

    args = parser.parse_args()
    if not args.output_file.endswith('.npz'):
        parser.error('Output file must end with .npz')

    times_post = np.linspace(args.time_start, args.time_start + args.duration, args.nb_sampled)

    c_art, c_myo, c_inf = get_concentrations(
        times_post=times_post,
        alpha_myo=args.alpha_myo,
        beta_myo=args.beta_myo,
        alpha_inf=args.alpha_inf,
        beta_inf=args.beta_inf,
        lv_myo_shift=args.lv_myo_shift,
        lv_inf_shift=args.lv_inf_shift,
        contrast_dose=args.contrast_dose,
        is_stress=args.is_stress,
        time_step=args.time_step,
        plot=args.plot,
    )
    np.savez_compressed(args.output_file, c_art=c_art, c_myo=c_myo, c_inf=c_inf, times_post=times_post)
    print(f"Concentrations saved to {args.output_file}")


if __name__ == "__main__":
    main()
