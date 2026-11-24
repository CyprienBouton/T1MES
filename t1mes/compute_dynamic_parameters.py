"""Time-varying T1/T2 of the T1MES tubes for the dynamic (contrast-injection) phantom.

Two ways of describing the dynamics are supported:

1. **Concentration curves** (``concentrations_file``, e.g. produced by
   ``t1mes_simulate_fermi``): the five post-GBCA tubes are given the native
   T1/T2 of their matching native tube, shortened by the contrast agent with
   the fast-exchange relaxivity model

       1/T_i(t) = 1/T_i,native + r_i * C(t),   i = 1, 2

   where ``r_i`` is the relaxivity [s^-1 mM^-1] and ``C(t)`` the gadolinium
   concentration [mM]. Blood tubes follow the arterial curve, the
   "Long post-GBCA myocardium" tube follows the infarct curve and the other
   post-GBCA myocardium tubes follow the myocardium curve.

2. **Explicit time series in the JSON file** (``times_post``): any ``T1`` or
   ``T2`` field may be a list with one value per entry of ``times_post``;
   scalar fields are kept constant.
"""
from numbers import Real

import numpy as np

from t1mes.convert_arguments import get_tubes_ID_from_B0
from t1mes.phantom_utils import JSON_PARAMETERS, load_phantom_parameters

# Default relaxivities, kept at the historical values used to generate all existing
# dynamic phantoms (same defaults as pymrzeroxcat).
# KNOWN ISSUE (to be decided): times are converted to seconds, so these values act
# as s^-1 mM^-1 and are ~1000x smaller than published Gd-DTPA relaxivities
# (r1 ~ 4.1, r2 ~ 4.6 s^-1 mM^-1 in plasma at 1.5 T, Rohrer et al., Invest Radiol 2005).
# With the published values the Fermi concentration curves give unrealistic myocardial
# T1 (~50 ms): the concentration model must be calibrated before changing them.
# Override from the CLI with --r1 / --r2.
DEFAULT_R1 = 4.1e-3
DEFAULT_R2 = 4.6e-3


def compute_current_relaxation(T_native, relaxivity, concentration, time_unit='ms'):
    """Relaxation time in presence of a contrast agent (fast-exchange model)::

        1/T_post = 1/T_native + relaxivity * concentration

    Args:
        T_native (float or np.ndarray): native relaxation time, in ``time_unit``.
        relaxivity (float): relaxivity of the contrast agent [s^-1 mM^-1].
        concentration (float or np.ndarray): contrast agent concentration [mM].
        time_unit (str): unit of ``T_native`` and of the result, 'ms' or 's'.

    Returns:
        float or np.ndarray: relaxation time after injection, in ``time_unit``.
    """
    assert time_unit in ['ms', 's'], "time_unit should be either 'ms' or 's'"
    time_factor = 1e-3 if time_unit == 'ms' else 1
    T_native_s = np.asarray(T_native, dtype=float) * time_factor
    T_post_s = 1 / (1 / T_native_s + relaxivity * np.asarray(concentration, dtype=float))
    return T_post_s / time_factor


def get_matching_native_tube(tube_description, phantom_data, field_strength=1.5):
    """Return the ID of the native tube matching a post-GBCA tube.

    Post-GBCA blood tubes map to "Normal native blood <B0>"; post-GBCA
    myocardium tubes map to the native tube with the same short/medium/long
    qualifier at the given field strength.
    """
    field_str = ' 1.5 T' if field_strength == 1.5 else ' 3 T'
    if 'blood' in tube_description:
        native_tube_str = 'Normal native blood' + field_str
    else:
        native_tube_str = tube_description.replace('post-GBCA', 'native') + field_str
    native_tube_ID = next(
        (k for k, v in phantom_data.items() if v.get('description') == native_tube_str), None)
    if native_tube_ID is None:
        raise KeyError(f"No native tube with description '{native_tube_str}' in the JSON file")
    return native_tube_ID


def get_concentration_dict(concentrations_file, tubes_ID_injected, phantom_data):
    """Assign a concentration curve to each injected tube.

    Args:
        concentrations_file (str): .npz file with ``c_art``, ``c_myo``, ``c_inf`` [mM]
            and ``times_post`` [s] (output of ``t1mes_simulate_fermi``).
        tubes_ID_injected (list): IDs of the post-GBCA tubes.
        phantom_data (dict): tissue parameters loaded from the JSON file.

    Returns:
        tuple: ({tube_ID: concentration over time [mM]}, times_post [s]).
    """
    with np.load(concentrations_file) as data:
        c_art, c_myo, c_inf = data['c_art'], data['c_myo'], data['c_inf']
        times_post = data['times_post']
    concentration_tubes = {}
    for tube_ID in tubes_ID_injected:
        description = phantom_data[tube_ID]['description'].lower()
        if 'blood' in description:
            concentration_tubes[tube_ID] = c_art
        elif 'long' in description:
            concentration_tubes[tube_ID] = c_inf
        else:
            concentration_tubes[tube_ID] = c_myo
    return concentration_tubes, times_post


def get_T1_T2_over_time(
    concentrations_file=None,
    field_strength=1.5,
    times_post=None,
    r1=DEFAULT_R1,
    r2=DEFAULT_R2,
    param_json=JSON_PARAMETERS,
):
    """T1 and T2 of every tube (and of the outer matrix) over time.

    Args:
        concentrations_file (str, optional): .npz file with contrast agent
            concentrations over time (see :func:`get_concentration_dict`).
        field_strength (float): 1.5 or 3 [T].
        times_post (np.ndarray, optional): time after injection of each frame [s].
            Used only if ``concentrations_file`` is None; the T1/T2 time series
            are then read from ``param_json``.
        r1 (float or dict): longitudinal relaxivity [s^-1 mM^-1], or ``{tube_ID: r1}``.
        r2 (float or dict): transverse relaxivity [s^-1 mM^-1], or ``{tube_ID: r2}``.
        param_json (str): JSON file with the tissue parameters.

    Returns:
        dict: ``{tube_ID: T1 over time [ms]}`` (including 'OuterMatrix').
        dict: ``{tube_ID: T2 over time [ms]}`` (including 'OuterMatrix').
        np.ndarray: time after injection of each frame [s].
    """
    phantom_data = load_phantom_parameters(param_json)
    tubes_ID = get_tubes_ID_from_B0(field_strength)

    if concentrations_file is None:
        if times_post is None:
            raise ValueError('Either concentrations_file or times_post must be provided')
        return read_T1_T2_over_time(tubes_ID, times_post, param_json=param_json)

    tubes_ID_injected = [ID for ID in tubes_ID if 'post-GBCA' in phantom_data[ID]['description']]
    concentrations, times_post = get_concentration_dict(concentrations_file, tubes_ID_injected, phantom_data)
    time_frames = len(times_post)

    if isinstance(r1, Real):
        r1 = {ID: r1 for ID in tubes_ID}
    if isinstance(r2, Real):
        r2 = {ID: r2 for ID in tubes_ID}

    T1_over_time, T2_over_time = {}, {}
    for tube_ID in tubes_ID + ['OuterMatrix']:
        if tube_ID in tubes_ID_injected:
            # Post-GBCA tube: native values of the matching tube shortened by the contrast agent
            native = phantom_data[get_matching_native_tube(
                phantom_data[tube_ID]['description'], phantom_data, field_strength)]
            c = concentrations[tube_ID]
            T1_over_time[tube_ID] = compute_current_relaxation(native['T1'], r1[tube_ID], c)
            T2_over_time[tube_ID] = compute_current_relaxation(native['T2'], r2[tube_ID], c)
        else:
            # Native tubes and outer matrix: constant
            T1_over_time[tube_ID] = np.full(time_frames, phantom_data[tube_ID]['T1'], dtype=float)
            T2_over_time[tube_ID] = np.full(time_frames, phantom_data[tube_ID]['T2'], dtype=float)

    return T1_over_time, T2_over_time, times_post


def read_T1_T2_over_time(tubes_ID, times_post, param_json=JSON_PARAMETERS):
    """Read T1/T2 time series from a JSON file.

    Scalar ``T1``/``T2`` fields are repeated over time; list fields must have
    one value per entry of ``times_post``.

    Args:
        tubes_ID (list): IDs of the tubes in the phantom.
        times_post (np.ndarray): time after injection of each frame [s].
        param_json (str): JSON file with the tissue parameters.

    Returns:
        dict: ``{tube_ID: T1 over time [ms]}`` (including 'OuterMatrix').
        dict: ``{tube_ID: T2 over time [ms]}`` (including 'OuterMatrix').
        np.ndarray: ``times_post`` [s].
    """
    phantom_data = load_phantom_parameters(param_json)
    times_post = np.asarray(times_post, dtype=float)
    time_frames = len(times_post)

    def as_time_series(tube_ID, key):
        value = phantom_data[tube_ID][key]
        if np.isscalar(value):
            return np.full(time_frames, value, dtype=float)
        if len(value) != time_frames:
            raise ValueError(
                f"'{key}' of tube '{tube_ID}' has {len(value)} values but {time_frames} "
                f"time points were given with times_post")
        return np.asarray(value, dtype=float)

    T1_over_time = {ID: as_time_series(ID, 'T1') for ID in tubes_ID + ['OuterMatrix']}
    T2_over_time = {ID: as_time_series(ID, 'T2') for ID in tubes_ID + ['OuterMatrix']}
    return T1_over_time, T2_over_time, times_post


#: Backward-compatible name of :func:`read_T1_T2_over_time`.
read_dynamic_concentration = read_T1_T2_over_time
