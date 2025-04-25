import json
import numpy as np

from t1mes.concentration_models.fermi_concentration import get_concentrations
from t1mes.convert_arguments import get_tubes_ID_from_B0

JSON_PARAMETERS = 't1mes/T1MES_values.json'

def compute_current_relaxation(T_native, relaxivity, concentration, time_unit='ms'):
    """Compute the current relaxation from a native with the following formula:
        1/T_post = 1/T_native + relaxivity * concentration

    Args:
        T_native (float): native relaxation time.
        relaxivity (float): relaxivity of the contrast agent [s^-1.mM^-1].
        concentration (float):  contrast agent concentration [mM].
        relaxivity (float): relaxivity of the contrast agent [s^-1.mM^-1].

    Returns:
        float: relaxation time after injection.
    """
    assert time_unit in ['ms', 's'], "time_unit should be either 'ms' or 's'"
    time_factor = 1e-3 if time_unit=='ms' else 1
    T_native *= time_factor
    T_post = 1 / (1 / T_native + relaxivity * concentration)
    T_post /= time_factor
    return T_post


def get_matching_native_tube(tube_description, phantom_data, field_strength=1.5):
    field_str = ' 1.5 T' if field_strength==1.5 else ' 3 T'
    if 'blood' in tube_description:
        native_tube_str = 'Normal native blood' + field_str 
    else:
        native_tube_str = tube_description.replace('post-GBCA', 'native') + field_str
    native_tube_ID = next((k for k, v in phantom_data.items() if v.get('description')==native_tube_str), None)
    return native_tube_ID


def get_T1_T2_over_time( 
    concentrations_file,
    field_strength=1.5,
    r1=4.1e-3,
    r2=4.6e-3,
    ):
    """Get T1 and T2 relaxation times by label over time.

    Args:
        concentrations_file (str): file containing contrast agent concentrations for injected tube over time [mM].
        field_strength (float, optional): static field strength. Defaults to 1.5 [T].
        r1 (float, dict, optional): T1 relaxation time. Defaults to 4.1e-3 [s^-1.mM^-1].
        r2 (float, dict, optional): T2 relaxation time. Defaults to 4.6e-3 [s^-1.mM^-1].

    Raises:
        ValueError: if field strength is different than 1.5 or 3

    Returns:
        np.ndarray: shape (time_frames, nb_label) T1 relaxation times for each label over time [ms].
        np.ndarray: shape (time_frames, nb_label) T2 relaxation times for each label over time [ms].
    """
    # Load T1MES phantom tissue parameters from JSON
    phantom_data = json.load(open(JSON_PARAMETERS, 'r'))  # Values from Bruker minispec (1.4 T, 22°C)
            
    # Define tube labels based on field strength
    tubes_ID = get_tubes_ID_from_B0(field_strength)
    tubes_ID_injected = [tube_ID for tube_ID in tubes_ID if 'post-GBCA' in phantom_data[tube_ID]['description']]
    
    # Initialize relaxation times 
    concentrations = get_concentration_dict(concentrations_file, tubes_ID_injected, phantom_data)
    time_frames = len(concentrations['G'])
    nb_labels = 10
    T1_over_time = np.zeros((time_frames, nb_labels))
    T2_over_time = np.zeros((time_frames, nb_labels))

    
    # Convert Inputs
    if isinstance(r1, float):
        r1 = {tube_ID:r1 for tube_ID in tubes_ID}
    if isinstance(r2, float):
        r2 = {tube_ID:r2 for tube_ID in tubes_ID}    
    
    for idx, tube_ID in enumerate(tubes_ID): 
        # tubes with contrast agent
        if tube_ID in tubes_ID_injected:
            native_tube_ID = get_matching_native_tube(phantom_data[tube_ID]['description'], phantom_data)
            for t_idx, c in enumerate(concentrations[tube_ID]):
                T1_over_time[t_idx, idx] = compute_current_relaxation(phantom_data[native_tube_ID]['T1'], r1[tube_ID], c)
                T2_over_time[t_idx, idx] = compute_current_relaxation(phantom_data[native_tube_ID]['T2'], r2[tube_ID], c)
        # native tubes
        else:
            T1_over_time[:, idx] = np.repeat( np.expand_dims(phantom_data[tube_ID]['T1'], axis=0), time_frames, axis=0)
            T2_over_time[:, idx] = np.repeat( np.expand_dims(phantom_data[tube_ID]['T2'], axis=0), time_frames, axis=0)
    
    # Outer Matrix
    T1_over_time[:, -1] = np.repeat( np.expand_dims(phantom_data['OuterMatrix']['T1'], axis=0), time_frames, axis=0)
    T2_over_time[:, -1] = np.repeat( np.expand_dims(phantom_data['OuterMatrix']['T2'], axis=0), time_frames, axis=0)
    
    return T1_over_time, T2_over_time            


def get_concentration_dict(concentrations_file, tubes_ID_injected, phantom_data):
    concentration_tubes = {}
    c_art, c_myo, c_inf = np.load(concentrations_file).values()
    for tube_ID in tubes_ID_injected:
        if 'blood' in phantom_data[tube_ID]['description'].lower():
            concentration_tubes[tube_ID] = c_art
        elif 'long' in phantom_data[tube_ID]['description'].lower():
            concentration_tubes[tube_ID] = c_inf
        else:
            concentration_tubes[tube_ID] = c_myo
    return concentration_tubes