import json
import numpy as np
import MRzeroCore as mr0
from mrtwin import b1field, b0field

from generate_mask import create_t1mes_phantom_mask

JSON_PARAMETERS = 'T1MES_values.json'

def get_phantom_values(
    field_strength=1.5, 
    PD=None,
    chi=None, 
    D=None, 
    phantom_file=None, 
    plot=True, 
    T2dash_map=3e-2,
):
    """
    Create and save a T1MES phantom using values from doi:10.1186/s12968-016-0280-z.

    Args:
        field_strength (float): Field strength in Tesla. Only 1.5 or 3.0 supported.
        PD (dict): Proton density values for 'tubes' and 'OuterMatrix'. Defaults to {1.0, 0.5}.
        D (dict): Diffusion coefficient values for 'tubes' and 'OuterMatrix'. Defaults to {.0, 0.}.
        chi (dict): magnetic susceptibility for 'tubes' and 'OuterMatrix'. Defaults to {-9., -8.}.
        phantom_file (str): Output .npz file to save phantom data. Defaults to 'T1MES_default.npz'.
        plot (bool): If True, displays a center slice of the phantom using MRzeroCore.
        T2dash_map (float or array): Value or array for T2* map (used by MRzero). Defaults to 30 ms.
    """
    
    # Set default phantom file name
    if phantom_file is None:
        phantom_file = 'T1MES_default.npz'

    # Set default PD if not provided
    if PD is None:
        PD = {
            'tubes': 1.0,
            'OuterMatrix': 0.5,
        }
        
    # Set default chi_map if not provided
    if chi is None:
        chi = {
            'tubes': -9.,
            'OuterMatrix': -8.,
        }
    if D is None:
        D = {
            'tubes': 0.,
            'OuterMatrix': 0.,
        }
    # Load T1MES phantom tissue parameters from JSON
    phantom_data = json.load(open(JSON_PARAMETERS, 'r'))  # Values from Bruker minispec (1.4 T, 22°C)

    # Generate labeled mask (3D array with labels: 1–9 for tubes, 10 for resin)
    mask = create_t1mes_phantom_mask()

    # Define tube labels based on field strength
    if field_strength == 1.5:
        tubes_ID = [
            'I', 'F', 'C', 
            'H', 'E', 'B',
            'G', 'D', 'A',
        ]
    elif field_strength == 3.0:
        tubes_ID = [
            'I', 'M', 'C', 
            'H', 'L', 'J',
            'G', 'K', 'A',
        ]
    else:
        raise ValueError('field strength must be 1.5 or 3')

    # Initialize property maps
    T1_map = np.zeros_like(mask, dtype=float)
    T2_map = np.zeros_like(mask, dtype=float)
    PD_map = np.zeros_like(mask, dtype=float)
    chi_map = np.zeros_like(mask, dtype=float)
    D_map = np.zeros_like(mask, dtype=float)

    # If T2dash_map is a float, fill a full array with that value
    if isinstance(T2dash_map, float):
        T2dash_map = np.full_like(T2_map, T2dash_map)

    # Fill values for each tube (labels 1 to 9)
    for idx, tube_id in enumerate(tubes_ID):
        label = idx + 1
        T1_map[mask == label] = phantom_data[tube_id]['T1']
        T2_map[mask == label] = phantom_data[tube_id]['T2']
        PD_map[mask == label] = PD['tubes']
        chi_map[mask == label] = chi['tubes']
        D_map[mask == label] = D['tubes']

    # Fill values for the outer resin matrix (label 10)
    T1_map[mask == 10] = phantom_data['OuterMatrix']['T1']
    T2_map[mask == 10] = phantom_data['OuterMatrix']['T2']
    PD_map[mask == 10] = PD['OuterMatrix']
    chi_map[mask == 10] = chi['OuterMatrix']
    D_map[mask == label] = D['OuterMatrix']

    # compute B0 and B1
    B0_map = b0field(chi_map, (-200,200), field_strength)
    # B1_map = b1field(B0_map.shape)
    B1_map = b1field(B0_map.shape, b1range=(0.8, 1.2), nrings=1)
    
    # Save maps to compressed .npz file (convert T1, T2 to seconds)
    np.savez_compressed(
        phantom_file,
        PD_map=PD_map,
        T1_map=T1_map * 1e-3,
        T2_map=T2_map * 1e-3,
        T2dash_map=T2dash_map,
        D_map=D_map,  # Diffusion is set to 0 everywhere
        B0_map=B0_map,
        B1_map=B1_map,
    )

    # Optional plot: display center slice using MRzeroCore
    if plot:
        mr0.VoxelGridPhantom.load(phantom_file).plot(plot_slice='center')
    
if __name__ == "__main__":
    get_phantom_values(3)
    