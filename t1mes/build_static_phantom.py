import argparse
import os
import json

import numpy as np
import MRzeroCore as mr0
from mrtwin import b1field, b0field, sensmap

from t1mes.generate_mask import create_t1mes_phantom_mask
from t1mes.convert_arguments import load_T2dash_map, parse_key_value, get_tubes_ID_from_B0

JSON_PARAMETERS = 't1mes/T1MES_values.json'


def get_phantom_values(
    FOV=(200, 200, 200),
    resolution=(128,128,5),
    field_strength=1.5, 
    PD=None,
    chi=None, 
    D=None, 
    phantom_file=None, 
    plot=True, 
    T2dash_map=3e-2,
    ncoils=8,
    b0field_kwargs={},
    b1field_kwargs={},
    sensmap_kwargs={},
):
    """
    Create and save a T1MES phantom using values from doi:10.1186/s12968-016-0280-z.

    Args:
        FOV (list, tuple): field of view of the phantom [mm].
        resolution (list, tuple): Resolution matrix.
        field_strength (float): Field strength in Tesla. Only 1.5 or 3.0 supported.
        PD (dict): Proton density for 'tubes' and 'OuterMatrix'. Defaults to {1.0, 0.5}.
        D (dict): Diffusion coefficient for 'tubes' and 'OuterMatrix'. Defaults to {.0, 0.}.
        chi (dict): magnetic susceptibility for 'tubes' and 'OuterMatrix'. Defaults to {-9., -8.}.
        phantom_file (str): Output .npz file to save phantom data. Defaults to 'T1MES_default.npz'.
        plot (bool): If True, displays a center slice of the phantom using MRzeroCore.
        T2dash_map (float or array): T2' map. Defaults to 30 ms.
        ncoils (int): Number of receiver coils for sensitivity map generation. Defaults to 8.
        b0field_kwargs (dict, optional): Keyword arguments passed to `mrtwin.b0field()`.
        b1field_kwargs (dict, optional): Keyword arguments passed to `mrtwin.b1field()`.
        sensmap_kwargs (dict, optional): Keyword arguments passed to `mrtwin.sensmap()`.
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
    mask = create_t1mes_phantom_mask(FOV, resolution)

    tubes_ID = get_tubes_ID_from_B0(field_strength)

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
    D_map[mask == 10] = D['OuterMatrix']

    # compute B0
    b0field_kwargs['chi'] = chi_map
    b0field_kwargs['B0'] = field_strength
    if 'b0range' not in b0field_kwargs:
        b0field_kwargs['b0range'] = (-200,200) # Default range of -200/+200Hz
    B0_map = b0field(**b0field_kwargs)
    
    # Compute B1
    b1field_kwargs['shape'] = B0_map.squeeze().shape
    b1field_kwargs['mask'] = PD_map.squeeze()>0
    B1_map = b1field(**b1field_kwargs)
    if B1_map.ndim==2:
        B1_map = np.expand_dims(B1_map, axis=2)
    
    # Compute sensitivity map
    sensmap_kwargs['shape'] = (ncoils,) + PD_map.shape
    coil_sens = sensmap(**sensmap_kwargs)
    
    # Add tissue masks
    tubes_mask = {phantom_data[tube_id]['description']: mask==i+1  for i, tube_id in enumerate(tubes_ID)}
    tubes_mask['OuterMatrix'] = mask==10
    tubes_mask = {'tissue_'+k: v for k,v in tubes_mask.items()}
    
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
        FOV=(0.2, 0.2, 0.2),
        coil_sens=coil_sens,
        **tubes_mask,
    )

    # Optional plot: display center slice using MRzeroCore
    if plot:
        mr0.DynamicVoxelPhantom.load(phantom_file).plot(time_unit='ms', display_units=True)


def main():
    parser = argparse.ArgumentParser('Build Static T1MES Phantom')
    parser.add_argument('-FOV', help='Field of View of the phantom in (mm). Default to (200,200,200)', type=float, default=(200,200,200), nargs='+')
    parser.add_argument('-r', '--resolution', help='Resolution of the phantom. Default to (128,128,5)', type=float, default=(128,128,5), nargs='+')
    parser.add_argument('-B0', '--field_strength', help='Static Field strength in Tesla. Default to 1.5', type=float, default=1.5)
    parser.add_argument('-PD', help="Proton density for 'tubes' and 'OuterMatrix' as comma-separated values (e.g., '1.0,0.5')")
    parser.add_argument('-chi', help="Susceptibility for 'tubes' and 'OuterMatrix' as comma-separated values (e.g., -chi='-9.,-8.')"+\
        "Use = followed by quotes to prevent negative values from being interpreted as flags.")
    parser.add_argument('-D', help="Diffusion coefficient for 'tubes' and 'OuterMatrix' as comma-separated values (e.g., '0.,0.')")
    parser.add_argument('-p', '--phantom_file', help='Output MRzero phantom file (.npz)')
    parser.add_argument('--plot', help="Whether to plot resulting phantom", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('-T2dash', help="T2' value. Currently supports uniform T2′ only.", default=30e-3)
    parser.add_argument('--ncoils', help="Number of receiver coils for sensitivity map generation. Default to 8", type=int, default=8)
    parser.add_argument('--b0_kwargs', help="Keyword arguments passed to `mrtwin.b0field()`.", type=parse_key_value, default={}, nargs='+')
    parser.add_argument('--b1_kwargs', help="Keyword arguments passed to `mrtwin.b1field()`.", type=parse_key_value, default={}, nargs='+')
    parser.add_argument('--sensmap_kwargs', help="Keyword arguments passed to `mrtwin.sensmap()`.", type=parse_key_value, default={}, nargs='+')

    args = parser.parse_args()
    
    assert len(args.FOV)==3, 'You must provide the FOV for the 3 dimensions'
    assert len(args.resolution)==3, 'You must provide the resolution for the 3 dimensions'
    locations = ['tubes', 'OuterMatrix']
    if args.PD:
        assert len(args.PD.split(','))==2, \
            "You should provide two comma-separated values for 'tubes' and 'OuterMatrix' (e.g., 1.0,0.5)"
        args.PD = dict(zip(locations, args.PD.split(',')))
    if args.chi:
        assert len(args.chi.split(','))==2, \
            "You should provide two comma-separated values for 'tubes' and 'OuterMatrix' (e.g., ' -9.,-8.')"
        args.chi = dict(zip(locations, args.chi.split(',')))
    if args.D:
        assert len(args.D.split(','))==2, \
            "You should provide two comma-separated values for 'tubes' and 'OuterMatrix' (e.g., 0.,0.)"
        args.D = dict(zip(locations, args.D.split(',')))
    T2dash_map = load_T2dash_map(args.T2dash)
    
    get_phantom_values(
        args.FOV,
        args.resolution,
        args.field_strength, 
        args.PD,
        args.chi,
        args.D,
        args.phantom_file,
        args.plot,
        T2dash_map,
        args.ncoils,
        dict(args.b0_kwargs),
        dict(args.b1_kwargs),
        dict(args.sensmap_kwargs),
    )


if __name__ == "__main__":
    main()