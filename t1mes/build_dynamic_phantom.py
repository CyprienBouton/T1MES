import json
import argparse
import numpy as np
import MRzeroCore as mr0
from mrtwin import b1field, b0field

from t1mes.generate_mask import create_t1mes_phantom_mask
from t1mes.convert_arguments import load_T2dash_map, parse_key_value, get_tubes_ID_from_B0
from t1mes.compute_dynamic_parameters import get_T1_T2_over_time


JSON_PARAMETERS = 't1mes/T1MES_values.json'


def str_to_seconds(time_str: str) -> float:
    """Convert a time string (M or M:S) to total seconds."""
    if ':' in time_str:
        try:
            minutes, seconds = map(float, time_str.split(':'))
            return minutes * 60 + seconds
        except ValueError:
            raise ValueError(f"Invalid time format: '{time_str}'. Use M or M:S")
    else:
        try:
            return float(time_str) * 60
        except ValueError:
            raise ValueError(f"Invalid time format: '{time_str}'. Use M or M:S")


def get_phantom_values(
    concentrations_file=None,
    FOV=(200, 200, 200),
    resolution=(128,128,5),
    field_strength=1.5, 
    PD=None,
    chi=None, 
    D=None, 
    phantom_file=None, 
    plot=True, 
    T2dash_map=3e-2,
    times_post=None,
    b0field_kwargs={},
    b1field_kwargs={},  
    plot_kwargs={},
):
    """
    Create and save a T1MES phantom using values from doi:10.1186/s12968-016-0280-z.

    Args:
        concentrations_file (str): file containing contrast agent concentrations for injected tube over time [mM].
        FOV (list, tuple): field of view of the phantom [mm].
        resolution (list, tuple): Resolution matrix.
        field_strength (float): Field strength in Tesla. Only 1.5 or 3.0 supported.
        PD (dict): Proton density for 'tubes' and 'OuterMatrix'. Defaults to {1.0, 0.5}.
        D (dict): Diffusion coefficient for 'tubes' and 'OuterMatrix'. Defaults to {.0, 0.}.
        chi (dict): magnetic susceptibility for 'tubes' and 'OuterMatrix'. Defaults to {-9., -8.}.
        phantom_file (str): Output .npz file to save phantom data. Defaults to 'T1MES_default.npz'.
        plot (bool): If True, displays a center slice of the phantom using MRzeroCore.
        T2dash_map (float or array): T2' map. Defaults to 30 ms.
        times_post (np.ndarray): Timing of each concentration after the injection [s]. Default to None. Use if concentration file is None.
        b0field_kwargs (dict, optional): Keyword arguments passed to `mrtwin.b0field()`.
        b1field_kwargs (dict, optional): Keyword arguments passed to `mrtwin.b1field()`.
        plot_kwargs (dict): Keyword arguments for the plot_dynamic.
    """
    # Load T1MES phantom tissue parameters from JSON
    phantom_data = json.load(open(JSON_PARAMETERS, 'r'))  # Values from Bruker minispec (1.4 T, 22°C)
    
    # Set default phantom file name
    if phantom_file is None:
        phantom_file = 'T1MES_dynamic_default.npz'

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

    # Generate labeled mask (3D array with labels: 1–9 for tubes, 10 for resin)
    mask = create_t1mes_phantom_mask(FOV, resolution)
    T1, T2, times_post = get_T1_T2_over_time(concentrations_file, field_strength, times_post)
    time_frames = len(times_post)
    
    # Initialize property maps
    T1_map = np.zeros((time_frames, *mask.shape))
    T2_map = np.zeros((time_frames, *mask.shape))
    PD_map = np.zeros_like(mask, dtype=float)
    chi_map = np.zeros_like(mask, dtype=float)
    D_map = np.zeros_like(mask, dtype=float)

    # If T2dash_map is a float, fill a full array with that value
    if isinstance(T2dash_map, float):
        T2dash_map = np.full_like(T2_map, T2dash_map)
    tubes_ID = get_tubes_ID_from_B0(field_strength)

    # Fill values for each tube (labels 1 to 9)
    for idx, tube_id in enumerate(tubes_ID):
        label = idx + 1
        PD_map[mask == label] = PD['tubes']
        chi_map[mask == label] = chi['tubes']
        D_map[mask == label] = D['tubes']
        for t_idx in range(time_frames):
            T1_map[t_idx][mask == label] = T1[tube_id][t_idx]
            T2_map[t_idx][mask == label] = T2[tube_id][t_idx]

    # Fill values for the outer resin matrix (label 10)
    PD_map[mask == 10] = PD['OuterMatrix']
    chi_map[mask == 10] = chi['OuterMatrix']
    D_map[mask == 10] = D['OuterMatrix']
    T1_map[:, mask == 10] = T1['OuterMatrix'][:, None]
    T2_map[:, mask == 10] = T2['OuterMatrix'][:, None]


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
        B1_map = np.expand_dims(B1_map, -1)
    
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
        time_points=times_post,
        **tubes_mask,
    )

    # Optional plot: display center slice using MRzeroCore
    if plot:
        phantom = mr0.DynamicVoxelPhantom.load(phantom_file)
        if 'time_unit' not in plot_kwargs:
            plot_kwargs['time_unit'] = 'ms'
        if 'display_units' not in plot_kwargs:
            plot_kwargs['display_units'] = True
        phantom.plot_dynamic(**plot_kwargs)


def main():
    parser = argparse.ArgumentParser('Build Static T1MES Phantom')
    parser.add_argument('-c', '--concentrations_file', help='File containing contrast agent concentrations over time [mM] for injected tubes', default=None)
    parser.add_argument('-FOV', help='Field of View of the phantom in (mm). Default to (200,200,200)', type=float, default=(200,200,200), nargs='*')
    parser.add_argument('-r', '--resolution', help='Resolution of the phantom. Default to (128,128,5)', type=float, default=(128,128,5), nargs='*')
    parser.add_argument('-B0', '--field_strength', help='Static Field strength in Tesla', type=float, default=1.5)
    parser.add_argument('-PD', help="Proton density for 'tubes' and 'OuterMatrix' as comma-separated values (e.g., '1.0,0.5')")
    parser.add_argument('-chi', help="Susceptibility for 'tubes' and 'OuterMatrix' as comma-separated values (e.g., -chi='-9.,-8.')"+\
        "Use = followed by quotes to prevent negative values from being interpreted as flags.")
    parser.add_argument('-D', help="Diffusion coefficient for 'tubes' and 'OuterMatrix' as comma-separated values (e.g., '0.,0.')")
    parser.add_argument('-p', '--phantom_file', help="Output MRzero phantom file (.npz). Default to 'T1MES_dynamic_default.npz'")
    parser.add_argument('--plot', help="Whether to plot resulting phantom", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('-T2dash', help="T2' value. Currently supports uniform T2′ only.", default=30e-3)
    parser.add_argument('--times_post', help='Timing post injection of the T1 and T2 relaxation times [s].', nargs='*', type=str_to_seconds,)
    
    parser.add_argument('--b0_kwargs', help="Keyword arguments passed to `mrtwin.b0field()`.", type=parse_key_value, default={}, nargs='+')
    parser.add_argument('--b1_kwargs', help="Keyword arguments passed to `mrtwin.b1field()`.", type=parse_key_value, default={}, nargs='+')
    parser.add_argument('--plot_kwargs', help="Keyword arguments passed to `DynamicVoxelPhantom.plot_dynamic()`.", type=parse_key_value, default={}, nargs='+')
    
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
        concentrations_file=args.concentrations_file,
        FOV=args.FOV,
        resolution=args.resolution,
        field_strength=args.field_strength,
        PD=args.PD,
        chi=args.chi,
        D=args.D,
        phantom_file=args.phantom_file,
        plot=args.plot,
        T2dash_map=T2dash_map,
        times_post=args.times_post,
        b0field_kwargs=args.b0_kwargs,
        b1field_kwargs=args.b1_kwargs,
        plot_kwargs=args.plot_kwargs,
    )


if __name__ == "__main__":
    main()