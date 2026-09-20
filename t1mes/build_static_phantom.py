import argparse
import json
import importlib.resources as pkg_resources

import numpy as np
import MRzeroCore as mr0
from mrtwin import b1field, b0field, sensmap

from t1mes.generate_mask import create_t1mes_phantom_mask
from t1mes.convert_arguments import get_tubes_ID_from_B0, parse_key_value_or_json_file, complete_imaging_args

JSON_PARAMETERS = pkg_resources.files("t1mes").joinpath("T1MES_values.json").as_posix()


def get_phantom_values(
    FOV=(200, 200, 200),
    matrix=(128,128,5),
    field_strength=1.5,
    ncoils=8,
    param_json=JSON_PARAMETERS,
    b0field_kwargs={},
    b1field_kwargs={},
    sensmap_kwargs={}, 
    phantom_file='T1MES_default.npz',
    plot=True, 
    plot_kwargs={},
    rotate_LR=False,
):
    """
    Create and save a T1MES phantom using values from doi:10.1186/s12968-016-0280-z.

    Args:
        FOV (list, tuple): field of view of the phantom [mm].
        matrix (list, tuple): Resolution matrix.
        field_strength (float): Field strength in Tesla. Only 1.5 or 3.0 supported.
        ncoils (int): Number of receiver coils for sensitivity map generation. Defaults to 8.
        param_json (str): Path to JSON file with phantom parameters. Defaults to 't1mes/T1MES_values.json'.
        b0field_kwargs (dict, optional): Keyword arguments passed to `mrtwin.b0field()`.
        b1field_kwargs (dict, optional): Keyword arguments passed to `mrtwin.b1field()`.
        sensmap_kwargs (dict, optional): Keyword arguments passed to `mrtwin.sensmap()`.
        phantom_file (str): Path to save the generated phantom .npz file. Defaults to 'T1MES_default.npz'.
        plot (bool): If True, displays the phantom using MRzeroCore. Defaults to True.
        plot_kwargs (dict, optional): Keyword arguments passed to `mr0.DynamicVoxelPhantom.plot()`.
        rotate_LR (bool): Rotate the phantom by 90 degrees in the LR plane.
    """
    # Load T1MES phantom tissue parameters from JSON
    phantom_data = json.load(open(param_json, 'r'))  # Values from Bruker minispec (1.4 T, 22°C)

    # Generate labeled mask (3D array with labels: 1–9 for tubes, 10 for resin)
    mask = create_t1mes_phantom_mask(FOV, matrix)
    
    # Optional 90° rotation
    if rotate_LR:
        mask = np.rot90(mask, k=1, axes=(0, 1))

    tubes_ID = get_tubes_ID_from_B0(field_strength)

    # Initialize property maps
    T1_map = np.zeros_like(mask, dtype=float)
    T2_map = np.zeros_like(mask, dtype=float)
    PD_map = np.zeros_like(mask, dtype=float)
    chi_map = np.zeros_like(mask, dtype=float)
    D_map = np.zeros_like(mask, dtype=float)
    T2dash_map = np.zeros_like(mask, dtype=float)

    # Fill values for each tube (labels 1 to 9)
    for idx, tube_id in enumerate(tubes_ID):
        label = idx + 1
        T1_map[mask == label] = phantom_data[tube_id]['T1']
        T2_map[mask == label] = phantom_data[tube_id]['T2']
        T2dash_map[mask == label] = phantom_data[tube_id]['T2dash']
        PD_map[mask == label] = phantom_data[tube_id]['PD']
        chi_map[mask == label] = phantom_data[tube_id]['chi']
        D_map[mask == label] = phantom_data[tube_id]['D']

    # Fill values for the outer resin matrix (label 10)
    T1_map[mask == 10] = phantom_data['OuterMatrix']['T1']
    T2_map[mask == 10] = phantom_data['OuterMatrix']['T2']
    T2dash_map[mask == 10] = phantom_data['OuterMatrix']['T2dash']
    PD_map[mask == 10] = phantom_data['OuterMatrix']['PD']
    chi_map[mask == 10] = phantom_data['OuterMatrix']['chi']
    D_map[mask == 10] = phantom_data['OuterMatrix']['D']

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
    sensmap_kwargs['shape'] = (ncoils,) + PD_map.squeeze().shape
    coil_sens = sensmap(**sensmap_kwargs)
    # Ensure coil_sens has 4 dimensions
    if coil_sens.ndim == 3:
        coil_sens = np.expand_dims(coil_sens, axis=-1)
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
        T2dash_map=T2dash_map * 1e-3,
        D_map=D_map,  # Diffusion is set to 0 everywhere
        B0_map=B0_map,
        B1_map=B1_map,
        FOV=np.array(FOV, dtype=float)*1e-3, # in meters
        coil_sens=coil_sens,
        **tubes_mask,
    )

    # Optional plot: display center slice using MRzeroCore
    if plot:
        phantom = mr0.VoxelGridPhantom.load(phantom_file)
        if 'time_unit' not in plot_kwargs:
            plot_kwargs['time_unit'] = 'ms'
        if 'display_units' not in plot_kwargs:
            plot_kwargs['display_units'] = True
        phantom.plot(**plot_kwargs)


def main():
    parser = argparse.ArgumentParser(description="Build Static T1MES Phantom")

    # Imaging Parameters
    imaging_group = parser.add_argument_group("Imaging Parameters")
    imaging_group.add_argument('-FOV', type=float, nargs=3, help='Field of view in mm. Default to (200,200,50)')   
    imaging_group.add_argument('-r', '--resolution', type=float, nargs=3, help='Voxel resolution. Default to 1x1x5 mm')
    imaging_group.add_argument('-m', '--matrix', type=int, nargs=3, help='Matrix size. Default to 200x200x10')
    imaging_group.add_argument('-B0', '--field_strength', type=float, help='Main magnetic field strength in Tesla. Default to 1.5', default=1.5)
    imaging_group.add_argument('--ncoils', type=int, help='Number of receiver coils. Default to 8', default=8)
    imaging_group.add_argument('--LR', action=argparse.BooleanOptionalAction, help="Whether to rotate the phantom by 90 degrees. Default to False", default=False)

    # Phantom Configuration
    phantom_group = parser.add_argument_group("Phantom Configuration")
    phantom_group.add_argument('--param_json', type=str, help='Path to JSON file used to generate phantom', default=JSON_PARAMETERS)
    phantom_group.add_argument('--b0_kwargs', type=parse_key_value_or_json_file, nargs='+', help='JSON string or file path for B0 simulation params', default={})
    phantom_group.add_argument('--b1_kwargs', type=parse_key_value_or_json_file, nargs='+', help='JSON string or file path for B1 simulation params', default={})
    phantom_group.add_argument('--sensmap_kwargs', type=parse_key_value_or_json_file, nargs='+', help='JSON string or file path for sensitivity map params', default={})

    # Output
    output_group = parser.add_argument_group("Output")
    output_group.add_argument('-p', '--phantom_file', type=str, help='Path to output phantom file', default='T1MES_default.npz')

    # Visualization
    vis_group = parser.add_argument_group("Visualization")
    vis_group.add_argument('--plot', action=argparse.BooleanOptionalAction, help="Whether to plot resulting phantom. Default to True", default=True)
    vis_group.add_argument('--plot_kwargs', type=parse_key_value_or_json_file, nargs='+', help='JSON string or file path with plot customization options', default={})

    args = parser.parse_args()
    
    args = complete_imaging_args(args)
    
    print(f"FOV: {args.FOV}")
    print(f"Resolution: {args.resolution}")
    print(f"Matrix: {args.matrix}")

    for attr in ['b0_kwargs', 'b1_kwargs', 'sensmap_kwargs', 'plot_kwargs']:
        val = getattr(args, attr)
        if val:
            setattr(args, attr, {k: v for d in val for k, v in d.items()})

    get_phantom_values(
        FOV=args.FOV,
        matrix=args.matrix,
        field_strength=args.field_strength,
        ncoils=args.ncoils,
        param_json=args.param_json,
        b0field_kwargs=args.b0_kwargs,
        b1field_kwargs=args.b1_kwargs,
        sensmap_kwargs=args.sensmap_kwargs,
        phantom_file=args.phantom_file,
        plot=args.plot,
        plot_kwargs=args.plot_kwargs,
        rotate_LR=args.LR,
    )


if __name__ == "__main__":
    main()