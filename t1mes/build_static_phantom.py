"""Build a static digital T1MES phantom and save it as an MRzero-compatible .npz file.

Command line: ``t1mes_build_static --help``.
"""
import argparse

import numpy as np

from t1mes.convert_arguments import add_phantom_arguments, finalize_phantom_arguments, get_tubes_ID_from_B0
from t1mes.phantom_utils import (
    JSON_PARAMETERS,
    build_label_mask,
    build_tissue_masks,
    compute_field_maps,
    fill_static_map,
    load_phantom_parameters,
    plot_phantom,
)


def build_static_phantom(
    FOV=(200, 200, 50),
    matrix=(200, 200, 10),
    field_strength=1.5,
    ncoils=8,
    param_json=JSON_PARAMETERS,
    b0field_kwargs=None,
    b1field_kwargs=None,
    sensmap_kwargs=None,
    phantom_file='T1MES_default.npz',
    plot=True,
    plot_kwargs=None,
    rotate_LR=False,
):
    """Create and save a static T1MES phantom (Captur et al. 2016, doi:10.1186/s12968-016-0280-z).

    Args:
        FOV (sequence of 3 float): field of view [mm].
        matrix (sequence of 3 int): matrix size (number of voxels).
        field_strength (float): 1.5 or 3 [T]; selects the tubes of the 1.5 T or 3 T phantom.
        ncoils (int): number of receive coils of the sensitivity maps.
        param_json (str): JSON file with the tissue parameters. Defaults to the packaged T1MES_values.json.
        b0field_kwargs (dict, optional): extra keyword arguments for ``mrtwin.b0field()``.
        b1field_kwargs (dict, optional): extra keyword arguments for ``mrtwin.b1field()``.
        sensmap_kwargs (dict, optional): extra keyword arguments for ``mrtwin.sensmap()``.
        phantom_file (str): output .npz file.
        plot (bool): plot the phantom with ``MRzeroCore.VoxelGridPhantom.plot``.
        plot_kwargs (dict, optional): extra keyword arguments for the plot method.
        rotate_LR (bool): rotate the phantom by 90 degrees in the xy-plane.

    Returns:
        dict: the arrays written to ``phantom_file`` (times in s, FOV in m).
    """
    phantom_data = load_phantom_parameters(param_json)
    tubes_ID = get_tubes_ID_from_B0(field_strength)
    mask, FOV = build_label_mask(FOV, matrix, rotate_LR)

    # Tissue property maps (times in ms as in the JSON file)
    T1_map, T2_map, T2dash_map, PD_map, chi_map, D_map = (
        fill_static_map(mask, tubes_ID, phantom_data, key)
        for key in ('T1', 'T2', 'T2dash', 'PD', 'chi', 'D')
    )

    B0_map, B1_map, coil_sens = compute_field_maps(
        PD_map, chi_map, field_strength, ncoils, b0field_kwargs, b1field_kwargs, sensmap_kwargs,
    )

    # MRzero expects times in s and the FOV in m
    maps = dict(
        PD_map=PD_map,
        T1_map=T1_map * 1e-3,
        T2_map=T2_map * 1e-3,
        T2dash_map=T2dash_map * 1e-3,
        D_map=D_map,
        B0_map=B0_map,
        B1_map=B1_map,
        FOV=np.array(FOV, dtype=float) * 1e-3,
        coil_sens=coil_sens,
        **build_tissue_masks(mask, tubes_ID, phantom_data),
    )
    np.savez_compressed(phantom_file, **maps)
    print(f"Static phantom saved to {phantom_file}")

    if plot:
        plot_phantom(phantom_file, plot_kwargs, dynamic=False)
    return maps


#: Backward-compatible name of :func:`build_static_phantom`.
get_phantom_values = build_static_phantom


def main():
    parser = argparse.ArgumentParser(
        description="Build a static T1MES phantom (.npz) for the MRzero simulator.")
    add_phantom_arguments(parser, JSON_PARAMETERS)
    args = finalize_phantom_arguments(parser, parser.parse_args())

    build_static_phantom(
        FOV=args.FOV,
        matrix=args.matrix,
        field_strength=args.field_strength,
        ncoils=args.ncoils,
        param_json=args.param_json,
        b0field_kwargs=args.b0_kwargs,
        b1field_kwargs=args.b1_kwargs,
        sensmap_kwargs=args.sensmap_kwargs,
        phantom_file=args.output_file,
        plot=args.plot,
        plot_kwargs=args.plot_kwargs,
        rotate_LR=args.rotate_LR,
    )


if __name__ == "__main__":
    main()
