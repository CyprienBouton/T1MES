"""Build a dynamic digital T1MES phantom (time-varying T1/T2 after contrast injection).

The output .npz has the same content as the static phantom, except that
``T1_map`` and ``T2_map`` have shape ``(n_frames, x, y, z)`` and an extra
``time_points`` array [s] gives the time after injection of each frame.
Loading, plotting and simulating it requires ``MRzeroCore.DynamicVoxelPhantom``,
which is only provided by the MRzeroCore fork (see README).

Command line: ``t1mes_build_dynamic --help``.
"""
import argparse

import numpy as np

from t1mes.compute_dynamic_parameters import DEFAULT_R1, DEFAULT_R2, get_T1_T2_over_time
from t1mes.convert_arguments import (
    add_phantom_arguments,
    finalize_phantom_arguments,
    get_tubes_ID_from_B0,
    str_to_seconds,
)
from t1mes.generate_mask import OUTER_MATRIX_LABEL
from t1mes.phantom_utils import (
    JSON_PARAMETERS,
    build_label_mask,
    build_tissue_masks,
    compute_field_maps,
    fill_static_map,
    load_phantom_parameters,
    plot_phantom,
    require_dynamic_voxel_phantom,
)


def build_dynamic_phantom(
    concentrations_file=None,
    times_post=None,
    FOV=(200, 200, 50),
    matrix=(200, 200, 10),
    field_strength=1.5,
    ncoils=8,
    param_json=JSON_PARAMETERS,
    b0field_kwargs=None,
    b1field_kwargs=None,
    sensmap_kwargs=None,
    phantom_file='T1MES_dynamic_default.npz',
    plot=True,
    plot_kwargs=None,
    rotate_LR=False,
    r1=DEFAULT_R1,
    r2=DEFAULT_R2,
):
    """Create and save a dynamic T1MES phantom.

    The dynamics are given either by ``concentrations_file`` (gadolinium
    concentration curves, e.g. from ``t1mes_simulate_fermi``) or by T1/T2 time
    series in ``param_json`` sampled at ``times_post``
    (see :mod:`t1mes.compute_dynamic_parameters`).

    Args:
        concentrations_file (str, optional): .npz file with ``c_art``, ``c_myo``,
            ``c_inf`` [mM] and ``times_post`` [s].
        times_post (sequence of float, optional): time after injection of each
            frame [s]. Used only if ``concentrations_file`` is None.
        FOV (sequence of 3 float): field of view [mm].
        matrix (sequence of 3 int): matrix size (number of voxels).
        field_strength (float): 1.5 or 3 [T].
        ncoils (int): number of receive coils of the sensitivity maps.
        param_json (str): JSON file with the tissue parameters.
        b0field_kwargs (dict, optional): extra keyword arguments for ``mrtwin.b0field()``.
        b1field_kwargs (dict, optional): extra keyword arguments for ``mrtwin.b1field()``.
        sensmap_kwargs (dict, optional): extra keyword arguments for ``mrtwin.sensmap()``.
        phantom_file (str): output .npz file.
        plot (bool): animate the phantom with ``MRzeroCore.DynamicVoxelPhantom.plot_dynamic``.
        plot_kwargs (dict, optional): extra keyword arguments for the plot method.
        rotate_LR (bool): rotate the phantom by 90 degrees in the xy-plane.
        r1 (float or dict): longitudinal relaxivity [s^-1 mM^-1] (concentration mode only).
        r2 (float or dict): transverse relaxivity [s^-1 mM^-1] (concentration mode only).

    Returns:
        dict: the arrays written to ``phantom_file`` (times in s, FOV in m).
    """
    if plot:
        # Fail early (before the computation) if the dynamic viewer is unavailable
        require_dynamic_voxel_phantom()

    phantom_data = load_phantom_parameters(param_json)
    T1, T2, times_post = get_T1_T2_over_time(
        concentrations_file, field_strength, times_post, r1=r1, r2=r2, param_json=param_json)
    time_frames = len(times_post)

    tubes_ID = get_tubes_ID_from_B0(field_strength)
    mask, FOV = build_label_mask(FOV, matrix, rotate_LR)

    # Time-varying T1/T2 maps [ms], shape (time_frames, x, y, z)
    T1_map = np.zeros((time_frames, *mask.shape))
    T2_map = np.zeros((time_frames, *mask.shape))
    for label, tube_ID in [*enumerate(tubes_ID, start=1), (OUTER_MATRIX_LABEL, 'OuterMatrix')]:
        T1_map[:, mask == label] = T1[tube_ID][:, None]
        T2_map[:, mask == label] = T2[tube_ID][:, None]

    # Static maps
    T2dash_map, PD_map, chi_map, D_map = (
        fill_static_map(mask, tubes_ID, phantom_data, key)
        for key in ('T2dash', 'PD', 'chi', 'D')
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
        time_points=np.asarray(times_post, dtype=float),
        **build_tissue_masks(mask, tubes_ID, phantom_data),
    )
    np.savez_compressed(phantom_file, **maps)
    print(f"Dynamic phantom ({time_frames} frames) saved to {phantom_file}")

    if plot:
        plot_phantom(phantom_file, plot_kwargs, dynamic=True)
    return maps


#: Backward-compatible name of :func:`build_dynamic_phantom`.
get_phantom_values = build_dynamic_phantom


def main():
    parser = argparse.ArgumentParser(
        description="Build a dynamic (contrast-injection) T1MES phantom (.npz) for the MRzero simulator.")
    dynamic_group = parser.add_argument_group(
        "Dynamic phantom parameters (give exactly one of -c or --times_post)")
    source = dynamic_group.add_mutually_exclusive_group(required=True)
    source.add_argument('-c', '--concentrations_file', type=str,
                        help='.npz file with gadolinium concentrations over time [mM] '
                             '(output of t1mes_simulate_fermi)')
    source.add_argument('--times_post', type=str_to_seconds, nargs='+', metavar='M[:S]',
                        help='Time after injection of each frame, in minutes (M) or minutes:seconds (M:S). '
                             'The T1/T2 time series are read from --param_json')
    dynamic_group.add_argument('--r1', type=float, default=DEFAULT_R1,
                               help=f'Longitudinal relaxivity [s^-1 mM^-1] (with -c). Default: {DEFAULT_R1}')
    dynamic_group.add_argument('--r2', type=float, default=DEFAULT_R2,
                               help=f'Transverse relaxivity [s^-1 mM^-1] (with -c). Default: {DEFAULT_R2}')
    add_phantom_arguments(parser, JSON_PARAMETERS, default_output_file='T1MES_dynamic_default.npz')
    args = finalize_phantom_arguments(parser, parser.parse_args())

    build_dynamic_phantom(
        concentrations_file=args.concentrations_file,
        times_post=args.times_post,
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
        r1=args.r1,
        r2=args.r2,
    )


if __name__ == "__main__":
    main()
