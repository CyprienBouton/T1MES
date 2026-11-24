"""Command-line helpers shared by the ``t1mes_*`` commands."""
from ast import literal_eval
import argparse
import json
import os

# Defaults used when the imaging geometry is not fully specified
DEFAULT_FOV = (200., 200., 50.)      # [mm]
DEFAULT_RESOLUTION = (1., 1., 5.)    # voxel size [mm]
DEFAULT_MATRIX = (200, 200, 10)      # number of voxels


# ---------------------------------------------------------------------------
# Imaging geometry: FOV = resolution * matrix
# ---------------------------------------------------------------------------
def is_close(a, b, tol=1e-3):
    """Return True if two sequences are element-wise equal within ``tol``."""
    return all(abs(x - y) < tol for x, y in zip(a, b))


def compute_fov(res, mat):
    """Field of view [mm] from voxel size [mm] and matrix size."""
    return tuple(float(r * m) for r, m in zip(res, mat))


def compute_resolution(fov, mat):
    """Voxel size [mm] from field of view [mm] and matrix size."""
    return tuple(float(f / m) for f, m in zip(fov, mat))


def compute_matrix(fov, res):
    """Matrix size from field of view [mm] and voxel size [mm] (rounded)."""
    return tuple(int(round(f / r)) for f, r in zip(fov, res))


def complete_imaging_args(args):
    """Fill in the missing values among ``args.FOV``, ``args.resolution`` and
    ``args.matrix`` so that ``FOV = resolution * matrix``.

    - all three given: they are checked for consistency;
    - two given: the third one is computed;
    - only one given: ``DEFAULT_MATRIX`` (or ``DEFAULT_FOV`` if only the matrix
      is given) is used to compute the remaining one;
    - none given: the defaults (200x200x50 mm, 1x1x5 mm, 200x200x10) are used.

    Raises:
        ValueError: if the three values are given but inconsistent.
    """
    FOV, res, mat = args.FOV, args.resolution, args.matrix

    if FOV is None and res is None and mat is None:
        FOV, res, mat = DEFAULT_FOV, DEFAULT_RESOLUTION, DEFAULT_MATRIX
    elif FOV is not None and res is None and mat is None:
        mat = DEFAULT_MATRIX
    elif FOV is None and res is not None and mat is None:
        mat = DEFAULT_MATRIX
    elif FOV is None and res is None and mat is not None:
        FOV = DEFAULT_FOV

    if FOV is None:
        FOV = compute_fov(res, mat)
    elif res is None:
        res = compute_resolution(FOV, mat)
    elif mat is None:
        mat = compute_matrix(FOV, res)

    expected_fov = compute_fov(res, mat)
    if not is_close(FOV, expected_fov):
        raise ValueError(
            f"Inconsistent geometry: FOV {tuple(FOV)} mm != resolution {tuple(res)} mm "
            f"x matrix {tuple(mat)} = {expected_fov} mm"
        )

    args.FOV = tuple(float(f) for f in FOV)
    args.resolution = tuple(float(r) for r in res)
    args.matrix = tuple(int(m) for m in mat)
    return args


# ---------------------------------------------------------------------------
# Parsing of free-form keyword arguments
# ---------------------------------------------------------------------------
def parse_key_value_or_json_file(arg):
    """argparse ``type`` accepting either a JSON file path or ``key=value``.

    ``value`` is parsed with :func:`ast.literal_eval`, so numbers, tuples,
    lists, booleans and quoted strings are supported, e.g.
    ``b0range=(-300,300)`` or ``title='T1MES'``.

    Returns:
        dict: the content of the JSON file, or ``{key: value}``.
    """
    if os.path.isfile(arg):
        try:
            with open(arg, 'r') as f:
                return json.load(f)
        except Exception as e:
            raise argparse.ArgumentTypeError(f"Could not parse JSON file '{arg}': {e}")

    try:
        key, value = arg.split('=', 1)
        return {key: literal_eval(value)}
    except Exception:
        raise argparse.ArgumentTypeError(
            f"Each argument must be key=value or a valid JSON file path. Got '{arg}'"
        )


def merge_dicts(dicts):
    """Merge a list of dicts (as produced by ``nargs='+'``) into a single dict."""
    if not dicts:
        return {}
    if isinstance(dicts, dict):
        return dict(dicts)
    return {k: v for d in dicts for k, v in d.items()}


def str_to_seconds(time_str: str) -> float:
    """argparse ``type`` converting a time given in minutes to seconds.

    Accepted formats: ``M`` (e.g. ``"16"`` or ``"2.5"``) or ``M:S``
    (e.g. ``"16:30"``).
    """
    try:
        if ':' in time_str:
            minutes, seconds = map(float, time_str.split(':'))
            return minutes * 60 + seconds
        return float(time_str) * 60
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid time format: '{time_str}'. Use M or M:S (minutes[:seconds])")


# ---------------------------------------------------------------------------
# Tube layout
# ---------------------------------------------------------------------------
def get_tubes_ID_from_B0(field_strength):
    """Return the tube identifiers placed at labels 1..9 of the mask.

    The T1MES phantom exists in a 1.5 T and a 3 T version. Both share the
    post-contrast tubes (A, C, G, H, I) and differ by their native tubes
    (B, D, E, F at 1.5 T; J, K, L, M at 3 T). The list is ordered row by row,
    from the top-left tube, i.e. ``tubes_ID[label - 1]`` is the tube at
    ``label`` in :func:`t1mes.generate_mask.create_t1mes_phantom_mask`.

    Args:
        field_strength (float): 1.5 or 3 [T].

    Raises:
        ValueError: for any other field strength.
    """
    if field_strength == 1.5:
        return [
            'I', 'F', 'C',
            'H', 'E', 'B',
            'G', 'D', 'A',
        ]
    if field_strength == 3.0:
        return [
            'I', 'M', 'C',
            'H', 'L', 'J',
            'G', 'K', 'A',
        ]
    raise ValueError(f'field strength must be 1.5 or 3 T, got {field_strength}')


# ---------------------------------------------------------------------------
# Arguments shared by t1mes_build_static and t1mes_build_dynamic
# ---------------------------------------------------------------------------
def add_phantom_arguments(parser, default_param_json, default_output_file='T1MES_default.npz'):
    """Add the imaging / configuration / output / plot arguments shared by the
    static and dynamic phantom builders to ``parser``."""
    imaging_group = parser.add_argument_group("Imaging parameters")
    imaging_group.add_argument('-FOV', '--FOV', dest='FOV', type=float, nargs=3, metavar=('X', 'Y', 'Z'),
                               help=f'Field of view [mm]. Default: {DEFAULT_FOV}')
    imaging_group.add_argument('-r', '--resolution', type=float, nargs=3, metavar=('X', 'Y', 'Z'),
                               help=f'Voxel size [mm]. Default: {DEFAULT_RESOLUTION}')
    imaging_group.add_argument('-m', '--matrix', type=int, nargs=3, metavar=('X', 'Y', 'Z'),
                               help=f'Matrix size (number of voxels). Default: {DEFAULT_MATRIX}')
    imaging_group.add_argument('-B0', '--field_strength', type=float, default=1.5, choices=[1.5, 3.0],
                               help='Main magnetic field strength [T], selects the 1.5 T or 3 T phantom. Default: 1.5')
    imaging_group.add_argument('--ncoils', type=int, default=8,
                               help='Number of receive coils of the sensitivity maps. Default: 8')
    imaging_group.add_argument('--rotate_LR', '--LR', dest='rotate_LR', action=argparse.BooleanOptionalAction,
                               default=False, help='Rotate the phantom by 90 degrees in the xy-plane')

    config_group = parser.add_argument_group("Phantom configuration")
    config_group.add_argument('--param_json', type=str, default=default_param_json,
                              help='JSON file with the tissue parameters of each tube. Default: packaged T1MES_values.json')
    config_group.add_argument('--b0_kwargs', type=parse_key_value_or_json_file, nargs='+', default={},
                              metavar='KEY=VALUE|FILE.json', help='Extra keyword arguments for mrtwin.b0field()')
    config_group.add_argument('--b1_kwargs', type=parse_key_value_or_json_file, nargs='+', default={},
                              metavar='KEY=VALUE|FILE.json', help='Extra keyword arguments for mrtwin.b1field()')
    config_group.add_argument('--sensmap_kwargs', type=parse_key_value_or_json_file, nargs='+', default={},
                              metavar='KEY=VALUE|FILE.json', help='Extra keyword arguments for mrtwin.sensmap()')

    output_group = parser.add_argument_group("Output")
    output_group.add_argument('-o', '--output_file', '-p', '--phantom_file', dest='output_file', type=str,
                              default=default_output_file,
                              help=f'Output phantom file (.npz). Default: {default_output_file}')

    vis_group = parser.add_argument_group("Visualization")
    vis_group.add_argument('--plot', action=argparse.BooleanOptionalAction, default=True,
                           help='Plot the resulting phantom with MRzeroCore')
    vis_group.add_argument('--plot_kwargs', type=parse_key_value_or_json_file, nargs='+', default={},
                           metavar='KEY=VALUE|FILE.json', help='Extra keyword arguments for the MRzeroCore plot method')
    return parser


def finalize_phantom_arguments(parser, args):
    """Complete the geometry and merge the ``*_kwargs`` lists into dicts."""
    try:
        args = complete_imaging_args(args)
    except ValueError as e:
        parser.error(str(e))
    for attr in ['b0_kwargs', 'b1_kwargs', 'sensmap_kwargs', 'plot_kwargs']:
        setattr(args, attr, merge_dicts(getattr(args, attr)))
    print(f"FOV: {args.FOV} mm | resolution: {args.resolution} mm | matrix: {args.matrix}")
    return args
