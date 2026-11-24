"""Building blocks shared by the static and dynamic T1MES phantom builders."""
import importlib.resources as pkg_resources
import inspect
import json

import numpy as np

from t1mes.generate_mask import OUTER_MATRIX_LABEL, create_t1mes_phantom_mask

#: Default tissue parameters shipped with the package (see README, "JSON configuration").
JSON_PARAMETERS = pkg_resources.files("t1mes").joinpath("T1MES_values.json").as_posix()

#: Git requirement providing ``MRzeroCore.DynamicVoxelPhantom``.
MRZERO_DYNAMIC_REQUIREMENT = (
    "MRzeroCore @ git+https://github.com/CyprienBouton/MRzero-Core.git@T1_T2_variable_v1.6"
)

# Default B0 inhomogeneity range passed to mrtwin.b0field() [Hz]
DEFAULT_B0_RANGE = (-200, 200)


def load_phantom_parameters(param_json=JSON_PARAMETERS):
    """Load the tissue parameters of the T1MES tubes from a JSON file.

    Returns:
        dict: ``{tube_ID: {"description", "PD", "T1", "T2", "T2dash", "D", "chi", ...}}``
        with times in ms (see README for all fields and units).
    """
    with open(param_json, 'r') as f:
        return json.load(f)


def build_label_mask(FOV, matrix, rotate_LR=False):
    """Labelled T1MES mask, optionally rotated by 90 degrees in the xy-plane.

    Returns:
        tuple: (mask, FOV) where FOV has its x and y entries swapped if the
        phantom was rotated, so that it stays consistent with the mask shape.
    """
    mask = create_t1mes_phantom_mask(FOV, matrix)
    FOV = tuple(FOV)
    if rotate_LR:
        mask = np.rot90(mask, k=1, axes=(0, 1)).copy()
        FOV = (FOV[1], FOV[0], FOV[2])
    return mask, FOV


def fill_static_map(mask, tubes_ID, phantom_data, key):
    """Return a map with ``phantom_data[tube][key]`` in each tube and in the outer matrix."""
    out = np.zeros(mask.shape, dtype=float)
    for label, tube_ID in [*enumerate(tubes_ID, start=1), (OUTER_MATRIX_LABEL, 'OuterMatrix')]:
        value = phantom_data[tube_ID][key]
        if not np.isscalar(value):
            raise ValueError(
                f"'{key}' of tube '{tube_ID}' is a time series; "
                "use t1mes_build_dynamic with --times_post for such JSON files."
            )
        out[mask == label] = value
    return out


def compute_field_maps(PD_map, chi_map, field_strength, ncoils=8,
                       b0field_kwargs=None, b1field_kwargs=None, sensmap_kwargs=None):
    """Simulate B0, B1 and receive-coil sensitivity maps with ``mrtwin``.

    - B0: off-resonance [Hz] computed by ``mrtwin.b0field`` from the
      susceptibility map ``chi_map`` [ppm]; rescaled to ``b0range``
      (default -200..200 Hz).
    - B1: relative transmit field (1 = nominal flip angle) from ``mrtwin.b1field``.
    - coil_sens: ``ncoils`` complex receive sensitivities from ``mrtwin.sensmap``.

    The user kwargs dicts are copied, never modified in place.

    Returns:
        tuple: (B0_map, B1_map, coil_sens) with shapes ``PD_map.shape``,
        ``PD_map.shape`` (or ``(nmodes, *PD_map.shape)``) and ``(ncoils, *PD_map.shape)``.
    """
    from mrtwin import b0field, b1field, sensmap

    shape = PD_map.shape
    squeezed_shape = PD_map.squeeze().shape

    b0_kwargs = dict(b0field_kwargs or {})
    b0_kwargs['chi'] = chi_map
    b0_kwargs['B0'] = field_strength
    b0_kwargs.setdefault('b0range', DEFAULT_B0_RANGE)
    B0_map = np.asarray(b0field(**b0_kwargs))
    if B0_map.size == PD_map.size:
        B0_map = B0_map.reshape(shape)

    b1_kwargs = dict(b1field_kwargs or {})
    b1_kwargs['shape'] = squeezed_shape
    b1_kwargs['mask'] = PD_map.squeeze() > 0
    B1_map = np.asarray(b1field(**b1_kwargs))
    B1_map = B1_map.reshape(shape) if B1_map.size == PD_map.size else B1_map.reshape(-1, *shape)

    sens_kwargs = dict(sensmap_kwargs or {})
    sens_kwargs['shape'] = (ncoils, *squeezed_shape)
    coil_sens = np.asarray(sensmap(**sens_kwargs)).reshape(ncoils, *shape)

    return B0_map, B1_map, coil_sens


def build_tissue_masks(mask, tubes_ID, phantom_data):
    """Boolean mask of each tube and of the outer matrix, keyed ``tissue_<description>``.

    MRzeroCore's ``VoxelGridPhantom.load`` reads every ``tissue_*`` entry of
    the .npz file into ``phantom.tissue_masks``.
    """
    masks = {
        f"tissue_{phantom_data[tube_ID]['description']}": mask == label
        for label, tube_ID in enumerate(tubes_ID, start=1)
    }
    masks['tissue_OuterMatrix'] = mask == OUTER_MATRIX_LABEL
    return masks


def _filter_kwargs(func, kwargs):
    """Keep only the kwargs accepted by ``func`` (MRzeroCore plot signatures
    differ between versions, e.g. ``display_units`` only exists in the fork)."""
    params = inspect.signature(func).parameters
    if any(p.kind == p.VAR_KEYWORD for p in params.values()):
        return dict(kwargs)
    dropped = sorted(set(kwargs) - set(params))
    if dropped:
        print(f"Ignoring plot options not supported by this MRzeroCore version: {dropped}")
    return {k: v for k, v in kwargs.items() if k in params}


def plot_phantom(phantom_file, plot_kwargs=None, dynamic=False):
    """Load ``phantom_file`` with MRzeroCore and plot it.

    Static phantoms use ``VoxelGridPhantom.plot``; dynamic phantoms use
    ``DynamicVoxelPhantom.plot_dynamic`` (only available in the MRzeroCore fork).
    """
    import MRzeroCore as mr0

    kwargs = {'time_unit': 'ms', 'display_units': True}
    kwargs.update(plot_kwargs or {})
    if dynamic:
        phantom = require_dynamic_voxel_phantom().load(phantom_file)
        method = phantom.plot_dynamic
    else:
        phantom = mr0.VoxelGridPhantom.load(phantom_file)
        method = phantom.plot
    method(**_filter_kwargs(method, kwargs))


def require_dynamic_voxel_phantom():
    """Return ``MRzeroCore.DynamicVoxelPhantom`` or raise a helpful ImportError.

    ``DynamicVoxelPhantom`` (time-varying T1/T2 maps) is not part of the
    MRzeroCore release on PyPI; it is provided by a fork of MRzeroCore.
    """
    try:
        import MRzeroCore as mr0
    except ImportError as e:
        raise ImportError(
            "MRzeroCore is not installed. Install it with:\n"
            f"    pip install \"{MRZERO_DYNAMIC_REQUIREMENT}\""
        ) from e
    cls = getattr(mr0, 'DynamicVoxelPhantom', None)
    if cls is None:
        raise ImportError(
            "The installed MRzeroCore "
            f"({getattr(mr0, '__version__', 'unknown version')}) has no DynamicVoxelPhantom, "
            "which is needed to load, plot and simulate dynamic phantoms.\n"
            "Install the MRzeroCore fork with:\n"
            f"    pip install \"{MRZERO_DYNAMIC_REQUIREMENT}\"\n"
            "or install t1mes with the 'dynamic' extra: pip install \".[dynamic]\".\n"
            "Use --no-plot to only write the .npz file."
        )
    return cls
