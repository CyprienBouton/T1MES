"""Geometry of the T1MES phantom: a labelled 3D mask of the tubes and resin body.

The T1MES phantom (Captur et al., JCMR 2016, doi:10.1186/s12968-016-0280-z)
is a rounded-rectangle block of gel ("outer matrix") holding nine cylindrical
tubes on a 3x3 grid. The dimensions used here are approximations of the
physical device and can be changed through the keyword arguments.
"""
import numpy as np

# Labels used in the mask
BACKGROUND_LABEL = 0
OUTER_MATRIX_LABEL = 10
TUBE_LABELS = tuple(range(1, 10))  # 1..9, row by row from the top-left tube


def create_t1mes_phantom_mask(
    FOV=(200, 200, 200),
    matrix=(128, 128, 5),
    resin_size=(122, 122, 187),
    resin_radius=10,
    tube_diameter=28,
    tube_spacing=33,
):
    """Generate a 3D labelled mask of the T1MES phantom.

    The phantom is centred in the field of view. The tubes are parallel to the
    third (z) axis and are restricted to the z-extent of the resin body.

    Args:
        FOV (sequence of 3 float): field of view [mm].
        matrix (sequence of 3 int): number of voxels along each axis.
        resin_size (sequence of 3 float): outer size of the resin body [mm].
        resin_radius (float): radius of the rounded corners of the resin body
            in the xy-plane [mm].
        tube_diameter (float): diameter of each tube [mm].
        tube_spacing (float): centre-to-centre distance between neighbouring
            tubes [mm].

    Returns:
        np.ndarray: uint8 array of shape ``matrix`` with labels
            - 0: background,
            - 1-9: tubes, numbered row by row starting from the top-left tube
              (see ``t1mes.convert_arguments.get_tubes_ID_from_B0`` for the
              tube placed at each label),
            - 10: resin body (outer matrix).

    Note:
        Voxels are assumed square in the xy-plane: in-plane distances
        (corner radius, tube diameter and spacing) are converted to voxels
        with the voxel size along x.
    """
    matrix = tuple(int(n) for n in matrix)
    voxel_size = np.asarray(FOV, dtype=float) / np.asarray(matrix, dtype=float)
    resin_vox = np.asarray(resin_size, dtype=float) / voxel_size
    radius_vox = resin_radius / voxel_size[0]
    tube_radius_vox = tube_diameter / 2 / voxel_size[0]
    spacing_vox = tube_spacing / voxel_size[0]

    mask = np.zeros(matrix, dtype=np.uint8)
    center = np.asarray(matrix, dtype=float) / 2

    # Bounds of the resin body in voxel units (clipped to the volume)
    x0, y0, z0 = (center - resin_vox / 2).astype(int)
    x1, y1, z1 = (center + resin_vox / 2).astype(int)
    z0, z1 = max(z0, 0), min(z1, matrix[2] - 1)

    # In-plane coordinates, indexed as (x, y) to match the mask layout
    xx, yy = np.meshgrid(np.arange(matrix[0]), np.arange(matrix[1]), indexing="ij")

    # Rounded rectangle = two overlapping rectangles + four corner discs
    resin_xy = (
        ((xx >= x0 + radius_vox) & (xx < x1 - radius_vox) & (yy >= y0) & (yy < y1))
        | ((xx >= x0) & (xx < x1) & (yy >= y0 + radius_vox) & (yy < y1 - radius_vox))
    )
    corners = [
        (x0 + radius_vox, y0 + radius_vox),
        (x1 - radius_vox - 1, y0 + radius_vox),
        (x0 + radius_vox, y1 - radius_vox - 1),
        (x1 - radius_vox - 1, y1 - radius_vox - 1),
    ]
    for cx, cy in corners:
        resin_xy |= (xx - cx) ** 2 + (yy - cy) ** 2 <= radius_vox ** 2

    # Extrude the 2D shape along z
    mask[resin_xy, z0:z1 + 1] = OUTER_MATRIX_LABEL

    # Nine tubes (labels 1..9) on a 3x3 grid centred in the resin body
    grid_origin = center[:2] - spacing_vox
    label = 1
    for i in range(3):
        for j in range(3):
            cx = grid_origin[0] + j * spacing_vox
            # Invert the row index so that tubes are numbered from top to bottom
            cy = grid_origin[1] + (2 - i) * spacing_vox
            tube_xy = (xx - cx) ** 2 + (yy - cy) ** 2 <= tube_radius_vox ** 2
            mask[tube_xy, z0:z1 + 1] = label
            label += 1
    return mask
