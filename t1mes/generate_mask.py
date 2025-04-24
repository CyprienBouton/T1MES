import numpy as np

# ---------------------------------------------------------------------
# Function 1: Create T1MES phantom mask
# ---------------------------------------------------------------------
def create_t1mes_phantom_mask(
    FOV=(200, 200, 200),
    resolution=(128,128,5),
    resin_size=(122, 122, 187),
    resin_radius=10,
    tube_diameter=28,
    tube_spacing=33,
):
    """ Generate 3D labeled mask of the T1MES phantom.
    
    Args:
       FOV (list, tuple): field of view of the phantom [mm].
       resolution (list, tuple): Resolution matrix.
       resin_size (list, tuple) Size of the Phantom [mm].
       resin_radius (float): Radius used for the outer matrix corner [mm]. 
       tube diameter (float): diameter of the tubes [mm].
       tube spacing (float): Spacing between the center of neighbor tubes.
    
    Returns:
        numpy array: Mask with the following labels:
            - 0: background
            - 10: resin body
            - 1–9: individual tubes
    """
    # Convert physical sizes to voxel units
    voxel_size = np.array(FOV)/resolution
    volume_vox = tuple(np.round(np.array(FOV) / voxel_size).astype(int))
    resin_vox = np.array(resin_size, dtype=float) / voxel_size
    radius_vox = resin_radius / voxel_size[0]  # Assume isotropic in XY
    tube_radius_vox = tube_diameter / 2 / voxel_size[0]
    spacing_vox = tube_spacing / voxel_size[0]

    # Create empty volume
    mask = np.zeros(volume_vox, dtype=np.uint8)
    resin_center = np.array(volume_vox) / 2
    
    # Bounds of resin body
    x0, y0, z0 = (resin_center - resin_vox / 2).astype(int)
    x1, y1, z1 = (resin_center + resin_vox / 2).astype(int)

    # Build 2D XY mask for rounded rectangle
    yy, xx = np.ogrid[0:volume_vox[1], 0:volume_vox[0]]
    resin_xy = np.zeros((volume_vox[1], volume_vox[0]), dtype=bool)

    core_mask = (
        ((xx >= x0 + radius_vox) & (xx < x1 - radius_vox) & (yy >= y0) & (yy < y1)) |
        ((xx >= x0) & (xx < x1) & (yy >= y0 + radius_vox) & (yy < y1 - radius_vox))
    )
    resin_xy |= core_mask

    corners = [
        (x0 + radius_vox, y0 + radius_vox),
        (x1 - radius_vox - 1, y0 + radius_vox),
        (x0 + radius_vox, y1 - radius_vox - 1),
        (x1 - radius_vox - 1, y1 - radius_vox - 1)
    ]
    for cx, cy in corners:
        circle = (xx - cx)**2 + (yy - cy)**2 <= radius_vox**2
        resin_xy |= circle

    # Extrude 2D shape across Z
    mask[resin_xy, z0:z1+1] = 10

    # Create 9 tubes (labels 1 to 9)
    tube_radius = tube_diameter / 2

    xx, yy, zz = np.meshgrid(
        np.arange(volume_vox[0]),
        np.arange(volume_vox[1]),
        np.arange(volume_vox[2]),
        indexing="ij"
    )

    # Tube centers (3x3 grid)
    grid_origin = resin_center[:2] - spacing_vox 
    tube_idx = 1
    for i in range(3):
        for j in range(3):
            cx = grid_origin[0] + j * spacing_vox
            # Invert row index so tube grid is ordered top to bottom
            cy = grid_origin[1] + (2-i) * spacing_vox
            r2 = (xx - cx)**2 + (yy - cy)**2
            cylinder = r2 <= tube_radius_vox**2
            mask[..., cylinder] = tube_idx
            tube_idx += 1
    return mask