import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Function 1: Create T1MES phantom mask
# ---------------------------------------------------------------------
def create_t1mes_phantom_mask(
    volume_size=(200, 200, 200),
    resin_size=(122, 122, 187),
    resin_radius=10,
    tube_diameter=28,
    tube_spacing=33
):
    """
    Generate a 3D labeled mask of the T1MES phantom.
    Returns a numpy array with:
        - label 0: background
        - label 10: resin body
        - labels 1–9: individual tubes
    """
    volume_shape = tuple(np.array(volume_size).astype(int))
    resin_size = np.array(resin_size)
    resin_center = np.array(volume_size) / 2
    mask = np.zeros(volume_shape, dtype=np.uint8)

    # Bounds of resin body
    x0, y0, z0 = (resin_center - resin_size / 2).astype(int)
    x1, y1, z1 = (resin_center + resin_size / 2).astype(int)

    # Build 2D XY mask for rounded rectangle
    yy, xx = np.ogrid[0:volume_shape[1], 0:volume_shape[0]]
    resin_xy = np.zeros((volume_shape[1], volume_shape[0]), dtype=bool)

    core_mask = (
        ((xx >= x0 + resin_radius) & (xx < x1 - resin_radius) & (yy >= y0) & (yy < y1)) |
        ((xx >= x0) & (xx < x1) & (yy >= y0 + resin_radius) & (yy < y1 - resin_radius))
    )
    resin_xy |= core_mask

    corners = [
        (x0 + resin_radius, y0 + resin_radius),
        (x1 - resin_radius - 1, y0 + resin_radius),
        (x0 + resin_radius, y1 - resin_radius - 1),
        (x1 - resin_radius - 1, y1 - resin_radius - 1)
    ]
    for cx, cy in corners:
        circle = (xx - cx)**2 + (yy - cy)**2 <= resin_radius**2
        resin_xy |= circle

    # Extrude 2D shape across Z
    mask[resin_xy, z0:z1] = 10

    # Create 9 tubes (labels 1 to 9)
    tube_radius = tube_diameter / 2

    xx, yy, zz = np.meshgrid(
        np.arange(volume_shape[0]),
        np.arange(volume_shape[1]),
        np.arange(volume_shape[2]),
        indexing="ij"
    )

    # Tube centers (3x3 grid)
    grid_origin = resin_center[:2] - tube_spacing 
    tube_idx = 1
    for i in range(3):
        for j in range(3):
            cx = grid_origin[0] + j * tube_spacing
            # Invert row index so tube grid is ordered top to bottom
            cy = grid_origin[1] + (2-i) * tube_spacing
            r2 = (xx - cx)**2 + (yy - cy)**2
            cylinder = r2 <= tube_radius**2
            mask[..., cylinder] = tube_idx
            tube_idx += 1
    return mask