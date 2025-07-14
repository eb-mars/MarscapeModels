from landlab import RasterModelGrid
import numpy as np
import rasterio as rio
from rasterio.enums import Resampling
from rasterio.warp import reproject

# --- Functions to create various initial topographies with elevation patterns
def create_noisy_landscape(rows, cols, spacing=1000):
    """Creates a grid with simple random noise.
    
    Parameters:
    - rows (int): Number of rows in the grid.
    - cols (int): Number of columns in the grid.
    - spacing (float): Spacing between nodes in meters.

    Returns:
    - grid (RasterModelGrid): A Landlab grid with random noise added to the elevation
    """
    grid = RasterModelGrid((rows, cols), xy_spacing=spacing)
    z = grid.add_zeros('topographic__elevation', at='node')
    z += np.random.rand(grid.number_of_nodes) * 100  # Add 0-100m of noise
    return grid

def create_tilted_landscape(rows, cols, spacing=1000, tilt=0.01):
    """Creates a uniformly tilted grid.
    
    Parameters:
    - rows (int): Number of rows in the grid.
    - cols (int): Number of columns in the grid.
    - spacing (float): Spacing between nodes in meters.
    - tilt (float): Tilt factor for the landscape.
    
    Returns:
    - grid (RasterModelGrid): A Landlab grid with a tilted landscape
    """
    grid = RasterModelGrid((rows, cols), xy_spacing=spacing)
    z = grid.add_zeros('topographic__elevation', at='node')
    z += grid.x_of_node * tilt  # Tilt the landscape along the x-axis
    return grid

def create_step_landscape(rows, cols, step_height=100.0, spacing=1000):
    """
    Creates a grid with a single vertical step fault. The right half of the grid is elevated by the step_height.

    Parameters:
    - rows (int): Number of rows in the grid.
    - cols (int): Number of columns in the grid.
    - step_height (float): Height of the step in meters.
    - spacing (float): Spacing between nodes in meters.

    Returns:
    - grid (RasterModelGrid): A Landlab grid with a step landscape.
    """
    grid = RasterModelGrid((rows, cols), xy_spacing=spacing)
    z = grid.add_zeros('topographic__elevation', at='node')

    # Find the middle column of the grid
    middle_column_index = grid.number_of_node_columns // 2

    # Get all nodes on the right side of the grid (from the middle onward)
    right_side_nodes = grid.nodes[:, middle_column_index:]

    # Set the elevation of those nodes to the step height
    z[right_side_nodes] = step_height
    
    return grid

# --- Functions to add lithology patterns
def add_checkerboard_lithology(grid, erodibility_hard=1e-6, erodibility_soft=5e-6):
    """
    Adds a field for erodibility ('K_sp') in a checkerboard pattern.
    This function MODIFIES the grid it receives.

    Parameters:
    - grid (RasterModelGrid): The grid to modify.
    - erodibility_hard (float): Erodibility value for hard rock.
    - erodibility_soft (float): Erodibility value for soft rock.

    Returns:
    - grid (RasterModelGrid): The modified grid with the 'K_sp' field added
    """
    # Create the erodibility field, defaulting to the 'hard' value
    K_sp = grid.add_ones('K_sp', at='node') * erodibility_hard
    
    # Set every other 'column' of nodes to the 'soft' value
    # This is a simplified example
    for i in range(grid.number_of_node_columns):
        if i % 2 == 0:
            nodes_in_col = grid.nodes[:, i]
            K_sp[nodes_in_col] = erodibility_soft
            
    return grid # Return the modified grid

def import_topography(filename, cell_size):
    """
    Imports a topography from a GeoTIFF file and returns a grid at the same resolution as the synthetic landscapes
    
    Parameters:
    - filename (str): Path to the GeoTIFF file.
    - cell_size (float): Desired resolution of the grid in meters.
    
    Returns:
    - grid (RasterModelGrid): A Landlab grid with the imported topography.
    """

    # Read source DEM
    with rio.open(filename) as src:
        dem = src.read(1)
        dem[dem <= 0] = np.nan  
        src_transform = src.transform
        src_crs = src.crs
        bounds = src.bounds

    # Compute target shape from bounds and desired resolution
    new_width = int((bounds.right - bounds.left) / cell_size)
    new_height = int((bounds.top - bounds.bottom) / cell_size)

    # Build the new transform
    from rasterio.transform import from_origin
    new_transform = from_origin(bounds.left, bounds.top, cell_size, cell_size)

    # Destination array
    resampled = np.empty((new_height, new_width), dtype=np.float32)

    # Reproject and resample
    reproject(
        source=dem,
        destination=resampled,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=new_transform,
        dst_crs=src_crs,  # no CRS change, just resampling
        resampling=Resampling.bilinear  
    )
    grid = RasterModelGrid((new_height, new_width), xy_spacing=cell_size, at='node')
    return grid
