from landlab import RasterModelGrid
import numpy as np
import rasterio as rio
from rasterio.enums import Resampling
from rasterio.warp import reproject

# --- Functions to create various initial topographies with elevation patterns
def create_noisy_landscape(rows, cols, cell_size=1000, rf=0.1, seed=3):
    """Creates a grid with simple random noise.
    
    Parameters:
        rows (int): Number of rows in the grid.
        cols (int): Number of columns in the grid.
        cell_size (float): Spacing between nodes in meters.
        rf (float): Random factor to add noise to the landscape, as a percentage  of the cell size.
        seed (int): Random seed for reproducibility.

    Returns:
        grid (RasterModelGrid): A Landlab grid with random noise added to the elevation
    """
    grid = RasterModelGrid((rows, cols), xy_spacing=cell_size)
    z = grid.add_zeros('topographic__elevation', at='node')
    np.random.seed(seed)
    z += np.random.rand(grid.number_of_nodes) * rf * cell_size
    return grid

def create_tilted_landscape(rows, cols, cell_size=1000, tilt=0.01, rf=0.1, tilt_direction='South', seed=3):
    """Creates a uniformly tilted grid.
    
    Parameters:
        rows (int): Number of rows in the grid.
        cols (int): Number of columns in the grid.
        cell_size (float): Spacing between nodes in meters.
        rf (float): Random factor to add noise to the landscape, as a percentage  of the cell size.
        tilt (float): Tilt factor for the landscape.
        tilt_direction (str): Direction of the tilt ('North', 'South', 'East', 'West'). e.g., 'South' means the landscape slopes down towards the south.
        seed (int): Random seed for reproducibility.
        
    Returns:
        grid (RasterModelGrid): A Landlab grid with a tilted landscape
    """
    grid = RasterModelGrid((rows, cols), xy_spacing=cell_size)
    z = grid.add_zeros('topographic__elevation', at='node')
    if tilt_direction == 'North':
        z -= grid.y_of_node * tilt  # Tilt the landscape along the y-axis
    elif tilt_direction == 'South':
        z += grid.y_of_node * tilt  # Tilt the landscape along the y-axis
    elif tilt_direction == 'East':
        z -= grid.x_of_node * tilt  # Tilt the landscape along the x-axis
    elif tilt_direction == 'West':
        z += grid.x_of_node * tilt  # Tilt the landscape along the x-axis
    else:
        raise ValueError("Invalid tilt direction. Choose from 'North', 'South', 'East', or 'West'.")
    np.random.seed(seed)
    z += np.random.rand(grid.number_of_nodes) * rf * cell_size # Add some noise
    return grid

def create_step_landscape(rows, cols, cell_size=1000, step_height=100.0, rf=0.1, seed=3):
    """
    Creates a grid with a single vertical step fault. The right half of the grid is elevated by the step_height.

    Parameters:
        rows (int): Number of rows in the grid.
        cols (int): Number of columns in the grid.
        cell_size (float): Spacing between nodes in meters.
        step_height (float): Height of the step in meters.
        rf (float): Random factor to add noise to the landscape, as a percentage of the cell size.
        seed (int): Random seed for reproducibility.
        
    Returns:
        grid (RasterModelGrid): A Landlab grid with a step landscape.
    """
    grid = RasterModelGrid((rows, cols), xy_spacing=cell_size)
    z = grid.add_zeros('topographic__elevation', at='node')

    # Find the middle column of the grid
    middle_column_index = grid.number_of_node_columns // 2

    # Get all nodes on the right side of the grid (from the middle onward)
    right_side_nodes = grid.nodes[:, middle_column_index:]

    # Set the elevation of those nodes to the step height
    z[right_side_nodes] = step_height

    np.random.seed(seed)
    z += np.random.rand(grid.number_of_nodes) * rf * cell_size  # Add some noise
    
    return grid

# --- Functions to add lithology patterns
def add_uniform_lithology(grid, erodibility=1e-6):
    """
    Adds a uniform lithology field to the grid.

    Parameters:
        grid (RasterModelGrid): The grid to modify.
        erodibility (float): Erodibility value for the lithology.

    Returns:
        grid (RasterModelGrid): The modified grid with the 'K_sp' field added
    """
    # Create a uniform erodibility field
    K_sp = grid.add_ones('K_sp', at='node')
    K_sp *= erodibility
    return grid  # Return the modified grid

def add_checkerboard_lithology(grid, erodibility_1=1e-6, erodibility_2=5e-6):
    """
    Adds a field for erodibility ('K_sp') in a checkerboard pattern.
    This function MODIFIES the grid it receives.

    Parameters:
        grid (RasterModelGrid): The grid to modify.
        erodibility_hard (float): Erodibility value for hard rock.
        erodibility_soft (float): Erodibility value for soft rock.

    Returns:
        grid (RasterModelGrid): The modified grid with the 'K_sp' field added
    """
    K_sp = grid.add_ones('K_sp', at='node')
    K_sp *= erodibility_1 
    
    row_index = (grid.y_of_node // grid.dy).astype(int)
    col_index = (grid.x_of_node // grid.dx).astype(int)
    
    is_even_patch = (row_index + col_index) % 2 == 0
    
    K_sp[is_even_patch] = erodibility_2
            
    return grid

def add_striped_lithology(grid, layer_thickness=200.0, erodibility_1=1e-6, erodibility_2=5e-6):
    """
    Adds horizontal layers of two rock types to the grid.

    This function creates a 'K_sp' field on the grid, alternating between
    two erodibility values based on the y-coordinate.

    Parameters:
        grid (RasterModelGrid): The Landlab grid to modify.
        layer_thickness (float): The thickness of each rock layer in meters.
        K_rock_1 (float): Erodibility of the first rock type.
        K_rock_2 (float): Erodibility of the second rock type.

    Returns:
        The modified grid with the 'K_sp' field.
    """
    # Create the erodibility field, defaulting to the first rock type
    K_sp = grid.add_ones('K_sp', at='node') 
    K_sp *= erodibility_1
    
    # Calculate which layer each node belongs to by dividing its y-coordinate by the thickness
    # The modulo operator (%) determines if a layer is even or odd
    layer_number = (grid.y_of_node // layer_thickness) % 2
    
    # Assign the second rock type to all nodes in the 'odd' numbered layers (layer_number == 1)
    K_sp[layer_number == 1] = erodibility_2
    
    return grid

def add_split_lithology(grid, erodibility_left=1e-6, erodibility_right=5e-6):
    """
    Adds a vertically split lithology to the grid.

    This function creates a 'K_sp' field with one rock type on the left
    half of the grid and another on the right half.

    Parameters:
        grid (RasterModelGrid): The Landlab grid to modify.
        K_left (float): Erodibility of the left-side rock type.
        K_right (float): Erodibility of the right-side rock type.

    Returns:
        The modified grid with the 'K_sp' field.
    """
    # 1. Create the field, then get a reference to it to modify in-place.
    K_sp = grid.add_ones('K_sp', at='node')
    K_sp[:] = erodibility_left # Set the default value for all nodes

    # 2. Split the grid by column INDEX, not by coordinate.
    # This ensures an even split of node columns.
    mid_col_index = grid.number_of_node_columns // 2
    
    # Get the column index for every node.
    col_indices = (grid.x_of_node // grid.dx).astype(int)
    
    # Find nodes whose column index is in the right half.
    right_side_nodes = np.where(col_indices >= mid_col_index)[0]

    # Assign the right-side erodibility to those nodes.
    K_sp[right_side_nodes] = erodibility_right
    
    return grid

def import_topography(filename, cell_size):
    """
    Imports a topography from a GeoTIFF file and returns a grid at the same resolution as the synthetic landscapes
    
    Parameters:
        filename (str): Path to the GeoTIFF file.
        cell_size (float): Desired resolution of the grid in meters.
    
    Returns:
        grid (RasterModelGrid): A Landlab grid with the imported topography.
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
    grid = RasterModelGrid((new_height, new_width), xy_spacing=cell_size)
    resampled = np.flipud(resampled)
    grid.add_field("topographic__elevation", resampled.flatten(), at="node")

    return grid
