# initial_conditions.py
from landlab import RasterModelGrid
import numpy as np

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
    - grid (RasterModelGrid): A Landlab grid with a tilted landscape"""
    grid = RasterModelGrid((rows, cols), xy_spacing=spacing)
    z = grid.add_zeros('topographic__elevation', at='node')
    z += grid.x_of_node * tilt  # Tilt the landscape along the x-axis
    return grid
    
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