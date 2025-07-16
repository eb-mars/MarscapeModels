from landlab import RasterModelGrid
import numpy as np

def add_uniform_precip(grid, precipitation_rate=1.0):
    """
    Adds a uniform precipitation field to the grid.

    This function adds the field 'water__unit_flux_in' and sets its value
    to the precipitation rate only for nodes where 'topographic__elevation'
    is not NaN.

    Parameters:
        grid (RasterModelGrid): The Landlab grid to modify.
        precipitation_rate (float): Precipitation rate (e.g., in m/year).

    Returns:
        grid (RasterModelGrid): The modified grid.
    """
    # Ensure the 'water__unit_flux_in' field exists, initializing to zero.
    # Using add_zeros is robust; it creates the field if it doesn't exist.
    if 'water__unit_flux_in' not in grid.at_node:
        grid.add_zeros('node', 'water__unit_flux_in')
    
    # Set the entire field to zero first to clear previous values.
    grid.at_node['water__unit_flux_in'][:] = 0.0

    # Find the node IDs where topographic elevation is a finite number (i.e., not NaN).
    valid_nodes = np.where(np.isfinite(grid.at_node['topographic__elevation']))[0]
    
    # Assign the precipitation rate to these valid nodes.
    grid.at_node['water__unit_flux_in'][valid_nodes] = precipitation_rate
    
    return grid

def add_random_precip(grid, min_rate=0.5, max_rate=2.0):
    """
    Adds a random precipitation field to the grid.

    This function adds the field 'water__unit_flux_in' and sets its value
    to a random precipitation rate between min_rate and max_rate for nodes
    where 'topographic__elevation' is not NaN.

    Parameters:
        grid (RasterModelGrid): The Landlab grid to modify.
        min_rate (float): Minimum precipitation rate (e.g., in m/year).
        max_rate (float): Maximum precipitation rate (e.g., in m/year).

    Returns:
        grid (RasterModelGrid): The modified grid.
    """
    # Ensure the 'water__unit_flux_in' field exists, initializing to zero.
    if 'water__unit_flux_in' not in grid.at_node:
        grid.add_zeros('node', 'water__unit_flux_in')
    
    # Set the entire field to zero first to clear previous values.
    grid.at_node['water__unit_flux_in'][:] = 0.0

    # Find the node IDs where topographic elevation is a finite number (i.e., not NaN).
    valid_nodes = np.where(np.isfinite(grid.at_node['topographic__elevation']))[0]
    
    # Generate random precipitation rates for these valid nodes.
    random_rates = np.random.uniform(min_rate, max_rate, size=len(valid_nodes))
    
    # Assign the random precipitation rates to these valid nodes.
    grid.at_node['water__unit_flux_in'][valid_nodes] = random_rates
    
    return grid
