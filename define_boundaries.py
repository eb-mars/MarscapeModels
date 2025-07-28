import numpy as np

def define_boundaries(grid, slope_direction):

    """Defines the boundaries of the grid.
    Sets all boundaries/edges to closed, except the one downstream of slope, 
    which is set as a fixed gradient boundary (BC_NODE_IS_FIXED_GRADIENT)
    
    Parameters:
        grid (RasterModelGrid): The grid to set boundaries for.
        slope_direction (str): Direction of the slope, which determines which boundary is open.
    """
    #Slope direction tells us which boundary should be open. False means open.
    North = "North" != slope_direction
    East = "East" != slope_direction
    South = "South" != slope_direction
    West = "West" != slope_direction

    grid.set_closed_boundaries_at_grid_edges(right_is_closed=East, top_is_closed=North, left_is_closed=West, bottom_is_closed=South)
    # self.grid.status_at_node[self.grid.fixed_value_boundary_nodes] = self.grid.BC_NODE_IS_FIXED_GRADIENT
    grid.status_at_node[grid.fixed_value_boundary_nodes] = grid.BC_NODE_IS_FIXED_VALUE
    return grid

def define_boundaries_outlet_fixed_value(grid):
    """Defines the boundaries of a single cell grid.
    
    Parameters:
        grid (RasterModelGrid): The grid to set boundaries for.
    """
    # Create a grid and the elevation field (with ties on the boundary)
    z = grid.at_node['topographic__elevation']

    # Automate the random selection of an outlet node from ties
    boundary_nodes = grid.boundary_nodes
    boundary_elevs = z[boundary_nodes]
    min_elev = np.min(boundary_elevs)

    # Find all boundary nodes that have the minimum elevation
    tied_nodes = boundary_nodes[boundary_elevs == min_elev]

    # Randomly choose one of these nodes to be the outlet
    random_outlet_node = np.random.choice(tied_nodes)

    # Set the watershed boundary using the randomly chosen ID
    grid.set_watershed_boundary_condition_outlet_id(random_outlet_node, z)
    return grid
        
def define_boundaries_outlet_fixed_gradient(grid, outlet_gradient=0.002):
    """
    Defines watershed boundaries with a single fixed-gradient outlet.
    
    This function selects the lowest boundary node as the outlet, sets all
    other boundaries to closed, and adjusts the outlet's elevation to ensure
    it connects properly to the watershed interior.

    Parameters:
        grid (RasterModelGrid): The grid to set boundaries for.
        outlet_gradient (float): The gradient to enforce at the outlet node.
                                Default is 0.002 (0.2%).
    """
    # 1. Select the lowest boundary node as the outlet
    z = grid.at_node['topographic__elevation']
    boundary_nodes = grid.boundary_nodes
    outlet_node = boundary_nodes[np.argmin(z[boundary_nodes])]

    # 2. Set all boundary nodes to CLOSED
    grid.status_at_node[grid.boundary_nodes] = grid.BC_NODE_IS_CLOSED

    # 3. Set the chosen outlet node status to FIXED_GRADIENT
    grid.status_at_node[outlet_node] = grid.BC_NODE_IS_FIXED_GRADIENT

    # Find the active (core) neighbor node of the outlet
    neighbors = grid.active_adjacent_nodes_at_node[outlet_node]
    if neighbors.size == 0:
        raise ValueError("Outlet node has no active neighbors. Check grid setup.")
    
    inland_neighbor = neighbors[0]

    # Get the link connecting the neighbor and the outlet
    link_to_outlet = np.intersect1d(grid.links_at_node[inland_neighbor], grid.links_at_node[outlet_node])[0]
    link_length = grid.length_of_link[link_to_outlet]

    # Calculate the new elevation for the outlet
    neighbor_elevation = grid.at_node['topographic__elevation'][inland_neighbor]
    new_outlet_elevation = neighbor_elevation - (outlet_gradient * link_length)

    # Set the new elevation in the grid field
    grid.at_node['topographic__elevation'][outlet_node] = new_outlet_elevation

    # 5. Set the gradient value for the outlet node for other components to use
    grid.at_node['topographic__steepest_slope'][outlet_node] = outlet_gradient

    return grid
    
def define_boundaries_outlet_center(grid, slope_direction):
    """Defines the boundaries of the grid.

    Sets all boundaries to closed, then sets the single center node
    on the specified edge to be a fixed-gradient outlet.

    Parameters:
        slope_direction (str): Direction of the slope ('North', 'East',
                            'South', or 'West'), which determines the
                            outlet location.
    """
    # Get grid dimensions for calculating the center node
    nrows, ncols = grid.shape

    # 1. Set all boundaries to closed initially
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)

    # 2. Find the outlet node ID at the center of the specified edge
    if slope_direction == "South":
        outlet_node = grid.nodes_at_bottom_edge[(ncols - 1) // 2]
    elif slope_direction == "North":
        outlet_node = grid.nodes_at_top_edge[(ncols - 1) // 2]
    elif slope_direction == "West":
        outlet_node = grid.nodes_at_left_edge[(nrows - 1) // 2]
    elif slope_direction == "East":
        outlet_node = grid.nodes_at_right_edge[(nrows - 1) // 2]
    else:
        raise ValueError("slope_direction must be one of: 'North', 'East', 'South', 'West'")

    # 3. Set the identified outlet node to a fixed-gradient boundary
    grid.at_node['topographic__elevation'][outlet_node] -= 5.
    grid.status_at_node[outlet_node] = grid.BC_NODE_IS_FIXED_VALUE
    return grid

def define_boundaries_corner(grid, slope_direction):
    nrows, ncols = grid.shape

    if slope_direction == "Northwest":
        # Top-left corner
        outlet_node = (nrows - 1) * ncols
    elif slope_direction == "Northeast":
        # Top-right corner (the very last node)
        outlet_node = grid.number_of_nodes - 1
    elif slope_direction == "Southwest":
        # Bottom-left corner
        outlet_node = 0
    elif slope_direction == "Southeast":
        # Bottom-right corner
        outlet_node = ncols - 1
    else:
        raise ValueError("slope_direction must be one of: 'Northwest', 'Northeast', 'Southwest', 'Southeast'")

    # Set the identified outlet node to a fixed-gradient boundary
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    grid.status_at_node[outlet_node] = grid.BC_NODE_IS_FIXED_GRADIENT
    grid.at_node['topographic__elevation'][outlet_node] -= 5.
    return grid
