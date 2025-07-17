import numpy as np

def get_channel_erosion_and_discharge(grid, profiler,
                                      initial_topo='initial_topographic__elevation',
                                      discharge_field='surface_water__discharge'):
    """
    Calculates erosional depth and discharge along channels.

    NOTE: The grid must have the discharge field calculated before calling this
    function. This is typically done by running the FlowAccumulator after
    setting a 'water__unit_flux_in' field.

    Parameters:
        grid (RasterModelGrid): The grid containing topography and discharge data.
        profiler (ChannelProfiler): A ChannelProfiler instance that has been run on the grid.
        initial_topo (str): Grid field name for the initial topography.
        discharge_field (str): Grid field name for the discharge.

    Returns:
        dict: A dictionary keyed by outlet node ID. Each value contains arrays for:
        'nodes', 'distance_from_outlet', 'erosion_depth', and 'discharge'.
    """
    # Check for required fields
    if initial_topo not in grid.at_node:
        raise ValueError(f"Field '{initial_topo}' not found.")
    if 'topographic__elevation' not in grid.at_node:
        raise ValueError("Final 'topographic__elevation' field not found.")
    if discharge_field not in grid.at_node:
        raise ValueError(f"Discharge field '{discharge_field}' not found. "
                         "Did you run the FlowAccumulator?")

    # Calculate erosion depth across the entire grid
    erosion_depth_field = grid.at_node[initial_topo] - grid.at_node['topographic__elevation']

    # Iterate through channels and extract data
    # channel_data = {}
    nodes = []
    distances = []
    channel_erosion = []
    channel_discharge = []
    for outlet_node, segments in profiler.data_structure.items():
        # Simply extend the lists with data from each segment
        for segment_data in segments.values():
            segment_nodes = np.array(segment_data['ids'], dtype=int)
            nodes.extend(segment_nodes)
            distances.extend(segment_data['distances'])
            channel_erosion.extend(erosion_depth_field[segment_nodes])
            channel_discharge.extend(grid.at_node[discharge_field][segment_nodes])
        
        # Convert to numpy arrays
        # distance_from_outlet = np.array(all_distances)
        # all_nodes = np.array(all_nodes, dtype=int)

        # # Get erosion and discharge for the channel nodes
        # channel_erosion = erosion_depth_field[all_nodes]
        # channel_discharge = grid.at_node[discharge_field][all_nodes]

        # Store results
    channel_data = {
        'nodes': nodes,
        'distance_from_outlet': distances,
        'erosion_depth': channel_erosion,
        'discharge': channel_discharge,
    }
    return channel_data