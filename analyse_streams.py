import numpy as np
from scipy.stats import spearmanr, linregress
import matplotlib.pyplot as plt

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

def plot_erosion_vs_discharge(channel_data):
    """
    Plots erosion depth against discharge.

    Parameters:
        channel_data (dict): Dictionary containing 'distance_from_outlet',
                             'erosion_depth', and 'discharge'.
    """
    erosion = channel_data['erosion_depth']
    discharge = channel_data['discharge']

    # Create plots
    plt.figure(figsize=(5, 4))
    plt.scatter(discharge, erosion, alpha=0.6)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Discharge ($m^3/yr$)')
    plt.ylabel('Erosional Depth (m)')
    plt.title('Erosion vs. Discharge')

    res_s = spearmanr(discharge, erosion)
    res_l = linregress(discharge, erosion)
    graph_fit = linregress(np.log10(discharge), np.log10(erosion))

    # Plot the regression line
    x = np.linspace(np.nanmin(discharge), np.nanmax(discharge), 100)
    y = 10 ** (graph_fit.intercept + graph_fit.slope * np.log10(x))
    plt.plot(x, y, c='#698fd6', lw=2)  #ad664b

    textstr = '\n'.join((
        r"Spearman's rank correlation = %.3f" % (res_s.correlation, ),
        f'Pearson correlation = {res_l.rvalue:.3f}'))

    props = dict(facecolor='white', alpha=0.5)
    plt.text(0.17, 0.12, textstr, transform=plt.gca().transAxes,
            verticalalignment='top', bbox=props,)

    plt.grid(True)
    plt.tight_layout()
    plt.show()