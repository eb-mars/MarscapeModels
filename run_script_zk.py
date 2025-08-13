from model import TopoModel
from make_topography import *
from make_precip import *
from define_boundaries import *
from analyse_streams import *
from landlab.plot import imshow_grid
from landlab.components import ChannelProfiler
from load import load_params_txt
import matplotlib.pyplot as plt
from landlab import load_params

## LOAD PARAMETERS FROM PARAMETER FILE
# params = load_params_txt() ## parameter dictionary
params = load_params("params_zk.txt")  ## load the params file (using landlab load_params)
name = params['model_name']
seed = params['seed']
nrows = params['nrows']
ncols = params['ncols']
cell_size = params['cell_size']
landscape_type = params['landscape_type']
lithology_type = params['lithology_type']
slope = params['slope']
rf = params['rf']
tilt_direction = params['tilt_direction']  # Direction of the tilt
# xy = int(grid_size/cell_size)
flow_director = params['flow_director']
K_sp = params['K_sp']
m_sp = params['m_sp'] 
n_sp = params['n_sp'] 
runoff_rate = params['runoff_rate']
rain_variability = params['rain_variability']
precip_type = params['precip_type']
rainfall_rate = params['rainfall_rate']
dt = params['dt']
runtime = params['runtime']
## boundary conditions (True = closed, False = open);
## in mg.set_closed_boundaries(West, North, East, South);
West = params['West']
North = params['North']
East = params['East']
South = params['South']

# 1. Create the initial topography
# match landscape_type:
#     case 'tilted':
#         grid = create_tilted_landscape(rows=nrows, cols=ncols, cell_size=cell_size, rf=rf, tilt=slope, tilt_direction=tilt_direction, seed=seed)
grid = create_diagonal_tilted_landscape(rows=nrows, cols=ncols, cell_size=cell_size, rf=rf, tilt=slope, tilt_direction="Southeast", seed=seed)

# 2. Add a lithology pattern to it
# match lithology_type:
#     case 'uniform':
grid = add_uniform_lithology(grid, erodibility=K_sp)

# 2a. (optional) Create an array of rainfall to add at each step
##  find the value: -- GEL / time / grid number of nodes
##  create array of same size as grid, with every entry = the value above.
# match precip_type:
#     case 'uniform':
grid = add_uniform_precip(grid, precipitation_rate=rainfall_rate)

# 2b. Store the initial topography to calculate erosion later
grid.at_node['initial_topographic__elevation'] = grid.at_node['topographic__elevation'].copy()

# grid.set_watershed_boundary_condition('topographic__elevation',)
grid = define_boundaries_corner(grid, slope_direction="Southeast")
# grid = define_boundaries_outlet_fixed_gradient(grid, outlet_gradient=0.002)

# 3. Create a model instance with the prepared grid
model = TopoModel(grid, K_sp, m_sp, n_sp, flow_director, rain_variability = rain_variability)

# 4. Run the model
print("Starting model run...")
model.run_model(runtime, dt, name)
# model.run_model_with_animation(runtime, dt)
print("Model run complete.")

# node_status = np.zeros(grid.number_of_nodes)
# node_status[grid.core_nodes] = 1
# imshow_grid(
#     grid,
#     node_status,
#     cmap='viridis',
#     show_elements=True, # Optional: shows grid lines
# )
# plt.show()

# 5. Visualize the result
imshow_grid(
    model.grid, 
    'topographic__elevation',
    cmap='terrain',
    grid_units=('m', 'm')
)
plt.title("Final Topography")
plt.show()

# 6. Optional Checks plotting other data on the grid
imshow_grid(
    model.grid,
    'drainage_area',
    cmap='terrain',
    grid_units=('m', 'm')
)
plt.title("Final Drainage Area")
plt.show()

# imshow_grid(
#     model.grid,
#     'surface_water__discharge',
#     cmap='terrain',
#     grid_units=('m', 'm')
# )
# plt.title("Final surface Water Discharge")
# plt.show()

## CHANNEL PROFILER - handy visualization and channel node ID / channel measurement tool
# plot with channels shown
prf = ChannelProfiler(
    model.grid,
    main_channel_only=False,
    number_of_watersheds=None,
    minimum_channel_threshold=1000,
)
prf.run_one_step()

plt.figure(1)
prf.plot_profiles_in_map_view()
plt.show()

plt.figure(2)
prf.plot_profiles()
plt.show()

# Run FlowAccumulator to calculate discharge on FINAL topography
model.fr.run_one_step() 
# model.df.run_one_step()  # Ensure depressions are handled
channel_data = get_channel_erosion_and_discharge(model.grid, prf)
plot_erosion_vs_discharge(channel_data)
