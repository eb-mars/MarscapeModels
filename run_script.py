from model import TopoModel
from make_topography import *
from make_precip import *
from landlab.plot import imshow_grid
from landlab.components import ChannelProfiler
from load import load_params_txt
import matplotlib.pyplot as plt

## LOAD PARAMETERS FROM PARAMETER FILE
params = load_params_txt() ## parameter dictionary
name = params['model_name']
seed = params['seed']
nrows = params['nrows']
ncols = params['ncols']
cell_size = params['cell_size']
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
rainfall_rate = params['rainfall_rate']
dt = params['dt']
runtime = params['runtime']
## boundary conditions (True = closed, False = open);
## in mg.set_closed_boundaries(West, North, East, South);
West = params['West']
North = params['North']
East = params['East']
South = params['South']

# --- Experiment 1: Tilted landscape with two rock types ---
# 1. Create the initial topography
grid = create_tilted_landscape(rows=nrows, cols=ncols, cell_size=cell_size, rf=rf, tilt=slope, tilt_direction=tilt_direction, seed=seed)

# 2. Add a lithology pattern to it
grid = add_uniform_lithology(grid, erodibility=K_sp)

# 2a. (optional) Create an array of rainfall to add at each step
##  find the value: -- GEL / time / grid number of nodes
##  create array of same size as grid, with every entry = the value above. 
grid = add_uniform_precip(grid, precipitation_rate=rainfall_rate)

# 3. Create a model instance with the prepared grid
model_run = TopoModel(grid, K_sp, m_sp, n_sp, flow_director, rain_variability = rain_variability)

# 3a. Define the boundaries of the grid
model_run.define_boundaries(grid, tilt_direction)

# 4. Run the model
print("Starting model run...")
model_run.run_model(runtime, dt, name)
print("Model run complete.")

# 5. Visualize the result
imshow_grid(
    model_run.grid, 
    'topographic__elevation',
    cmap='terrain',
    grid_units=('m', 'm')
)
plt.title("Final Topography")
plt.show()


## CHANNEL PROFILER - handy visualization and channel node ID / channel measurement tool
# plot with channels shown
prf = ChannelProfiler(
    model_run.grid,
    main_channel_only=True,
    number_of_watersheds=None,
    minimum_channel_threshold=model_run.grid.dy*model_run.grid.dx,
)
prf.run_one_step()

plt.figure(1)
prf.plot_profiles_in_map_view()
plt.show()

## to access data about the channels, e.g. where are node IDs;
outlet_nodes = prf.data_structure.keys() # node IDs of the outlet nodes
prf.data_structure[list(prf.data_structure.keys())[0]].keys(); # node IDs representing the upstream & downstream segments of the river channel, for the first outlet node
prf.data_structure[list(prf.data_structure.keys())[0]][list(prf.data_structure[list(prf.data_structure.keys())[0]].keys())[0]]["ids"] ## node IDs representing the whole channel, for the channel of the FIRST outlet node
prf.data_structure[list(prf.data_structure.keys())[0]][list(prf.data_structure[list(prf.data_structure.keys())[0]].keys())[0]]["distances"] ## distances between segements for the whole channel, for the channel of the  outlet node

prf.data_structure[list(prf.data_structure.keys())[1]].keys() # node IDs representing the upstream & downstream segments of the river channel, for the SECOND outlet node



# 6. Optional Checks plotting other data on the grid
# imshow_grid(
#     model_run.grid, at='node',
#     'drainage_area',
#     cmap='terrain',
#     grid_units=('m', 'm')
# )
# plt.title("Final Drainage Area")
# plt.show()

# imshow_grid_at_node(
#     model_run.grid,  at='node',
#     'surface_water__discharge',
#     cmap='terrain',
#     grid_units=('m', 'm')
# )
# plt.title("Final surface Water Discharge")
# plt.show()

# imshow_grid_at_node(
#     model_run.grid,  at='node',
#     'K_sp',
#     cmap='terrain',
#     grid_units=('m', 'm')
# )
# plt.title("Final K_Sp")
# plt.show()

# imshow_grid_at_node(
#     model_run.grid,  at='node',
#     'water__unit_flux_in',
#     cmap='terrain',
#     grid_units=('m', 'm')
# )
# plt.title("Final water unit flux in")
# plt.show()





