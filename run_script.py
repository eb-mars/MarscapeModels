from model import TopoModel
from make_topography import *
from landlab.plot import imshow_grid
from load import load_params_txt
import matplotlib.pyplot as plt

# --- Experiment 1: Tilted landscape with two rock types ---

#### LOAD PARAMETERS FROM PARAMETER FILE
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
initial_grid = create_tilted_landscape(rows=nrows, cols=ncols, cell_size=cell_size, rf=rf, tilt=slope, tilt_direction=tilt_direction)

# 2. Add a lithology pattern to it
final_grid = add_uniform_lithology(initial_grid)

# 2a. (optional) Create an array of rainfall to add at each step
##  find the value: -- GEL / time / grid number of nodes
##  create array of same size as grid, with every entry = the value above. 

# 3. Create a model instance with the prepared grid
model_run = TopoModel(final_grid, K_sp, m_sp, n_sp, flow_director)

# 3a. Define the boundaries of the grid
model_run.define_boundaries(final_grid, tilt_direction)

# 4. Run the model
print("Starting model run...")
model_run.run_model(runtime=runtime, dt=dt) # Run for 500,000 years
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
