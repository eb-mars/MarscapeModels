# run_simulation.py
from model import TopoModel
from make_topography import create_tilted_landscape, add_checkerboard_lithology
from landlab.plot import imshow_grid
from load import load_params
import matplotlib.pyplot as plt

# --- Experiment 1: Tilted landscape with two rock types --- 🏗️🧱

#### LOAD PARAMETERS FROM PARAMETER FILE
params = load_params(txt); ## parameter dictionary
name = params['model_name'];
seed = params['seed']; 
grid_size = params['grid_size']; 
cell_size = params['cell_size'];
slope = params['slope']; 
rf = params['rf'];
xy = int(grid_size/cell_size);
K_sp = params['K_sp']; 
m_sp = params['m_sp']; 
n_sp = params['n_sp']; 
runoff_rate = params['runoff_rate']; 
dt = params['dt']; 
steps = params['steps']; 
## boundary conditions (True = closed, False = open);
## in mg.set_closed_boundaries(West, North, East, South);
West = params['West'];
North = params['North'];
East = params['East'];
South = params['South']; 


# 1. Create the initial topography
initial_grid = create_tilted_landscape(rows=50, cols=100)

# 2. Add a lithology pattern to it
final_grid = add_checkerboard_lithology(initial_grid)

# 3. Create a model instance with the prepared grid
model_run = TopoModel(final_grid, K_sp_field='K_sp')

# 4. Run the model
print("Starting model run...")
model_run.run_model(runtime=500000) # Run for 500,000 years
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