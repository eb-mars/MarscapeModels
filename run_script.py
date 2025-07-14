# run_simulation.py
from model import TopoModel
from make_topography import create_tilted_landscape, add_checkerboard_lithology
from landlab.plot import imshow_grid
import matplotlib.pyplot as plt

# --- Experiment 1: Tilted landscape with two rock types --- 🏗️🧱

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