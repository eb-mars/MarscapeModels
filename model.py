# model.py
from landlab.components import FastscapeEroder, FlowAccumulator

class TopoModel:
    def __init__(self, grid, K_sp_field='K_sp'):
        """
        Initializes the model using a PRE-CONFIGURED grid.
        It doesn't create the grid itself.

        Parameters:
        - grid (RasterModelGrid): The grid to use for the model.
        - K_sp_field (str): The name of the field that contains erodibility values
                            (default is 'K_sp').
        """
        self.grid = grid
        
        # Instantiate components, telling the eroder to use the specified field
        self.flow_accumulator = FlowAccumulator(self.grid)
        self.eroder = FastscapeEroder(self.grid, K_sp=K_sp_field)

    def run_one_step(self, dt):
        """Runs the model for a single timestep.
        
        Parameters:
        - dt (float): The duration of the timestep in years.
        """
        self.flow_accumulator.run_one_step()
        self.eroder.run_one_step(dt)

    def run_model(self, runtime, dt=1000):
        """Runs the model for a specified duration.
        
        Parameters:
        - runtime (float): Total time to run the model in years.
        - dt (float): Duration of each timestep in years (default is 1000)."""
        num_steps = int(runtime / dt)
        for i in range(num_steps):
            self.run_one_step(dt)
            if i % 10 == 0: # Print progress
                print(f"Step {i} of {num_steps}")