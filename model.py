from landlab.components import FastscapeEroder, FlowAccumulator, DepressionFinderAndRouter, LinearDiffuser, LakeMapperBarnes, LinearDiffuser
import numpy as np
from landlab.io import write_esri_ascii
import os
import matplotlib.pyplot as plt
from landlab.plot import imshow_grid
from matplotlib.animation import FuncAnimation

class TopoModel:
    def __init__(self, grid, K_sp, m_sp, n_sp, flow_director, diffusivity = 0, rain_variability=False):
        """
        Initializes the model using a PRE-CONFIGURED grid.
        It doesn't create the grid itself.

        Parameters:
            grid (RasterModelGrid): The grid to use for the model.
            K_sp: erodability.
            m_sp: stream power exponent for drainage area.
            n_sp: stream powerexponent for slope.
            flow_director: The flow director to use for flow routing.
            rain_variability (bool): If True, rainfall rate varies across the grid.
        """
        self.grid = grid
        
        # Instantiate flow accumulator, dpression finder, and fastscape eroder components
        self.fr = FlowAccumulator(self.grid, flow_director = flow_director)
        # self.df = DepressionFinderAndRouter(self.grid)
        self.df = LakeMapperBarnes(self.grid, method=flow_director, redirect_flow_steepest_descent=True, reaccumulate_flow=True)
        self.ld = LinearDiffuser(self.grid, linear_diffusivity = diffusivity)  

        if rain_variability == False: 
            self.fsc = FastscapeEroder(self.grid, K_sp, m_sp, n_sp, discharge_field='drainage_area')
        elif rain_variability == True:
            self.fsc = FastscapeEroder(self.grid, K_sp, m_sp, n_sp, discharge_field='surface_water__discharge')
            # if rainfall_rate is not None:
            # self.grid.at_node['water__unit_flux_in'] = np.ones_like(self.grid.at_node['drainage_area']) * rainfall_rate

    def run_one_step(self, dt):
        """Runs the model for a single timestep.
        
        Parameters:
            dt (float): Duration of the timestep in years.
        """
        # Route flow and handle depressions
        self.fr.run_one_step()
        self.df.run_one_step()
        self.fsc.run_one_step(dt=dt)
        self.ld.run_one_step(dt=dt)

    def run_model(self, runtime, dt, name):
        """Runs the model for a specified duration.
        
        Parameters:
            runtime (float): Total time to run the model in years.
            dt (float): Duration of each timestep in years (default is 1000).
            name (str): Name of the model run, used for saving output files.
        """
        
        num_steps = int(runtime / dt)
        for i in range(num_steps):
            self.run_one_step(dt)
            # # this is how you could save grid to ascii file, every step (or also put below for every 10th step)
            # filename = str(os.getcwd() + '/' + name + '_step{}.asc'.format(i))
            # write_esri_ascii(filename, self.grid, clobber = True)
            if i % 10 == 0: # Print progress
                print(f"Step {i} of {num_steps}")
                
        filename = str(os.getcwd() + '/' + name + '_t{}.asc'.format(runtime))
        write_esri_ascii(filename, self.grid, clobber = True)

    def update_frame(self, frame_num, dt, fig, steps_per_frame):
        """
        Updates the plot. For frame 0, it shows the initial state.
        For all other frames, it runs the model first, then draws.
        """
        # --- THIS IS THE FIX ---
        # Special case for the very first frame
        if frame_num == 0:
            # Don't run the model, just set the time for the title
            current_time = 0
        else:
            # For all other frames, run the simulation first
            for _ in range(steps_per_frame):
                self.run_one_step(dt)
            # Calculate time based on the frame number
            current_time = frame_num * steps_per_frame * dt

        # --- The drawing logic is now the same for all frames ---
        fig.clear()
        ax = fig.add_subplot(1, 1, 1)
        imshow_grid(
            self.grid,
            'topographic__elevation',
            colorbar_label='Elevation (m)',
            cmap='terrain',
            grid_units=('m', 'm'),
        )
        ax.set_title(f'Time: {int(current_time)} years')

    def run_model_with_animation(self, runtime, dt, steps_per_frame=1, filename=None):
        # --- FIX: Add 1 to the frame count for the initial state ---
        total_frames = int(runtime / (dt * steps_per_frame)) + 1
        
        print(f"Generating animation with {total_frames} frames...")
        fig = plt.figure(figsize=(8, 6))

        anim = FuncAnimation(
            fig,
            func=self.update_frame,
            frames=total_frames,
            fargs=(dt, fig, steps_per_frame),
            interval=1,
            blit=False,
            repeat=False
        )

        if filename:
            print(f"Saving animation to {filename}...")
            anim.save(filename, writer='pillow', progress_callback=lambda i, n: print(f'Saving frame {i+1} of {n}'))
            print("Animation saved successfully.")
        else:
            print("Showing animation...")
            plt.show()
        
        plt.close(fig)

        


