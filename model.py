from landlab.components import FastscapeEroder, FlowAccumulator, DepressionFinderAndRouter
import numpy as np
from landlab import RasterModelGrid
from landlab.io import write_esri_ascii
import os


class TopoModel:
    def __init__(self, grid, K_sp, m_sp, n_sp, flow_director, rain_variability=False):
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
        self.df = DepressionFinderAndRouter(self.grid)

        if rain_variability == False: 
            self.fsc = FastscapeEroder(self.grid, K_sp, m_sp, n_sp, discharge_field='drainage_area')
        elif rain_variability == True:
            self.fsc = FastscapeEroder(self.grid, K_sp, m_sp, n_sp, discharge_field='surface_water__discharge')
            # if rainfall_rate is not None:
            # self.grid.at_node['water__unit_flux_in'] = np.ones_like(self.grid.at_node['drainage_area']) * rainfall_rate

    def define_boundaries(self, grid, slope_direction):

        """Defines the boundaries of the grid.
        Sets all boundaries/edges to closed, except the one downstream of slope, 
        which is set as a fixed gradient boundary (BC_NODE_IS_FIXED_GRADIENT)
        
        Parameters:
        grid (RasterModelGrid): The grid to set boundaries for.
        slope_direction (str): Direction of the slope, which determines which boundary is open.
        """
        self.grid = grid

        #Slope direction tells us which boundary should be open. False means open.
        North = "North" != slope_direction
        East = "East" != slope_direction
        South = "South" != slope_direction
        West = "West" != slope_direction

        self.grid.set_closed_boundaries_at_grid_edges(right_is_closed=East, top_is_closed=North, left_is_closed=West, bottom_is_closed=South)
        # self.grid.status_at_node[self.grid.fixed_value_boundary_nodes] = self.grid.BC_NODE_IS_FIXED_GRADIENT
        self.grid.status_at_node[self.grid.fixed_value_boundary_nodes] = self.grid.BC_NODE_IS_FIXED_VALUE
            

    def run_one_step(self, dt):
        """Runs the model for a single timestep.
        
        Parameters:
            dt (float): Duration of the timestep in years.
        """
        # Route flow and handle depressions
        self.fr.run_one_step()
        self.df.map_depressions()
        self.fsc.run_one_step(dt=dt)


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
        
        


