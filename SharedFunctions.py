import landlab
from landlab.components import FlowAccumulator
from landlab.components import DepressionFinderAndRouter
from landlab.components import ErosionDeposition
from landlab.components import FastscapeEroder

import landlab
from landlab import RasterModelGrid
import numpy as np
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import DepressionFinderAndRouter
from landlab.components import ErosionDeposition
from landlab.components import FastscapeEroder
np.random.seed(seed=5000)
import matplotlib.pyplot as plt
import rasterio
from rasterio.transform import from_origin
grid = RasterModelGrid((500, 500), xy_spacing=10.0)
np.random.seed(seed=5000)
grid.at_node["topographic__elevation"] = (grid.y_of_node / 100 + grid.x_of_node / 100 + np.random.rand(grid.number_of_nodes) * 10 + 10)

grid.set_closed_boundaries_at_grid_edges(
    bottom_is_closed=True,
    left_is_closed=True,
    right_is_closed=True,
    top_is_closed=True )

grid.set_watershed_boundary_condition_outlet_id(0, grid.at_node["topographic__elevation"], -9999.0)
fsc_dt = 100.0
ed_dt = 1.0


grid.at_node["topographic__elevation"].reshape(grid.shape)
plt.imshow(grid.at_node["topographic__elevation"].reshape(grid.shape), origin='lower')
plt.colorbar()

flow_director = "D8"
K_sp = 0.0001
m_sp = 0.7
n_sp = 0.5


def erode(grid, flow_director, K_sp, m_sp, n_sp):
    fr = FlowAccumulator(grid, flow_director="D8")
    df = DepressionFinderAndRouter(grid)
    fsc = FastscapeEroder(grid, K_sp, m_sp, n_sp)

