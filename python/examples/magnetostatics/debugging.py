import sys
import importlib
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from pathlib import Path
from magtense.magstatics import Tiles, run_simulation
from magtense.utils import create_plot, validation, run_simulation
#To use my local version of utils.py and not the installed one. Is there a better solution?
# Path to local magtense/src folder
#local_path = Path(r"C:\Users\eroso\Desktop\Magtense\python\src").resolve()

# Temporarily insert path and import
#sys.path.insert(0, str(local_path))
#try:
    # Force reload in case a version is already imported
    #import magtense.utils as utils
    #importlib.reload(utils)

#from magtense.utils import Tiles, create_plot, validation, run_simulation
#finally:
    # Remove the path after importing to avoid messing with other imports
    #sys.path.pop(0)

mu0 = 4 * np.pi * 1e-7
tiles = Tiles(
    n=1,
    M_rem=1 / mu0,
    tile_type=[2],
    color=[[1, 0, 0]]
)

# 0: Prism 1
tiles.size = ([0.1, 0.1, 0.1], 0)
tiles.offset = ([0.0, 0.0, 0.0], 0)

# 1: Prism 2
#tiles.size = ([0.1, 0.3, 0.2], 1)
#tiles.offset = ([1.0, 1.0, 1.0], 1)

pts = np.array([[3.0, 2.0, 2.0], [2.0,2.0,2.0]], dtype=np.float64, order='F')
obs_size = np.array([1.0,1.0,1.0])
# n_points = [10, 100, 80]
# seg = [0.9, 1, 0.8] / np.asarray(n_points)
# pts = [[(i + 0.5) * seg[0], (j + 0.5) * seg[1], (k + 0.5) * seg[2]]
#         for i in range(n_points[0])
#         for j in range(n_points[1])
#         for k in range(n_points[2])]
# pts = np.asarray(pts, dtype=np.float64, order='F')

updated_tiles, H_demag = run_simulation(tiles, pts, obs_size=obs_size)

print(H_demag)