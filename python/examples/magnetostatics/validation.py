#%%
import math

from magtense.magstatics import Tiles
from magtense.utils.eval import validation


def valid_prism():
    tile_val = Tiles(
        n=1,
        size=[0.6, 0.1, 0.3],
        offset=[0.5, 0.4, 0.1],
        rot=[math.pi/2, -math.pi/3, math.pi/4],
        tile_type=2,
        M_rem=1.2/(4*math.pi*1e-7),
        easy_axis=[0.35355339, 0.61237244, 0.70710678],
        color=[1, 0, 0]
    )
    validation('prism', tile_val, [0.5, 0.4, 0.1])


def valid_circpiece():
    tile_val = Tiles(
        n=1,
        center_pos=[0.45, math.pi/3, 0.35],
        dev_center=[0.5, math.pi/7, 0.15],
        offset=[0.1, 0.3, 0.2],
        rot=[0, 0, 0],
        tile_type=3,
        M_rem=1.2/(4*math.pi*1e-7),
        easy_axis=[-0.3095974 , -0.22493568,  0.92387953],
        color=[1, 0, 0]
    )
    validation('circpiece', tile_val,
               [0.406328622087633, 0.8631363102808783, 0.55])


def valid_circ_inv():
    tile_val = Tiles(
        n=1,
        center_pos=[0.3, math.pi/0.55, 0.6],
        dev_center=[0.15, math.pi/6, 0.4],
        offset=[0.3, 0.5, 0.1],
        rot=[0, 0, 0],
        tile_type=4,
        M_rem=1.2/(4*math.pi*1e-7),
        easy_axis=[0.41562694, 0.41562694, 0.80901699],
        color=[1, 0, 0]
    )
    validation('circpiece_inverted', tile_val,
               [0.6271937452259475, 0.27251823835641853, 0.7])


def valid_sphere():
    tile_val = Tiles(
        n=1,
        size=[0.25, 0, 0],
        offset=[0.3, 0.5, 0.1],
        tile_type=6,
        M_rem=1.2/(4*math.pi*1e-7),
        easy_axis=[1, 0, -1],
        color=[1, 0, 0]
    )
    validation('sphere', tile_val, [0.3, 0.4, 0.15],
               [0.25, 0.25, 0.75], unit='T')


def valid_spheroid():
    tile_val = Tiles(
        n=1,
        size=[0.5, 2, 0.5],
        offset=[0.1, 0.3, 0.45],
        rot=[-2, 1, 2],
        tile_type=7,
        M_rem=1.2/(4*math.pi*1e-7),
        easy_axis=[0, 1, 0],
    )
    validation('spheroid', tile_val, [0.3, 0.45, 0.2],
               [1, 0.75, 1], unit='T')


def valid_tetrahedron():
    tile_val = Tiles(
        n=1,
        vertices=[[2.5,3,1],[2,1,4],[1.5,4,3],[4.5,5,2]],
        tile_type=5,
        M_rem=1.2 / (4*math.pi*1e-7),
        easy_axis=[0.324264068, 0.734846928, 0.891545179],
        color = [1, 0, 0]
    )    
    validation('tetrahedron', tile_val, [3, 3, 2.5])
# %%
