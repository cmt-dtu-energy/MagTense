import os
import subprocess
import numpy as np
import math
import random as rand

import util_plot

# Path of MagTense
DIRPATH_MAGTENSE = os.path.join(os.path.dirname(__file__), '../../executable')
# Input and output files required for MagTense defined in ../util/io.txt
# Settings for MagTense: settings_std_alone.txt
FILENAME_SETTINGS = 'settings_std_alone.txt'
# Required input file for magnetic tiles: test_tiles.txt
FILENAME_INPUT_TILES = 'test_tiles.txt'
# Required input file for evaluation points: std_alone_points.txt
FILENAME_INPUT_POINTS = 'std_alone_points.txt'
# Output file for iterated magnetic tiles: tiles_out.txt
FILENAME_OUTPUT_TILES = 'tiles_out.txt'
# Output file for calculated magentic field in the specified points: H_out.txt
FILENAME_OUTPUT_FIELD = 'H_out.txt'

def create_io():
        with open(os.path.join(DIRPATH_MAGTENSE, 'io.txt'), "w") as file:
                file.write('\'' +FILENAME_INPUT_TILES + '\'\n\'' + FILENAME_OUTPUT_TILES + '\'\n\'' + FILENAME_SETTINGS + '\'\n\'' \
                        + FILENAME_INPUT_POINTS + '\'\n\'' + FILENAME_OUTPUT_FIELD + '\'')

def set_settings(max_error=0.00001, max_it=500, iterate_solution=1, return_field=1):
        create_io()
        with open(os.path.join(DIRPATH_MAGTENSE, FILENAME_SETTINGS), "r") as file:
                rest_settings = file.readlines()[4:]
        with open(os.path.join(DIRPATH_MAGTENSE, FILENAME_SETTINGS), "w") as file:
                file.write(("{max_error:f} :max relative error \n"
                        "{max_it} \t:max no of iterations \n"
                        "{it_sol} \t:iterate solution (0 false, otherwise true) \n"
                        "{return_field} \t:return the field in the points specified in the appropriate file (0 false, otherwise true) \n"
                        ).format(max_error=max_error, max_it=max_it, it_sol=iterate_solution, return_field=return_field))
        with open(os.path.join(DIRPATH_MAGTENSE, FILENAME_SETTINGS), "a") as file:
                file.write(''.join(rest_settings))


class Tile():
        def __init__(self):
                self.grid_pos = np.zeros(3)
                self.center_pos = np.zeros(shape=(3), dtype=np.float64, order='F') # r0, theta0, z0
                self.dev_center = np.zeros(shape=(3), dtype=np.float64, order='F') # dr, dtheta, dz
                self.size = np.zeros(shape=(3), dtype=np.float64, order='F') # a, b, c
                self.M = np.zeros(shape=(3), dtype=np.float64, order='F') # Mx, My, Mz
                self.u_ea = np.array([1, 0, 0], dtype=np.float64, order='F') # Easy axis
                self.u_oa1 = np.array([0, 1, 0], dtype=np.float64, order='F')
                self.u_oa2 = np.array([0, 0, 1], dtype=np.float64, order='F')
                self.M_rem = 0.
                self.mu_r_ea = 1.
                self.mu_r_oa = 1.
                self.tile_type = 2 # 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv, 10 = ellipsoid
                self.magentic_type = 1 # 1 = hard magnet, 2 = soft magnet
                self.stfcn_index = 1 # default index into the state function
                self.incl_it = 1 # if equal to zero the tile is not included in the iteration
                self.offset = np.zeros(shape=(3), dtype=np.float64, order='F') # offset of global coordinates
                self.rot = np.zeros(shape=(3), dtype=np.float64, order='F')
                self.color = np.array([1, 0, 0], dtype=np.float64, order='F')
                self.use_sym = 0 # whether to exploit symmetry
                self.sym_op = np.array([1, 1, 1], dtype=np.float64, order='F') # 1 for symmetry and -1 for anti-symmetry respectively to the planes
                self.M_rel = 0.

        def __str__(self):
                return ("Tile at grid position ({grid_x},{grid_y},{grid_z}) with coordinates x={x}, y={y}, z={z}.").format(\
                        grid_x = self.grid_pos[0], grid_y = self.grid_pos[1], grid_z=self.grid_pos[2], \
                        x=self.offset[0], y = self.offset[1], z=self.offset[2])
        
        def set_grid_pos(self, grid_pos):
                self.grid_pos[:] = grid_pos[:]
        
        def set_size(self, size):
                self.size[:] = size[:]
        
        def get_size(self):
                return self.size

        def set_tile_type(self, tile_type):
                self.tile_type = tile_type

        def get_tile_type(self):
                return self.tile_type
        
        def set_center_pos(self, center_pos):
                self.center_pos[:] = center_pos[:]

        def get_center_pos(self):
                return self.center_pos

        def set_dev_center(self, dev_center):
                self.dev_center[:] = dev_center[:]

        def get_dev_center(self):
                return self.dev_center

        def set_offset(self, offset):
                self.offset[:] = offset[:]

        def get_offset(self):
                return self.offset
        
        def set_rotation(self, rotation):
                self.rot[:] = rotation[:]
        
        def get_rotation(self):
                return self.rot

        def set_easy_axis(self, easy_axis):
                self.u_ea[:] = easy_axis[:]
                self.M = self.M_rem * self.u_ea
        
        def set_oa1(self, other_axis):
                self.u_oa1[:] = other_axis[:]
        
        def set_oa2(self, other_axis):
                self.u_oa2[:] = other_axis[:]
        
        def set_remanence(self, M_rem):
                self.M_rem = M_rem
        
        def set_mag_angle(self, theta):
                azimuth = theta[0]
                polar_angle = theta[1]
                self.set_easy_axis([math.sin(polar_angle) * math.cos(azimuth), math.sin(polar_angle) * math.sin(azimuth), math.cos(polar_angle)])
                self.set_oa1([math.sin(polar_angle) * math.sin(azimuth), math.sin(polar_angle) * (-math.cos(azimuth)), 0])
                self.set_oa2([0.5*math.sin(2*polar_angle) * math.cos(azimuth), 0.5*math.sin(2*polar_angle) * math.sin(azimuth), -math.pow(math.sin(polar_angle),2)])
        
        def get_M(self):
                return self.M
        
        def get_color(self):
                return self.color
        
        # Input/outfiles files for magnetic tiles: test_tiles.txt, tiles_out.txt
        # Format(Number of tiles in first line followed by tiles with 13 lines of properties each):
        # 4 :no. of tiles
        # 0,0,0,0,0,0 :r0,theta0,z0,dr,dtheta,dz
        # 0.1,0.1,0.01 :a,b,c
        # 954929.6586,0,0 :Mx,My,Mz
        # 1,0,0 :u_ea_x,u_ea_y,u_ea_z
        # 0,-1,0 :u_oa1_x,u_oa1_y,u_oa1_z
        # 0,0,1 :u_oa2_x,u_oa2_y,u_oa2_z
        # 954929.6586,1,1 :Mrem,mu_r_ea,mu_r_oa
        # 2,1,1,1 :tileType,magnetType,stFcnIndex,inclIter
        # 0.2,0.4,0 :offset(1,2,3)
        # 0,0,0 :rot(1,2,3)
        # 0,1,1,1 :useSymm,symmOps(1,2,3)
        # 1,0,0 :color(1,2,3)
        # 0 :Mrel
        def tile2txt(self):
                txt = ("{r0},{theta0},{z0},{dr},{dtheta},{dz} \t\t:r0,theta0,z0,dr,dtheta,dz \n"
                        "{a},{b},{c} \t\t\t\t:a,b,c \n"
                        "{M_x},{M_y},{M_z} \t\t:Mx,My,Mz \n"
                        "{u_ea_x},{u_ea_y},{u_ea_z} \t\t\t\t:u_ea_x,u_ea_y,u_ea_z \n"
                        "{u_oa1_x},{u_oa1_y},{u_oa1_z} \t\t\t\t:u_oa1_x,u_oa1_y,u_oa1_z \n"
                        "{u_oa2_x},{u_oa2_y},{u_oa2_z} \t\t\t\t\t:u_oa2_x,u_oa2_y,u_oa2_z \n"
                        "{M_rem},{mu_r_ea},{mu_r_oa} \t\t\t:Mrem,mu_r_ea,mu_r_oa \n"
                        "{tile_type},{magentic_type},{stfcn_index},{incl_it} \t\t\t\t:tileType,magnetType,stFcnIndex,inclIter \n"
                        "{offset0},{offset1},{offset2} \t\t\t:offset(0,1,2) \n"
                        "{rot0},{rot1},{rot2} \t\t\t\t:rot(0,1,2) \n"
                        "{use_sym},{sym_op0},{sym_op1},{sym_op2} \t\t\t\t:useSymm,symmOps(0,1,2) \n"
                        "{color0},{color1},{color2} \t\t\t\t\t:color(0,1,2) \n"
                        "{M_rel} \t\t\t\t\t:Mrel \n").format(r0=self.center_pos[0], theta0=self.center_pos[1], z0=self.center_pos[0], \
                        dr=self.dev_center[0], dtheta=self.dev_center[1], dz=self.dev_center[2], \
                        a=self.size[0], b=self.size[1], c=self.size[2], M_x=self.M[0], M_y=self.M[1], M_z=self.M[2], \
                        u_ea_x=self.u_ea[0], u_ea_y=self.u_ea[1], u_ea_z=self.u_ea[2], \
                        u_oa1_x=self.u_oa1[0], u_oa1_y=self.u_oa1[1], u_oa1_z=self.u_oa1[2], \
                        u_oa2_x=self.u_oa2[0], u_oa2_y=self.u_oa2[1], u_oa2_z=self.u_oa2[2], \
                        M_rem=self.M_rem, mu_r_ea=self.mu_r_ea, mu_r_oa=self.mu_r_oa, \
                        tile_type=self.tile_type, magentic_type=self.magentic_type, stfcn_index=self.stfcn_index, incl_it=self.incl_it, \
                        offset0=self.offset[0], offset1=self.offset[1], offset2=self.offset[2], rot0=self.rot[0], rot1=self.rot[1], rot2=self.rot[2], \
                        use_sym=self.use_sym, sym_op0=self.sym_op[0], sym_op1=self.sym_op[1], sym_op2=self.sym_op[2], \
                        color0=self.color[0], color1=self.color[1], color2=self.color[2], M_rel=self.M_rel)
                return txt

        def txt2tile(self, l_update_txt):
                var = []
                for line in l_update_txt:
                        var.append(line.strip().split(","))
                self.center_pos = np.array([float(var[0][0]), float(var[0][1]), float(var[0][2])]) # r0, theta0, z0
                self.dev_center = np.array([float(var[0][3]), float(var[0][4]), float(var[0][5])]) # dr, dtheta, dz
                self.size = np.array([float(var[1][0]), float(var[1][1]), float(var[1][2])]) # a, b, c
                self.M = np.array([float(var[2][0]), float(var[2][1]), float(var[2][2])]) # Mx, My, Mz
                self.u_ea = np.array([float(var[3][0]), float(var[3][1]), float(var[3][2])]) # Easy axis
                self.u_oa1 = np.array([float(var[4][0]), float(var[4][1]), float(var[4][2])])
                self.u_oa2 = np.array([float(var[5][0]), float(var[5][1]), float(var[5][2])])
                self.M_rem = float(var[6][0])
                self.mu_r_ea = float(var[6][1])
                self.mu_r_oa = float(var[6][2])
                self.tile_type = int(var[7][0]) # 1 = cylinder, 2 = prism, 3 = ellipsoid
                self.magentic_type = int(var[7][1]) # 1 = hard magnet, 2 = soft magnet
                self.stfcn_index = int(var[7][2]) # default index into the state function
                self.incl_it = int(var[7][3]) # if equal to zero the tile is not included in the iteration
                self.offset = np.array([float(var[8][0]), float(var[8][1]), float(var[8][2])])
                self.rot = np.array([float(var[9][0]), float(var[9][1]), float(var[9][2])])
                self.use_sym = int(var[10][0]) # whether to exploit symmetry
                self.sym_op = np.array([float(var[10][1]), float(var[10][2]), float(var[10][3])]) # 1 for symmetry and -1 for anti-symmetry respectively to the planes
                self.color = np.array([float(var[11][0]), float(var[11][1]), float(var[11][2])])
                self.M_rel = float(var[12][0])                             

def create_tiles_txt(tiles):
        input = ("{n_tiles} \t\t\t\t\t:number of tiles \n").format(n_tiles=len(tiles))
        for tile in tiles:
                input = input + tile.tile2txt()
        with open(os.path.join(DIRPATH_MAGTENSE, FILENAME_INPUT_TILES), "w") as file:
                file.write(input)

def update_tiles(tiles):
        with open(os.path.join(DIRPATH_MAGTENSE, FILENAME_OUTPUT_TILES), "r") as file:
                # Read out the number of tiles in file
                n_tiles = int(file.readline().strip())
                if n_tiles != len(tiles):
                        print("Numbers of tiles are not matching to the update file!")
                        exit()
                lines = file.readlines()
                # Go through file with 13 lines for each tile
                for i in range(0,n_tiles):
                        l_update_txt = lines[(13*(i)):]
                        tiles[i].txt2tile(l_update_txt)
        return tiles                               

# class EvalPoint():
#         def __init__(self, pos=None):
#                 if pos is None:
#                         self.x = 0
#                         self.y = 0
#                         self.z = 0
#                 else:
#                         self.x = pos[0]
#                         self.y = pos[1]
#                         self.z = pos[2]

#         def __str__(self):
#                 s = ("Evaluation point at x: {x}, y: {y}, z: {z}").format(x=self.x, y=self.y, z=self.z) 
#                 return s
        
#         def set_pos(self, pos):
#                 self.x = pos[0]
#                 self.y = pos[1]
#                 self.z = pos[2]
        
#         def set_coordinates(self, x, y, z):
#                 self.x = x
#                 self.y = y
#                 self.z = z


class MagneticFieldIntensity():
        def __init__(self, points=None, H=None):
                self.field = {}
                if points is not None:
                        for i, point in enumerate(points):
                                point_hash = (point[0], point[1], point[2])
                                self.field[point_hash] = H[i]

        def get_average_magnetic_flux(self):
                norm = self.get_norm_magnetic_flux()
                average = sum(norm.values())/len(norm)
                return average

        def get_p2p(self):
                norm = self.get_norm_magnetic_flux()
                return max(norm.values())-min(norm.values())
        
        def get_norm_magnetic_flux(self):
                norm = {}
                # vacuum permeability
                mu0 = 4*math.pi*1e-7
                for point in self.field:
                        norm[point] = np.linalg.norm(self.field[point])*mu0
                return norm
        
        def get_intensity_from_eval(self, points):
                with open(os.path.join(DIRPATH_MAGTENSE, FILENAME_OUTPUT_FIELD), "r") as file:
                        for i, line in enumerate(file.readlines()):
                                values = line.strip().split(",")
                                values = [value.strip() for value in values]
                                point_hash = (points[i][0], points[i][1], points[i][2])
                                self.field[point_hash] = np.array([float(values[0]), float(values[1]), float(values[2])])
                                

def create_points_txt(points):
        input = ("{n_points} :number of points \n").format(n_points=len(points))
        for point in points:
                # input = input + ("{p_x}, {p_y}, {p_z} :x,y,z\n").format(p_x=point.x, p_y=point.y, p_z=point.z)
                input = input + ("{p_x}, {p_y}, {p_z} :x,y,z\n").format(p_x=point[0], p_y=point[1], p_z=point[2])  
        with open(os.path.join(DIRPATH_MAGTENSE, FILENAME_INPUT_POINTS), "w") as file:
                file.write(input)


class Grid():
        def __init__(self, places, area):
                self.places = np.asarray(places)
                self.area = np.asarray(area)
                self.size_tile = self.area/self.places
                self.grid = np.zeros(self.places.tolist())

        def get_tiles(self, grid_positions = [], n_tiles = None):
                tiles = []
                # Fill grid randomly
                if not grid_positions:
                        if n_tiles is None:
                                n_tiles = rand.randrange(np.prod(self.places))
                        else:
                                n_tiles = n_tiles
                        grid_positions = []
                        while n_tiles > 0:
                                new_pos = [rand.randrange(self.places[0]-1), rand.randrange(self.places[1]-1), rand.randrange(self.places[2]-1)]
                                if new_pos in grid_positions:
                                        continue
                                else:
                                        grid_positions.append(new_pos)
                                        n_tiles = n_tiles - 1
                for pos in grid_positions:
                        if True in np.greater_equal(np.asarray(pos), self.places):
                                print(("Desired position {} is not in the grid!").format(pos))
                                exit()
                        self.grid[pos] = 1
                        new_tile = Tile()
                        # Set grid position of tile
                        new_tile.set_grid_pos(pos)
                        # Set size of tile
                        new_tile.set_size(self.size_tile)
                        # Extract Cartesian coordinates of tile
                        new_tile.set_offset(np.around((new_tile.grid_pos * self.size_tile) + self.size_tile/2, decimals=9))
                        tiles.append(new_tile)
                return tiles
        
        def set_eval_points(self, n_points, mode):
                # points = []
                counter = 0
                points = np.zeros(shape=(n_points[0]*n_points[1]*n_points[2],3), dtype=np.float64, order='F')
                if mode == "uniform":
                        seg = self.area/np.asarray(n_points)
                        for i in range(0,n_points[0]):
                                for j in range(0, n_points[1]):
                                        for k in range(0, n_points[2]):
                                                # new_point = EvalPoint()
                                                # new_point.set_coordinates(i*seg[0]+seg[0]/2, j*seg[1]+seg[1]/2, k*seg[2]+seg[2]/2)
                                                # points.append(new_point)
                                                points[counter] = [i*seg[0]+seg[0]/2, j*seg[1]+seg[1]/2, k*seg[2]+seg[2]/2]
                                                counter = counter + 1
                elif mode == "center":
                        center = self.area/2
                        seg = (self.area[0]/10)/n_points[0]
                        seg_angle = math.pi/n_points[1]
                        seg_layer = self.area[2]/n_points[2]
                        for i in range(0,n_points[0]):
                                for j in range(0, n_points[1]):
                                        for k in range(0, n_points[2]):
                                                # new_point = EvalPoint()
                                                # new_pos = [i*seg, j*seg_angle, k*seg_layer] + center
                                                # new_point.set_pos(new_pos)
                                                # points.append(new_point)
                                                points[counter] = [i*seg, j*seg_angle, k*seg_layer] + center
                                                counter = counter + 1
                elif mode == "costumized":
                        pass
                else:
                        print("Please specify a valid area of interest!")
                        exit()
                return points

        def clear_grid(self):
                self.grid = np.zeros(self.places.tolist())

def setup(places, area, n_tiles=0, filled_positions=None, mag_angles=[], eval_points=[20, 20, 1], eval_mode="uniform",\
        B_rem=1.2, max_error=0.00001, max_it=500, iterate_solution=1, return_field=1):
        set_settings(max_error, max_it, iterate_solution, return_field)
        # Check format of input parameters
        if len(places) != 3:
                print("Format of number of possible magnets in each axis is not correct!")
                exit
        elif len(area) != 3:
                print("Format of area is not correct!")
                exit
        elif len(eval_points) != 3:
                print("Format of eval points is not correct!")
                exit
        meshgrid = Grid(places, area)
        # Extract coordinates of evaluation points
        points = meshgrid.set_eval_points(eval_points, eval_mode)
        # Fill grid with magnetic tiles
        if filled_positions is None:
                tiles = meshgrid.get_tiles(n_tiles)
        else:
                tiles = meshgrid.get_tiles(grid_positions=filled_positions)
        # Assign magnetization angles for tiles
        if not mag_angles:
                for i in range(0,len(tiles)):
                        mag_angles = mag_angles.append(2 * math.pi * rand.random())       
        for i, tile in enumerate(tiles):
                tile.set_remanence(B_rem / (4*math.pi*1e-7))
                tile.set_mag_angle(mag_angles[i])        
        return tiles, points, meshgrid

# Function for running standalone MagTense
def run_simulation(tiles, points, grid=None, plot=False, output=False):
        create_tiles_txt(tiles)
        create_points_txt(points)
        filepath = os.path.join(DIRPATH_MAGTENSE, 'MagTense_StandAlone.exe')
        result = subprocess.run(filepath, cwd=DIRPATH_MAGTENSE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if output==True:
                print(result.stdout)
        updated_tiles = update_tiles(tiles)
        H = MagneticFieldIntensity()
        H.get_intensity_from_eval(points)
        if plot==True:
                if grid is None:
                        util_plot.create_plot(updated_tiles, H)
                else:
                        util_plot.create_plot(updated_tiles, H, grid=grid)
        return updated_tiles, H
