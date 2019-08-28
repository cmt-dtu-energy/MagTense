from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm


def plot_cube(axes, size, offset, rotation, M, color):
    ax = axes

    # Define the vertices of the unit cubic and move them in order to center the cube on origo
    ver = np.array([[1, 1, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1], [1, 0, 1], [1, 0, 0], [0, 0, 0]]) - 0.5
    ver_cube = ver * size
    R = get_rotmat(rotation)
    ver_cube = (np.dot(R, ver_cube.T)).T
    ver_cube = ver_cube + offset

    print(ver_cube)

    # Define the faces of the unit cubic
    fac = np.array([[0, 1, 2, 3], [3, 2, 4, 5], [5, 6, 7, 4],
                    [0, 1, 7, 6], [5, 6, 0, 3], [1, 2, 4, 7]])
    surfaces = ver_cube[fac]

    # plot sides
    ax.add_collection3d(Poly3DCollection(surfaces, facecolors=color, linewidths=1, edgecolors=color, alpha=.25))

    # Plot vector of magnetization in the center of the cube
    ax.quiver(offset[0], offset[1], offset[2], M[0], M[1], M[2], color=color, length=np.linalg.norm(size)/4, pivot='middle', normalize=True)


def plot_cylindrical(axes, center_pos, dev_center, offset, rotation, M, color):
    ax = axes
    resolution = 100
    r, theta, z = center_pos[:]
    dr, dtheta, dz = dev_center[:]
    xc = center_pos[0]*math.cos(center_pos[1])
    yc = center_pos[0]*math.sin(center_pos[1]) 
    zc = center_pos[2]
    center = np.array([xc, yc, zc])

    # Define cylindrical tile regarding to its local spherical coordinate system
    ver_cyl = np.array([[(r - dr/2) * math.cos(theta - dtheta/2), (r - dr/2) * math.sin(theta - dtheta/2), z - dz/2],\
        [(r + dr/2) * math.cos(theta - dtheta/2), (r + dr/2) * math.sin(theta - dtheta/2), z - dz/2],\
        [(r - dr/2) * math.cos(theta + dtheta/2), (r - dr/2) * math.sin(theta + dtheta/2), z - dz/2],\
        [(r + dr/2) * math.cos(theta + dtheta/2), (r + dr/2) * math.sin(theta + dtheta/2), z - dz/2],\
        [(r - dr/2) * math.cos(theta - dtheta/2), (r - dr/2) * math.sin(theta - dtheta/2), z + dz/2],\
        [(r + dr/2) * math.cos(theta - dtheta/2), (r + dr/2) * math.sin(theta - dtheta/2), z + dz/2],\
        [(r - dr/2) * math.cos(theta + dtheta/2), (r - dr/2) * math.sin(theta + dtheta/2), z + dz/2],\
        [(r + dr/2) * math.cos(theta + dtheta/2), (r + dr/2) * math.sin(theta + dtheta/2), z + dz/2]])
    
    # Add rotation
    R = get_rotmat([rotation[0], rotation[1], 0])
    ver_cyl = (np.dot(R, ver_cyl.T)).T
    center = (np.dot(R, center.T)).T
    
    # Creating curves
    seg_theta = np.linspace(theta - dtheta/2, theta + dtheta/2, resolution)
    seg_curves = np.zeros(shape=(4,3,resolution))
    count = 0
    for radius in [r - dr/2, r + dr/2]:
        seg_x = radius * np.cos(seg_theta)
        seg_y = radius * np.sin(seg_theta)

        for height in [z - dz/2, z + dz/2]:
            seg_z = np.linspace(height, height, resolution)
            curve = np.asarray([seg_x, seg_y, seg_z])

            for k in range(curve.shape[1]):
                curve[:,k] = (np.dot(R, curve[:,k].T)).T + offset

            seg_curves[count] = curve
            count = count + 1
    
    # Moving tile with offset
    ver_cyl = ver_cyl + offset

    # Corner points
    # ax.plot(ver_cyl[:,0], ver_cyl[:,1], ver_cyl[:,2], 'ro')
    
    # Connecting lines
    fac = np.array([[0,1], [2,3], [0,4], [1,5], [2,6], [6,7], [4,5], [3,7]])
    lines = ver_cyl[fac]
    for line in lines:
        ax.plot(line[:,0], line[:,1], line[:,2], color=color)

    # Define the faces of the cylinder
    fac = np.array([[0, 4, 5, 1], [2, 6, 7, 3]])
    surfaces = ver_cyl[fac]
    # Plot rectangular sides
    ax.add_collection3d(Poly3DCollection(surfaces, facecolors=color, linewidths=1, edgecolors=color, alpha=.25))
    # Plot curves
    for seg_curve in seg_curves:
        ax.plot(seg_curve[0], seg_curve[1], seg_curve[2], color=color)
    # Plot curved surfaces
    curved_surfaces = np.zeros(shape=(4*(resolution-1),4,3))
    comb = np.array([[0, 1], [1, 3], [2, 3], [0, 2]])
    for j in range(comb.shape[0]):   
        for i in range(resolution-1):        
            curved_surfaces[i+(resolution-1)*j] = [seg_curves[comb[j][0]][:,i], seg_curves[comb[j][1]][:,i], seg_curves[comb[j][1]][:,i+1], seg_curves[comb[j][0]][:,i+1]]
    ax.add_collection3d(Poly3DCollection(curved_surfaces, facecolors=color, linewidths=1, alpha=.25))

    # Plot vector of magnetization in the center of the cube
    ax.quiver(offset[0] + center[0], offset[1] + center[1], offset[2] + center[2], M[0], M[1], M[2], color=color, length=np.linalg.norm(dr)/2, pivot='middle', normalize=True)

def plot_circpiece(axes, center_pos, dev_center, offset, rotation, M, color, inv=False):
    ax = axes
    resolution = 100
    r, theta, z = center_pos[:]
    dr, dtheta, dz = dev_center[:]
    
    # Calculating radius of inner edge
    p1 = np.array([r * math.cos(theta - dtheta/2), r * math.sin(theta - dtheta/2)])
    p2 = np.array([r * math.cos(theta + dtheta/2), r * math.sin(theta + dtheta/2)])
    r_mid = np.linalg.norm((p1 + p2)/2)
    
    # Difference between circpiece and circpiece_inv
    if inv is False:
        r_inner = r_mid - np.linalg.norm(p2-p1)/2
    else:
        r_inner = r_mid + np.linalg.norm(p2-p1)/2

    xc = (r_inner + (r-r_inner)/2)*math.cos(center_pos[1])
    yc = (r_inner + (r-r_inner)/2)*math.sin(center_pos[1]) 
    zc = center_pos[2]
    center = np.array([xc, yc, zc])


    # Define circ piece regarding to its local spherical coordinate system
    ver_cyl = np.array([[r * math.cos(theta - dtheta/2), r * math.sin(theta - dtheta/2), z - dz/2],\
        [r * math.cos(theta + dtheta/2), r * math.sin(theta + dtheta/2), z - dz/2],\
        [r * math.cos(theta - dtheta/2), r * math.sin(theta - dtheta/2), z + dz/2],\
        [r * math.cos(theta + dtheta/2), r * math.sin(theta + dtheta/2), z + dz/2],\
        [r_inner * math.cos(theta), r_inner * math.sin(theta), z + dz/2],\
        [r_inner * math.cos(theta), r_inner * math.sin(theta), z - dz/2]])
    
    # Add rotation
    R = get_rotmat([rotation[0], rotation[1], 0])
    ver_cyl = (np.dot(R, ver_cyl.T)).T
    center = (np.dot(R, center.T)).T
    
    # Creating curves
    seg_theta = np.linspace(theta - dtheta/2, theta + dtheta/2, resolution)
    seg_curves = np.zeros(shape=(4,3,resolution))
    count = 0
    for height in [z - dz/2, z + dz/2]:
        seg_x = r * np.cos(seg_theta)
        seg_y = r * np.sin(seg_theta)
        seg_z = np.linspace(height, height, resolution)
        curve = np.asarray([seg_x, seg_y, seg_z])

        for k in range(curve.shape[1]):
            curve[:,k] = (np.dot(R, curve[:,k].T)).T + offset

        seg_curves[count] = curve
        count = count + 1
    
    # Moving tile with offset
    ver_cyl = ver_cyl + offset

    # Corner points
    # ax.plot(ver_cyl[:,0], ver_cyl[:,1], ver_cyl[:,2], 'ro')

    # Define the faces of the cylinder
    fac = np.array([[0, 5, 4, 2], [1, 5, 4, 3]])
    surfaces = ver_cyl[fac]
    
    # Plot rectangular sides
    ax.add_collection3d(Poly3DCollection(surfaces, facecolors=color, linewidths=1, edgecolors=color, alpha=.25))
    
    # Plot curves
    for seg_curve in seg_curves:
        ax.plot(seg_curve[0], seg_curve[1], seg_curve[2], color=color)
    
    # Plot curved surfaces
    curved_surfaces = np.zeros(shape=(resolution-1,4,3))
    for i in range(resolution-1):        
        curved_surfaces[i] = [seg_curves[0][:,i], seg_curves[1][:,i], seg_curves[1][:,i+1], seg_curves[0][:,i+1]]
    ax.add_collection3d(Poly3DCollection(curved_surfaces, facecolors=color, linewidths=1, alpha=.25))

    # Plot triangular
    triangle_surfaces = np.zeros(shape=(2*(resolution-1),3,3))
    for j in range(2):   
        for i in range(resolution-1):        
            triangle_surfaces[i+(resolution-1)*j] = [seg_curves[j][:,i], ver_cyl[5-j], seg_curves[j][:,i+1]]
    ax.add_collection3d(Poly3DCollection(triangle_surfaces, facecolors=color, linewidths=1, alpha=.25))

    # Plot vector of magnetization in the center of the cube
    ax.quiver(offset[0] + center[0], offset[1] + center[1], offset[2] + center[2], M[0], M[1], M[2], color=color, length=np.linalg.norm(r-r_inner)/2, pivot='middle', normalize=True)


def get_rotmat(rot):
    rot_x = [1, 0, 0], [0, math.cos(rot[0]), -math.sin(rot[0])], [0, math.sin(rot[0]), math.cos(rot[0])]
    rot_y = [math.cos(rot[1]), 0, math.sin(rot[1])], [0, 1, 0], [-math.sin(rot[1]), 0, math.cos(rot[1])]
    rot_z = [math.cos(rot[2]), -math.sin(rot[2]), 0], [math.sin(rot[2]), math.cos(rot[2]), 0], [0, 0, 1]
    R = np.matmul(np.matmul(rot_z, rot_y), rot_x)
    # cx,cy,cz = np.cos(theta) 
    # sx,sy,sz = np.sin(theta)
    # R = np.zeros((3,3))
    # R.flat = (cx*cz - sx*cy*sz, cx*sz + sx*cy*cz, sx*sy,
    #     -sx*cz - cx*cy*sz, -sx*sz + cx*cy*cz, 
    #     cx*sy, sy*sz, -sy*cz, cy)
    return R 

def plot_field(axes, points, H):
    ax = axes
    cmap = cm.get_cmap('Blues')
    # Max color is set to norm of 30000
    norm = colors.Normalize(vmin=0, vmax=30000)
    # Set length of field arrows depending on number of evaluated points
    if len(points) > 1000:
        len_arrow = 0.01
    elif len(points) > 500:
        len_arrow = 0.025
    elif len(points) > 100:
        len_arrow = 0.05
    else:
        len_arrow = 0.075
    for i, point in enumerate(points):
        ax.quiver(point[0], point[1], point[2], H[i][0], H[i][1], H[i][2], colors=cmap(norm(np.linalg.norm(H[i]))),
        pivot='middle', length=len_arrow, normalize=True)
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)


def plot_grid(ax, grid):
    # Draw grid lines
    ax.grid(False)
    x = np.linspace(0, grid.area[0], grid.places[0]+1)
    y = np.linspace(0, grid.area[1], grid.places[1]+1)
    z = np.linspace(0, grid.area[2], grid.places[2]+1)
    segs = []
    for mark_x in x:
        segs.append(np.array([[mark_x, 0, 0], [mark_x, grid.area[1], 0]]))
        segs.append(np.array([[mark_x, 1, 0], [mark_x, 1, grid.area[2]]]))
    for mark_y in y:
        segs.append(np.array([[0, mark_y, 0], [grid.area[0], mark_y, 0]]))
        segs.append(np.array([[1, mark_y, 0], [1, mark_y, grid.area[2]]]))
    for mark_z in z:
        if mark_z == 0:
            pass
        segs.append(np.array([[1, 0, mark_z], [1, 1, mark_z]]))
        segs.append(np.array([[1, 0, mark_z], [1, 1, mark_z]]))
        segs.append(np.array([[0, 1, mark_z], [1, 1, mark_z]]))
        segs.append(np.array([[0, 1, mark_z], [1, 1, mark_z]]))
    segs.append(np.array([[1, 1, grid.area[2]], [0, 1, grid.area[2]]]))
    line_segments = Line3DCollection(segs, colors='0.75')
    ax.add_collection(line_segments)


def zoom_factory(ax, data_xlim, data_ylim, data_zlim, data_lim_range, scale=0.1):
    def zoom_scroll(event):
        # get the current x and y limits
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_zlim = ax.get_zlim()
        cur_lim_range = max(cur_xlim[1] - cur_xlim[0], cur_ylim[1] - cur_ylim[0], cur_zlim[1] - cur_zlim[0])
        cur_center_x = sum(cur_xlim[:])/2
        cur_center_y = sum(cur_ylim[:])/2
        if event.button == 'up':
            # deal with zoom in
            cur_lim_range = max(cur_lim_range - scale, scale/5)
        elif event.button == 'down':
            # deal with zoom out
            cur_lim_range = min(cur_lim_range + scale, data_lim_range)
        else:
            # deal with something that should never happen
            print(event.button)
        # set new limits
        ax.set_xlim(cur_center_x - cur_lim_range/2, cur_center_x + cur_lim_range/2)
        ax.set_ylim(cur_center_y - cur_lim_range/2, cur_center_y + cur_lim_range/2)
        ax.set_zlim([cur_zlim[0], cur_zlim[0] + cur_lim_range])
        plt.draw()  # force re-draw

    def zoom_onpress(event):
        if event.button != 1:
            return
        mouse_x = event.xdata
        mouse_y = event.ydata
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_zlim = ax.get_zlim()
        cur_lim_range = max(
            cur_xlim[1] - cur_xlim[0], cur_ylim[1] - cur_ylim[0], cur_zlim[1] - cur_zlim[0])
        if mouse_x is None or mouse_y is None:
            ax.set_xlim([data_xlim[0], data_xlim[0] + data_lim_range])
            ax.set_ylim([data_ylim[0], data_ylim[0] + data_lim_range])
            ax.set_zlim([data_zlim[0], data_zlim[0] + data_lim_range])
        else:
            # format_coord function inside site-packages\mpl_toolkits\mplot3d\axes3d.py adjusted
            s = ax.format_coord(event.xdata, event.ydata)
            x = float(s[s.find('x')+2:s.find('y')-2])
            y = float(s[s.find('y')+2:s.find('z')-2])
            # z = float(s[s.find('z')+2:len(s)])
            if data_xlim[0] < x < data_xlim[1] and data_ylim[0] < y < data_ylim[1]:
                ax.set_xlim(x - cur_lim_range/2, x + cur_lim_range/2)
                ax.set_ylim(y - cur_lim_range/2, y + cur_lim_range/2)
                ax.set_zlim([data_zlim[0], data_zlim[0] + cur_lim_range])
        plt.draw()  # force re-draw

    fig = ax.get_figure()  # get the figure of interest
    # attach the call back
    fig.canvas.mpl_connect('scroll_event', zoom_scroll)
    fig.canvas.mpl_connect('button_press_event', zoom_onpress)
    # return the function
    return (zoom_scroll, zoom_onpress)


def create_plot(tiles, eval_points, H, grid=None):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    if tiles is not None:
        # Plotting for MagTenseStandalone
        if isinstance(tiles, list):
            for tile in tiles:
                if (tile.get_tile_type() == 2):
                    plot_cube(ax, tile.get_size(), tile.get_offset(), tile.get_rotation(), tile.get_M(), tile.get_color())
        # Plotting for MagTense
        else:
            for i in range(tiles.get_n()):
                # 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv, 10 = ellipsoid
                if (tiles.get_tile_type(i) == 1):
                    plot_cylindrical(ax, tiles.get_center_pos(i), tiles.get_dev_center(i), tiles.get_offset(i), tiles.get_rotation(i), tiles.get_M(i), tiles.get_color(i))
                elif (tiles.get_tile_type(i) == 2):
                    plot_cube(ax, tiles.get_size(i), tiles.get_offset(i), tiles.get_rotation(i), tiles.get_M(i), tiles.get_color(i))   
                elif (tiles.get_tile_type(i) == 3):
                    plot_circpiece(ax, tiles.get_center_pos(i), tiles.get_dev_center(i), tiles.get_offset(i), tiles.get_rotation(i), tiles.get_M(i), tiles.get_color(i))
                elif (tiles.get_tile_type(i) == 4):
                    plot_circpiece(ax, tiles.get_center_pos(i), tiles.get_dev_center(i), tiles.get_offset(i), tiles.get_rotation(i), tiles.get_M(i), tiles.get_color(i), inv=True)
                elif (tiles.get_tile_type(i) == 5):
                    pass
                else:
                    print("Tile type not supported!")

    if eval_points is not None and H is not None:
        plot_field(ax, eval_points, H)

    if grid is not None:
        plot_grid(ax, grid)
    # Workaround added get_proj function inside site-packages\mpl_toolkits\mplot3d\axes3d.py
    # Setting length of each axis individually - currrently not support for axes3d
    # https://stackoverflow.com/questions/10326371/setting-aspect-ratio-of-3d-plot
    # xmin, xmax = np.divide(self.get_xlim3d(), self.pbaspect[0])
    # ymin, ymax = np.divide(self.get_ylim3d(), self.pbaspect[1])
    # zmin, zmax = np.divide(self.get_zlim3d(), self.pbaspect[2])
    # ax.pbaspect = [ax.get_xlim()[1], ax.get_ylim()[1], ax.get_zlim()[1]]

    # Get limits from plotted data
    data_xlim = ax.get_xlim()
    data_ylim = ax.get_ylim()
    data_zlim = ax.get_zlim()
    # Scaling axis equally
    data_lim_range = max(data_xlim[1] - data_xlim[0], data_ylim[1] - data_ylim[0], data_zlim[1] - data_zlim[0])
    ax.set_xlim([data_xlim[0], data_xlim[0] + data_lim_range])
    ax.set_ylim([data_ylim[0], data_ylim[0] + data_lim_range])
    ax.set_zlim([data_zlim[0], data_zlim[0] + data_lim_range])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    zoom_factory(ax, data_xlim, data_ylim, data_zlim, data_lim_range)
    plt.show()
