from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm


def plot_cubes(axes, cubes):
    ax = axes
    # Define the vertices of the unit cubic and move them in order to center the cube on origo
    ver = np.array([[1, 1, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1], [0, 0, 1], [1, 0, 1], [1, 0, 0], [0, 0, 0]]) - 0.5

    # Plotting for MagTenseStandalone
    if isinstance(cubes, list):
        for cube in cubes:
            ver_cube = ver * cube.get_size() + cube.get_offset()
            # Define the faces of the unit cubic
            fac = np.array([[0, 1, 2, 3], [3, 2, 4, 5], [5, 6, 7, 4],
                            [0, 1, 7, 6], [5, 6, 0, 3], [1, 2, 4, 7]])
            surfaces = ver_cube[fac]
            # plot sides
            ax.add_collection3d(Poly3DCollection(surfaces, facecolors=cube.get_color(), linewidths=1, edgecolors=cube.get_color(), alpha=.25))
    # Plotting for MagTense
    else:
        for i in range(cubes.get_n()):
            # TODO Rotate the cubes according to the angles specified in cubes.getRot()
            # Rotation about x-axis: [[1, 0, 0], [0, cos(rot_x), -sin(rot_x)], [0, sin(rot_x), cos(rot_x)]]
            # Rotation about y-axis: [[cos(rot_y), 0, sin(rot_y)], [0, 1, 0], [-sin(rot_y), 0, cos(rot_y)]]
            # Rotation about z-axis: [[cos(rot_z), -sin(rot_z), 0], [sin(rot_z), cos(rot_z), 0], [0, 0, 1]]
            ver_cube = ver * cubes.get_size(i) + cubes.get_offset(i)
            # Define the faces of the unit cubic
            fac = np.array([[0, 1, 2, 3], [3, 2, 4, 5], [5, 6, 7, 4],
                            [0, 1, 7, 6], [5, 6, 0, 3], [1, 2, 4, 7]])
            surfaces = ver_cube[fac]
            # plot sides
            ax.add_collection3d(Poly3DCollection(surfaces, facecolors=cubes.get_color(i), linewidths=1, edgecolors=cubes.get_color(i), alpha=.25))


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


def create_plot(cubes, eval_points, H, grid=None):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plot_cubes(ax, cubes)
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
