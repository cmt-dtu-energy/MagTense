import os
import numpy as np

# Load reference points from COMSOL calculation
def load_COMSOL_eval(file_name, eval_offset, COMSOL_eval_path, pts_special=None):
    with open(os.path.join(COMSOL_eval_path, file_name), "r") as file:
                T = file.readlines()[8:]

    T_split = np.asarray([line.split() for line in T], dtype=np.float64)
    H_norm_COMSOL = T_split[:,1]

    # Preparing points along given axis
    if pts_special is None:
        pts_coor = T_split[:,0]
    else:
        pts_coor = pts_special
    struc = np.ones(len(pts_coor))
    if file_name[-5] == 'x':
        pts = np.c_[pts_coor, struc*eval_offset[1], struc*eval_offset[2]]
    elif file_name[-5] == 'y':
        pts = np.c_[struc*eval_offset[0], pts_coor, struc*eval_offset[2]]
    elif file_name[-5] == 'z':
        pts = np.c_[struc*eval_offset[0], struc*eval_offset[1], pts_coor]

    return (pts, H_norm_COMSOL)

# Adding subplot to figure
def add_subplot(ax, x, x_label, y_MagTense, y_COMSOL, plot_COMSOL=True, plot_error=False):    
    ax.plot(x, y_MagTense, 'r*', label='MagTense')
    if plot_COMSOL is True:
        ax.plot(x, y_COMSOL, 'bx', label='COMSOL')
    if plot_error is True:
        error = abs(y_COMSOL - y_MagTense)
        ax.plot(x, error, 'g', label= 'Error')
    ax.legend()
    ax.set_xlabel(x_label)
    ax.set_ylabel('H_norm')