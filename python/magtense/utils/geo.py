import math
import numpy as np

def euler_to_rot_axis(euler):
    '''
    Converting Euler angles to rotation axis.
    For spheroids the rotation axis has to be set rather than the specific Euler angles.
    '''
    assert np.asarray(euler).any()
    # Rotation in MagTense performed in local coordinate system:
    # (1) Rot_X_L, (2) Rot_Y_L', (3) Rot_Z_L''

    # symm-axis of geometry points in the direction of z-axis of L
    # Rotates given rotation axis with pi/2 around y_L''
    # Moves x-axis of L'' to z-axis
    # ax[0] = ax[2], ax[1] = ax[1], ax[2] = -ax[0]
    R_y = np.array([[np.cos(np.pi/2), 0, np.sin(np.pi/2)],
                    [0, 1, 0],
                    [-np.sin(np.pi/2), 0, np.cos(np.pi/2)]])
    ax = np.dot(R_y, np.asarray(euler).T)

    # Calculate the spherical coordinates: yaw and pitch
    # x_L'' has to match ax
    # Perform negative yaw around x_L and pitch around y_L'
    # The azimuthal angle is offset by pi/2 (zero position of x_L'')
    rot_x = -np.arctan2(ax[1],ax[0])
    rot_y = np.arccos(ax[2] / np.sqrt(ax[0]**2 + ax[1]**2 + ax[2]**2)) - np.pi/2
    
    return np.array([rot_x, rot_y, 0])


def get_rotmat(rot):
    rot_x = [1, 0, 0], [0, math.cos(rot[0]), -math.sin(rot[0])], [0, math.sin(rot[0]), math.cos(rot[0])]
    rot_y = [math.cos(rot[1]), 0, math.sin(rot[1])], [0, 1, 0], [-math.sin(rot[1]), 0, math.cos(rot[1])]
    rot_z = [math.cos(rot[2]), -math.sin(rot[2]), 0], [math.sin(rot[2]), math.cos(rot[2]), 0], [0, 0, 1]
    # TODO: Check rotation from local to global: (1) Rot_X, (2) Rot_Y, (3) Rot_Z
    # G to L in local coordinate system: 
    R = np.asarray(rot_x) @ np.asarray(rot_y) @ np.asarray(rot_z)

    return R
