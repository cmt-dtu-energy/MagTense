from scipy.io import loadmat
import numpy as np

def matlab_struct_to_dict(matobj):
    out = {}
    for name in matobj.dtype.names:
        elem = matobj[name][0, 0]
        if isinstance(elem, np.ndarray) and elem.dtype.names:
            out[name] = matlab_struct_to_dict(elem)
        else:
            out[name] = elem
    return out

def load_matlab_struct(fname):
    data = loadmat(fname) 


    mesh_cart = matlab_struct_to_dict(data["mesh_cart"])
    GridInfo  = matlab_struct_to_dict(data["GridInfo"])
    mesh_params = matlab_struct_to_dict(data["mesh_params"])    
    raw_iIn = mesh_cart["iIn"][0]   # shape (n_groups,) dtype=object
    iIn = [np.asarray(g, dtype=int).ravel() - 1 for g in raw_iIn]

    return mesh_cart, GridInfo, mesh_params, iIn