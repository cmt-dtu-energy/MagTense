#%%
import sys
import h5py
import numpy as np

from pathlib import Path
from typing import List, Optional
from datetime import datetime
from multiprocessing import Pool, Lock, Value, cpu_count
from magtense.magstatics import Tiles, run_simulation
from magtense.utils import plot_magfield


def gen_mag_sample(
    n_samples: int,
    data_range: List,
    pts_eval: np.ndarray,
    param_mat: np.ndarray,
    res: int,
    num_mag: int,
    n_params: int,
    dim: int,
    datapath: str,
    t_start: datetime,
    M_fixed: bool,
    shared: bool = False,
    check: bool = False,
    min_pos:float = 0,
    min_magnet_size: float = 0.5,
    max_magnet_size: float = 2,
    height: float = 1,
    area_x: float = 10,
): 
    sample_dict = {}

    for idx in range(data_range[0], data_range[1]): 
        pos = np.zeros(shape=(num_mag,3))
        size = np.zeros(shape=(num_mag,3))
        M_rem = np.zeros(shape=(num_mag,))
        mag_angle = np.zeros(shape=(num_mag,2))
        n_mag = num_mag

        for i in range(num_mag):
            glob_idx = idx * num_mag + i
            pos[i] = [param_mat[glob_idx,0] * (area_x - 2 * min_pos) + min_pos,
                      param_mat[glob_idx,1] * (area_x - 2 * min_pos) + min_pos,
                      height/2]
            size_x = param_mat[glob_idx,2] * (max_magnet_size - min_magnet_size) + min_magnet_size
            size[i] = [size_x, size_x, height]
            M_rem[i] = 1.28 if M_fixed else param_mat[glob_idx,3] + 0.5
            polar_angle = param_mat[glob_idx,-1] * np.pi if dim == 3 else np.pi/2
            mag_angle[i] = [param_mat[glob_idx,-2] * 2 * np.pi, polar_angle]

        overlap = True
        while overlap:
            overlap = False
            for i in range(pos.shape[0] - 1):
                for j in range(i + 1, pos.shape[0]):
                    d = np.sqrt(((pos[i,:2] - pos[j,:2])**2).sum())
                    if d < np.sqrt(2)*((size[i][0] + size[j][0]) / 2):
                        overlap = True
                        pos = np.delete(pos, j, 0)
                        size = np.delete(size, j, 0)
                        M_rem = np.delete(M_rem, j, 0)
                        mag_angle = np.delete(mag_angle, j, 0)
                        n_mag -= 1

        tile = Tiles(
            n=pos.shape[0],
            tile_type=2,
            offset=pos,
            size=size,
            M_rem=M_rem/(4*np.pi*1e-7),
            mag_angle=mag_angle, 
            color=[1,0,0]
        )
        _, H_out = run_simulation(tile, pts_eval)
        
        field = H_out.reshape((res,res,-1)).transpose((2,0,1)) * 4 * np.pi * 1e-7

        sample = np.zeros(shape=(num_mag, n_params), dtype=np.float32)
        for i in range(pos.shape[0]):
            sample[i,0] = pos[i][0]
            sample[i,1] = pos[i][1]
            sample[i,2] = size[i][0]
            sample[i,3] = M_rem[i]
            sample[i,4] = mag_angle[i][0]
            if dim == 3: sample[i,5] = mag_angle[i][1]

        if check and idx < 10: plot_magfield(field, magnet=sample, vmax=0.75)
        sample_dict[idx] = (field, sample, n_mag)

        if shared:
            with count.get_lock():
                count.value += 1
        k = count.value if shared else idx + 1
        
        # Progress bar
        dt = (datetime.utcnow() - t_start).total_seconds()
        sys.stdout.write("\r[PROGRESS] [%-50s] %d%%  |  %d of %d  |  %.2f samples/s" \
                         % ('='*int(50*(k/n_samples)), 100*(k/n_samples), k, n_samples, k/dt))
        sys.stdout.flush()

    # Saving sample
    if shared: lock.acquire()
    db = h5py.File(f"{datapath}/{n_samples}.h5", mode='a')
    for key in sample_dict:
        db['fields'][key] = sample_dict[key][0][:dim]
        db['samples'][key] = sample_dict[key][1]
        db['info'][sample_dict[key][2]] += 1
    db.close()
    if shared: lock.release()


def init(l, cnt):
    global lock
    global count
    lock = l
    count = cnt


def create_single_mags(
    n_samples: int,
    res: int,
    num_mag: int,
    dim: int = 2,
    M_fixed: bool = False,
    check: bool = False,
    n_proc: Optional[int] = None,
    seed: int = 42
):  
    datapath = Path(__file__).parent.resolve() / '..' / 'data' / f'{num_mag}_{res}'
    if M_fixed: datapath += '_M_fixed'
    if not datapath.exists(): datapath.mkdir(parents=True)
    if Path(f'{datapath}/{n_samples}.h5').is_file():
        input(f'Overwriting file "{num_mag}_{res}/{n_samples}.h5". Press Enter to continue')
    
    worker = cpu_count() if n_proc is None else n_proc
    n_params = 6 if dim == 3 else 5
    if M_fixed: n_params -= 1
    t_start = datetime.utcnow()

    rng = np.random.default_rng(seed)
    param_mat = rng.random(size=(n_samples * num_mag, n_params))

    height = 1
    area_x = 10
    x_eval = np.linspace(0, area_x, res) + (area_x / (2 * res))
    xv, yv = np.meshgrid(x_eval, x_eval)
    zv = np.zeros(res * res) + height / 2
    pts_eval = np.hstack([xv.reshape(-1,1), yv.reshape(-1,1), zv.reshape(-1,1)])

    print(f'[INFO {t_start.strftime("%d/%m %H:%M:%S")}] #Data: {n_samples} | #Worker: {worker} | #Path: {datapath}')

    db = h5py.File(f"{datapath}/{n_samples}.h5", libver="latest", mode='w')
    db.create_dataset("fields", shape=(n_samples, dim, res, res), dtype="float32")
    db.create_dataset("samples", shape=(n_samples, num_mag, n_params), dtype="float32")
    db.create_dataset("info", shape=(num_mag + 1,), dtype="float32")
    db.close()   
    
    if n_proc == 1:
        for idx in range(n_samples):
            gen_mag_sample(n_samples, [idx, idx+1], pts_eval, param_mat, res, num_mag,
                           n_params, dim, datapath, t_start, M_fixed, check=check)

    else:
        for i in reversed(range(min(100, int(n_samples / worker)) + 1)):
            if i == 0:
                print('[WARNING] Using batch size of 1. Might be slow due to I/O with hard disk.')
                batch_size = 1
            elif n_samples % i == 0:
                batch_size = i
                break

        l = Lock()
        cnt = Value('i', 0)
        with Pool(n_proc, initializer=init, initargs=(l, cnt,)) as p:
            p.starmap(gen_mag_sample, [(n_samples, [idx*batch_size, (idx+1)*batch_size], pts_eval, param_mat,
                                       res, num_mag, n_params, dim, datapath, t_start, M_fixed, True)
                                       for idx in range(n_samples // batch_size)])

    with h5py.File(f"{datapath}/{n_samples}.h5", mode='r') as db:
        n_mag_info = ''.join([f' #{i}:{int(n)}' for i, n in enumerate(db['info'])])
        print(f'\n[INFO {datetime.utcnow().strftime("%d/%m %H:%M:%S")}] Database created | n_mag {n_mag_info}')

# %%
