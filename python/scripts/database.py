#%%
import h5py
import numpy as np

from pathlib import Path
from typing import Optional, List
from multiprocessing import Process, cpu_count
from magtense.halbach import HalbachCylinder, EvaluationPoints

from tqdm import tqdm

from scripts.utils_halbach import load_demag_tensor, calc_demag_field, eval_shimming, get_shim_mat
from scripts.utils_std_prob_4 import gen_s_state, gen_seq


def db_std_prob_4(
    datapath: Path,
    res: List,
    grid_size: List = [500e-9, 500e-9, 3e-9],
    n_seq: int = 4,
    t_steps: int = 500,
    t_per_step: float = 4e-12,
    h_ext_a: List = [0, 360],
    h_ext_n: List = [0, 50],
    seed: int = 0,
    intv: Optional[List] = None,
    name: Optional[str] = None,
    empty: bool = False,
    cuda: bool = False,
) -> None:
    fname = name if name is not None else f'{n_seq}_{t_steps}_{res[0]}_{res[1]}'
    if intv is None: intv = [0, n_seq]
    if not empty: fname += f'_{intv[0]}_{intv[1]}'

    db = h5py.File(f'{datapath}/{fname}.h5', 'w')
    db.create_dataset('sequence', shape=(n_seq, t_steps, 3, res[0], res[1]), dtype='float32')
    db.create_dataset('field', shape=(n_seq, 3), dtype='float32')
    if not empty: db.attrs['intv'] = intv

    if empty:
        s_state = gen_s_state(res, grid_size)
        np.save(f'{datapath}/{res[0]}_{res[1]}_s_state.npy', s_state)
        db.attrs['res'] = res
        db.attrs['grid_size'] = grid_size
        db.attrs['t_steps'] = t_steps
        db.attrs['t_per_step'] = t_per_step
        db.attrs['h_ext_angle'] = h_ext_a
        db.attrs['h_ext_norm'] = h_ext_n
        db.attrs['seed'] = seed
        db.close()
        return fname, n_seq

    rng = np.random.default_rng(seed)
    rnd_mat = rng.random(size=(n_seq, 2))

    for i in tqdm(range(intv[1] - intv[0])):
        h_ext = np.zeros(3)
        d = (h_ext_n[1] - h_ext_n[0]) * rnd_mat[i,0] + h_ext_n[0]
        theta = np.deg2rad((h_ext_a[1] - h_ext_a[0]) * rnd_mat[i,1] + h_ext_a[0])
        h_ext[:2] = d * np.array([np.cos(theta), np.sin(theta)])
        db['field'][i] = h_ext

        seq = gen_seq(
            m0_state=np.load(f'{datapath}/{res[0]}_{res[1]}_s_state.npy'),
            res=res,
            grid_size=grid_size,
            h_ext=h_ext,
            t_steps=t_steps,
            t_per_step=t_per_step,
            cuda=cuda
        )

        # Output shape: (t, comp, res_x, res_y)
        db['sequence'][i] = np.moveaxis(seq.reshape(t_steps, res[1], res[0], 3).swapaxes(1,2), -1, 1)

    db.close()


def db_halbach(
    datapath: Path,
    n_halbach: int,
    n_mat: int,
    shim_segs: int = 8,
    shim_layers: int = 3,
    field_res: List[int] = [16, 16, 16],
    norm_var: float = 0.05,
    mu_r: float = 100,
    r_probe: float = 0.002,
    symm: bool = True,
    seed_shim: int = 0,
    seed_pertubation: int = 1,
    action: bool = False,
    no_shim: bool = True,
    intv: Optional[List] = None,
    name: Optional[str] = None,
    empty: bool = False
):
    ''' 
    Dataset creation.
    Name: '{available spots for shim magnets}_{# halbach configurations}_{# shim matrices}.h5'
    
    Args:
        datapath: Path to database.
        n_halbach: Number of different Halbach configuration.
        n_mat: Number of shim matrices per Halbach configuration.
        shim_segs: Number of shim magnets per layer.
        shim_layers: Layers of shim magnets in z-direction.
        field_res: Resolution of measured field in cylindrical bore.
        norm_var: Percentage pertubation of remanent magnetization.
        mu_r: Relative permeability.
        r_probe: Defines size of cube, where the magnetic field is measured.
        symm: If set, then symmetry along xy-plane is assumed.
        seed_shim: Seed for random number generator of shim matrix.
        seed_pertubation: Seed for random number generator of pertubation.
        action: If set, action as single placement of a shim magnet is stored.
        no_shim: If set, information of unshimmed field is stored.
        intv: Data range to iterate over.
        name: Optional filename for database.
        empty: If set, an empty database is created.
    '''
    halbach = HalbachCylinder()
    halbach.add_shim_magnets(shim_layers, shim_segs, mu_r=mu_r)
    pts_eval = EvaluationPoints(field_res, r_probe)
    load_demag_tensor(halbach, pts_eval, store=True)
    load_demag_tensor(halbach, halbach.pts_shim, store=True)

    rng_shim = np.random.default_rng(seed_shim)
    spots = shim_segs * (shim_layers // 2 + shim_layers % 2) if symm else halbach.n_shim
    symm_mat = rng_shim.choice([0, 1], size=(n_halbach, n_mat, spots))

    fname = name if name is not None else f'{spots}_{n_halbach}_{n_mat}_{field_res[0]}'
    if intv is None: intv = [0, n_halbach]
    n_config = intv[1] - intv[0]
    if not empty: fname += f'_{intv[0]}_{intv[1]}'

    db = h5py.File(f'{datapath}/{fname}.h5', libver='latest', mode='w')
    db.create_dataset("halbach_M_rem", shape=(n_config, halbach.n_hard_tiles), dtype="float32")
    db.create_dataset("halbach_ea", shape=(n_config, halbach.n_hard_tiles, 3), dtype="float32")
    db.create_dataset("field", shape=(n_config, n_mat, 3, field_res[0], field_res[1], field_res[2]), dtype="float32")
    db.create_dataset("shim_field", shape=(n_config, n_mat, halbach.pts_shim.coords.shape[0], 3), dtype="float32")
    db.create_dataset("p2p", shape=(n_config, n_mat, 1), dtype="float64")
    if not empty: db.attrs['intv'] = intv

    if no_shim:
        db.create_dataset("field_no_shim", shape=(n_config, 3, field_res[0], field_res[1], field_res[2]), dtype="float32")
        db.create_dataset("shim_field_no_shim", shape=(n_config, halbach.pts_shim.coords.shape[0], 3), dtype="float32")
        db.create_dataset("p2p_no_shim", shape=(n_config, 1), dtype="float64")

    if action:
        db.create_dataset("p2p_next", shape=(n_config, n_mat, 1), dtype="float64")
        rng_action = np.random.default_rng(seed_shim)
        empty_spots = np.where(symm_mat == 0)
        np_action = np.zeros(n_halbach, n_mat)
        
        for i in range(n_halbach):
            for j in range(n_mat):
                if (i in empty_spots[0]) and (j in empty_spots[1]):
                    np_action[i,j] = rng_action.integers(np.where(symm_mat[i,j] == 0)[0].shape[0])
                else:
                    np_action[i,j] = -1000

    if empty:
        db.create_dataset("shim_mat", shape=(n_halbach, n_mat, spots), dtype="int32")
        db['shim_mat'][:] = symm_mat

        if action:
            db.create_dataset("action", shape=(n_halbach, n_mat, 1), dtype="int32")
            db['action'][:] = np_action

        db.attrs['n_halbach'] = n_halbach
        db.attrs['n_mat'] = n_mat
        db.attrs['shim_pos'] = halbach.pts_shim.coords
        db.attrs['shim_segs'] = shim_segs
        db.attrs['shim_layers'] = shim_layers
        db.attrs['mu_r'] = mu_r
        db.attrs['r_probe'] = r_probe
        db.attrs['norm_var'] = norm_var
        db.attrs['symmetry'] = symm
        db.attrs['field_res'] = field_res
        db.attrs['seed_shim'] = seed_shim
        db.attrs['seed_pertubation'] = seed_pertubation
        db.close()

        return fname, n_halbach
   
    rng_per = np.random.default_rng(seed_pertubation)
    B_rem = rng_per.normal(loc=halbach.B_rem, scale=norm_var*halbach.B_rem,
                           size=(n_halbach, halbach.n_hard_tiles))

    with tqdm(total=n_config * n_mat) as pbar:
        for idx_per in range(n_config):
            halbach.perturb_config(B_rem[intv[0] + idx_per])
            db['halbach_M_rem'][idx_per] = halbach.remanence
            db['halbach_ea'][idx_per] = halbach.easy_axes

            if no_shim:
                shim_field_no_shim = calc_demag_field(halbach, halbach.pts_shim)
                B_field_no_shim, p2p_no_shim = eval_shimming(halbach, pts_eval)
                
                db['shim_field_no_shim'][idx_per] = shim_field_no_shim
                db['field_no_shim'][idx_per] = B_field_no_shim
                db['p2p_no_shim'][idx_per] = p2p_no_shim
            
            for idx_shim in range(n_mat):
                shim_mat = get_shim_mat(symm_mat[intv[0] + idx_per, idx_shim],
                                        shim_layers, shim_segs, symm)
                halbach.set_shim_matrix(shim_mat)
                shim_field = calc_demag_field(halbach, halbach.pts_shim)
                B_field, p2p = eval_shimming(halbach, pts_eval)

                db['shim_field'][idx_per, idx_shim] = shim_field
                db['field'][idx_per, idx_shim] = B_field
                db['p2p'][idx_per, idx_shim] = p2p
                print(p2p)

                if action:
                    if np_action[intv[0] + idx_per, idx_shim] < 0:
                        p2p_next = np_action[intv[0] + idx_per, idx_shim]

                    else:
                        a = np_action[intv[0] + idx_per, idx_shim]
                        halbach.set_one_shim_magnet(a)
                        if symm:
                            a_symm = halbach.n_shim - ((a // halbach.shim_segs + 1) \
                                     * halbach.shim_segs) + (a % halbach.shim_segs)
                            halbach.set_one_shim_magnet(a_symm)
                        _, p2p_next = eval_shimming(halbach, pts_eval)

                    db['p2p_next'][idx_per, idx_shim] = p2p_next

                pbar.update(1)
    
    db.close()


def db_halbach_field_stats(db_path):
    with h5py.File(db_path, mode='r') as db:
        arr = np.mean(db[f'field'], axis=(2,3,4))

    print(f'Mean: {np.mean(arr, axis=0)}')
    print(f'Std: {np.std(arr, axis=0)}')


def create_db_mp(
    data: str,
    n_workers: Optional[int] = None,
    datapath: Optional[Path] = None,
    **kwargs
) -> None:
    if data == 'halbach':
        target = db_halbach
    elif data == 'std_prob_4':
        target = db_std_prob_4
    else:
        raise NotImplementedError()

    if datapath is None: 
        datapath = Path(__file__).parent.absolute() / '..' / 'data'
    if not datapath.exists(): datapath.mkdir(parents=True)
    kwargs['datapath'] = datapath
    
    db_name, n_tasks = target(**kwargs, empty=True)
    
    if n_workers is None: n_workers = cpu_count()
    intv = n_tasks // n_workers
    if n_tasks % n_workers > 0: intv += 1
        
    l_p = []
    for i in range(n_workers):
        end_intv = min((i+1) * intv, n_tasks)
        kwargs['intv'] = [i*intv, end_intv]
        p = Process(target=target, kwargs=kwargs)
        p.start()
        l_p.append(p)
        if end_intv == n_tasks: break

    for proc in l_p:
        proc.join()
    
    path = datapath.glob('**/*')
    fnames = [x.name for x in path if x.is_file()
              and x.name[:len(db_name)] == db_name
              and x.name[:-3] != db_name]

    with h5py.File(f'{datapath}/{db_name}.h5', mode='a') as db_t:
        for name in fnames:
            with h5py.File(Path(datapath, name), mode='r') as db_s:
                intv = db_s.attrs['intv']
                for key in db_s.keys():
                    db_t[key][intv[0]:intv[1]] = db_s[key]
            Path(datapath, name).unlink()
    
    print('Database created')

# %%

if __name__ == '__main__':
    db_kwargs = {
        'n_halbach': 5,
        'n_mat': 2,
        'shim_segs': 8,
        'shim_layers': 3,
        'name': 'test1'
    }

    create_db_mp('halbach', n_workers=1, **db_kwargs)

    # db_kwargs = {
    #     'res': [16, 4, 1],
    #     'grid_size': [500e-9, 500e-9, 3e-9],
    #     'n_seq': 4,
    # }

    # create_db_mp('std_prob_4', n_workers=1, **db_kwargs)

