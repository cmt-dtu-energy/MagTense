#%%
import h5py
import numpy as np

from magtense.halbach import HalbachCylinder, EvaluationPoints
from magtense.magstatics import iterate_magnetization, get_demag_tensor, get_H_field
from typing import Optional, List
from pathlib import Path
from multiprocessing import Process, cpu_count

from tqdm import tqdm


def load_demag_tensor(
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
    store: bool = False
) -> np.ndarray:
    '''
    If store is set, then pre-calculation of the demagnetization tensor,
    which speeds up subsequent field calculation.
    Then, H = N * M becomes only a matrix multiplication.
    '''
    fname = f'N_{pts_eval.res[0]}_{pts_eval.res[1]}_{pts_eval.res[2]}_r_{pts_eval.r}_' \
            + f'hard_{halbach.n_layers}_{halbach.n_segs}_{halbach.n_axial}_' \
            + f'shim_{halbach.shim_layers}_{halbach.shim_segs}.npy'

    fpath = Path(__file__).parent.absolute() / '..' / 'demag' / fname

    if Path(fpath).is_file():
        demag_tensor = np.load(fpath)
    else:
        if not store:
            print('[WARNING] No demagnetization tensor found for this Halbach configuration!')
            print('[RECOMMENDATION] Run `load_demag_tensor(store=True)` first for a larger amount of samples.')
        print('Calculating demagnetization tensor...')
        demag_tensor = get_demag_tensor(halbach.tiles, pts_eval.coords)
        if store:
            if not fpath.parent.exists(): fpath.parent.mkdir()
            np.save(fpath, demag_tensor)

    return demag_tensor


def calculate_demag_field(
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
) -> np.ndarray:
    '''
    Calculation of the resulting demagnetizing field inside the Halbach cylinder bore.
    MagTense uses herefore three steps:
    (1) Demagnetization tensor N
        This might be pre-calculated from the existig structure and a full shim magnet matrix.
    (2) Iterating the magnetizations M of each single magnetic tile.
        This becomes especially computationally expensive, when a lof of shim magnets are present.
    (3) Demagnetizing field: H = N * M
    '''
    demag_tensor = load_demag_tensor(halbach, pts_eval)
    it_tiles = iterate_magnetization(halbach.tiles, max_error=1e-12, mu_r=halbach.mu_r)
    field = get_H_field(it_tiles, pts_eval.coords, demag_tensor)
    # Filter posssible nan values - TODO Check with magtense
    idx_nan = [i for i, H_vec in enumerate(field) if np.isnan(np.linalg.norm(H_vec))]
    field[idx_nan,:] = [0,0,0]

    return field


def get_p2p(H_field: np.ndarray, pts_eval: np.ndarray, log: bool = False) -> float:
    '''
    Args:
        H_field: Demagnetizing field strength.
        pts_eval: Coordinates of evaluation points.
        log: If set, then log10(p2p) is returned.

    Returns:
        Peak-to-peak (p2p) value in [parts per million].
    '''
    B_norm = [np.linalg.norm(H_vec) * 4 * np.pi * 1e-7 for i, H_vec in enumerate(H_field) 
        if (pts_eval.coords[i][0]**2 + pts_eval.coords[i][1]**2 + pts_eval.coords[i][2]**2) <= pts_eval.r**2]
    p2p = ((max(B_norm) - min(B_norm)) / max(B_norm))

    return np.log10(p2p) if log else p2p * 1e6


def get_shim_mat(symm_mat, shim_layers, shim_segs, symm) -> np.ndarray:
    n_shim = shim_layers * shim_segs
    shim_mat = np.zeros(n_shim)
    if symm:
        # Set bottom and center layers
        shim_mat[:shim_segs * (shim_layers // 2 + shim_layers % 2)] = symm_mat
        # Set top layers symmetrical to bottom layers
        for i in range(shim_layers // 2):
            shim_mat[n_shim - (i+1) * shim_segs:n_shim - i * shim_segs] \
                = symm_mat[i * shim_segs:(i+1) * shim_segs]
    else:
        shim_mat[:] = symm_mat

    return shim_mat


def evaluate_shimming(
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
    log_p2p: bool = False
) -> tuple[np.ndarray, float]:
    H_next = calculate_demag_field(halbach, pts_eval)

    B_next = H_next.reshape((*pts_eval.res,3)).transpose((3,0,1,2)) * 4 * np.pi * 1e-7
    p2p_next = get_p2p(H_next, pts_eval, log=log_p2p)

    return B_next, p2p_next


### API ###
def run_halbach_environment(
    solution: np.ndarray,
    halbach: HalbachCylinder,
    pts_eval: EvaluationPoints,
    symm: bool,
    log_p2p: bool = False,
):
    shim_mat = get_shim_mat(solution, halbach.shim_layers, halbach.shim_segs, symm)
    halbach.set_shim_matrix(shim_mat)
    B_field, p2p = evaluate_shimming(halbach, pts_eval, log_p2p=log_p2p)
    
    return B_field, p2p


def get_db_stats(db_name, no_shim=False):
    db_path = Path(__file__).parent.absolute() / '..' / 'data' / db_name
    db = h5py.File(db_path, mode='r')
    size = db['p2p_no_shim'].shape[0] if no_shim else db['action'].shape[0]
    arr = np.zeros([size, 3])

    for idx in range(size):
        field = db['field_no_shim'][idx] if no_shim else db['field'][idx]
        arr[idx] = field.mean(axis=(1,2,3))

    print(f'Mean: {arr.mean(axis=0)}')
    print(f'Std: {arr.std(axis=0)}')


def db_merge(db_name):
    datapath = Path(__file__).parent.absolute() / '..' / 'data'
    p = datapath.glob('**/*')
    filenames = [x.name for x in p if x.is_file() and x.name[:len(db_name)] == db_name and x.name[:-3] != db_name]
    
    if not Path(datapath, db_name).is_file():
        db = h5py.File(Path(datapath, filenames[0]), mode='r')
        create_halbach(
            n_halbach=db.attrs['n_halbach'],
            n_mat=db.attrs['n_mat'],
            shim_segs=db.attrs['shim_segs'],
            shim_layers=db.attrs['shim_layers'],
            field_res=db.attrs['field_res'],
            norm_var=db.attrs['norm_var'],
            mu_r=db.attrs['mu_r'],
            r_probe=db.attrs['r_probe'],
            symm=db.attrs['symmetry'],
            seed_shim=db.attrs['seed_shim'],
            seed_pertubation=db.attrs['seed_pertubation'],
            action="action" in db,
            no_shim="p2p_no_shim" in db,
            name=db_name,
            empty=True
        )
        db.close()
    
    db_final = h5py.File(f'{datapath}/{db_name}.h5', mode='a')
    for name in filenames:
        db_path = datapath / name
        db = h5py.File(db_path, mode='r')
        intv_0 = int(name[name[:name.rfind('_')].rfind('_') + 1:name.rfind('_')])
        intv_1 = int(name[name.rfind('_') + 1:-3])

        i_per = 0
        for idx_per in range(intv_0, intv_1):
            db_final['halbach_M_rem'][idx_per] = db['halbach_M_rem'][i_per]
            db_final['halbach_easy_axes'][idx_per] = db['halbach_easy_axes'][i_per]

            if "p2p_no_shim" in db_final:                
                db_final['shim_field_no_shim'][idx_per] = db['shim_field_no_shim'][i_per]
                db_final['field_no_shim'][idx_per] = db['field_no_shim'][i_per]
                db_final['p2p_no_shim'][idx_per] = db['p2p_no_shim'][i_per]
                db_final['sph_ampl_no_shim'][idx_per] = db['sph_ampl_no_shim'][i_per]
            
            for idx_shim in range(db_final['p2p'].shape[1]):
                db_final['shim_field'][idx_per, idx_shim] = db['shim_field'][i_per, idx_shim]
                db_final['field'][idx_per, idx_shim] = db['field'][i_per, idx_shim]
                db_final['shim_mat'][idx_per, idx_shim] = db['shim_mat'][i_per, idx_shim]
                db_final['p2p'][idx_per, idx_shim] = db['p2p'][i_per, idx_shim]
                db_final['sph_ampl'][idx_per, idx_shim] = db['sph_ampl'][i_per, idx_shim]

                if "action" in db_final:
                    db_final['action'][idx_per, idx_shim] = db['action'][i_per, idx_shim]
                    db_final['p2p_next'][idx_per, idx_shim] = db['p2p_next'][i_per, idx_shim]
                    db_final['sph_ampl_next'][idx_per, idx_shim] = db['sph_ampl_next'][i_per, idx_shim]

            i_per += 1

        db.close()
        db_path.unlink()

    db_final.close()


def create_halbach(
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
    interval: Optional[List] = None,
    append: bool = False,
    name: Optional[str] = None,
    empty: bool = False
):
    ''' 
    Dataset creation.
    Name: '{available spots for shim magnets}_{# halbach configurations}_{# shim matrices}.h5'
    
    Args:
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
        interval: Data range to iterate over.
        append: If set, appending to existing database.
        name: Optional filename for database.
        empty: If set, an empty database is created.
    '''
    halbach = HalbachCylinder()
    halbach.add_shim_magnets(shim_layers, shim_segs, mu_r=mu_r)
    pts_eval = EvaluationPoints(field_res, r_probe)
    load_demag_tensor(halbach, pts_eval, store=True)
    load_demag_tensor(halbach, halbach.shimming_points, store=True)
    spots = shim_segs * (shim_layers // 2 + shim_layers % 2) if symm else halbach.n_shim

    rng_shim = np.random.default_rng(seed_shim)
    symm_mat = rng_shim.choice([0, 1], size=(n_halbach, n_mat, spots))


    rng_per = np.random.default_rng(seed_pertubation)
    B_rem = rng_per.normal(
        loc=halbach.B_rem,
        scale=norm_var*halbach.B_rem,
        size=(n_halbach, halbach.n_hard_tiles)
    )

    if action:
        rng_action = np.random.default_rng(seed_shim)
        empty_spots = np.where(symm_mat == 0)
        np_action = np.zeros(n_halbach, n_mat)
        
        for i in range(n_halbach):
            for j in range(n_mat):
                if i in empty_spots[0] and j in empty_spots[1]:
                    np_action[i,j] = rng_action.integers(np.where(symm_mat[i,j] == 0)[0].shape[0])
                else:
                    np_action[i,j] = -1000
    
    datapath = Path(__file__).parent.absolute() / '..' / 'data'
    if not datapath.exists(): datapath.mkdir()
    fname = name if name is not None else f'{spots}_{n_halbach}_{n_mat}_{field_res[0]}'
    n_config = n_halbach

    if append:
        db = h5py.File(f'{datapath}/{fname}.h5', libver='latest', mode='a')
        assert db.attrs['n_halbach'] == n_halbach
        assert db.attrs['n_mat'] == n_mat
        assert db.attrs['seed_shim'] == seed_shim
        assert db.attrs['seed_pertubation'] == seed_pertubation
        assert db.attrs['shim_segs'] == shim_segs
        assert db.attrs['shim_layers'] == shim_layers
        assert db.attrs['mu_r'] == mu_r
        assert db.attrs['r_probe'] == r_probe
        assert db.attrs['norm_var'] == norm_var
        assert db.attrs['symmetry'] == symm
        assert (db.attrs['field_res'] == field_res).all()
        
        if interval is not None:
            start_idx = interval[0]
            end_idx = interval[1]

    else:
        if interval is not None:
            n_config = interval[1] - interval[0]
            fname += f'_{interval[0]}_{interval[1]}'

        start_idx = 0
        end_idx = n_config

        db = h5py.File(f'{datapath}/{fname}.h5', libver='latest', mode='w')
        db.create_dataset("halbach_M_rem", shape=(n_config, halbach.n_hard_tiles), dtype="float32")
        db.create_dataset("halbach_easy_axes", shape=(n_config, halbach.n_hard_tiles, 3), dtype="float32")
        db.create_dataset("field", shape=(n_config, n_mat, 3, field_res[0], field_res[1], field_res[2]), dtype="float32")
        db.create_dataset("shim_field", shape=(n_config, n_mat, halbach.shimming_points.coords.shape[0], 3), dtype="float32")
        db.create_dataset("shim_mat", shape=(n_config, n_mat, spots), dtype="int32")
        db['shim_mat'][:,:] = symm_mat[start_idx:end_idx]
        db.create_dataset("p2p", shape=(n_config, n_mat, 1), dtype="float64")

        db.attrs['n_halbach'] = n_halbach
        db.attrs['n_mat'] = n_mat
        db.attrs['shim_pos'] = halbach.shimming_points.coords
        db.attrs['shim_segs'] = shim_segs
        db.attrs['shim_layers'] = shim_layers
        db.attrs['mu_r'] = mu_r
        db.attrs['r_probe'] = r_probe
        db.attrs['norm_var'] = norm_var
        db.attrs['symmetry'] = symm
        db.attrs['field_res'] = field_res
        db.attrs['seed_shim'] = seed_shim
        db.attrs['seed_pertubation'] = seed_pertubation

    if no_shim and "p2p_no_shim" not in db:
        db.create_dataset("field_no_shim", shape=(n_config, 3, field_res[0], field_res[1], field_res[2]), dtype="float32")
        db.create_dataset("shim_field_no_shim", shape=(n_config, halbach.shimming_points.coords.shape[0], 3), dtype="float32")
        db.create_dataset("p2p_no_shim", shape=(n_config, 1), dtype="float64")

    if action and "action" not in db:
        db.create_dataset("action", shape=(n_config, n_mat, 1), dtype="int32")
        db['action'][:,:] = np_action[start_idx:end_idx]
        db.create_dataset("p2p_next", shape=(n_config, n_mat, 1), dtype="float64")

    if empty:
        db.close()
        return

    with tqdm(total=n_config * n_mat) as pbar:
        for idx_per in range(start_idx, end_idx):
            halbach.perturb_config(B_rem[idx_per])
            db['halbach_M_rem'][idx_per] = halbach.remanence
            db['halbach_easy_axes'][idx_per] = halbach.easy_axes

            if no_shim:
                shim_field_no_shim = calculate_demag_field(halbach, halbach.shimming_points)
                B_field_no_shim, p2p_no_shim = evaluate_shimming(halbach, pts_eval)
                
                db['shim_field_no_shim'][idx_per] = shim_field_no_shim
                db['field_no_shim'][idx_per] = B_field_no_shim
                db['p2p_no_shim'][idx_per] = p2p_no_shim
            
            for idx_shim in range(n_mat):
                shim_mat = get_shim_mat(symm_mat[idx_per, idx_shim], shim_layers, shim_segs, symm)
                halbach.set_shim_matrix(shim_mat)
                shim_field = calculate_demag_field(halbach, halbach.shimming_points)
                B_field, p2p = evaluate_shimming(halbach, pts_eval)

                db['shim_field'][idx_per, idx_shim] = shim_field
                db['field'][idx_per, idx_shim] = B_field
                db['p2p'][idx_per, idx_shim] = p2p

                if action:
                    if np_action[idx_per, idx_shim] < 0:
                        p2p_next = np_action[idx_per, idx_shim]
                        k_next_ampl = [np_action[idx_per, idx_shim] for _ in range(52)]

                    else:
                        a = np_action[idx_per, idx_shim]
                        halbach.set_one_shim_magnet(a)
                        if symm:
                            a_symm = halbach.n_shim - ((a // halbach.shim_segs + 1) * halbach.shim_segs) + (a % halbach.shim_segs)
                            halbach.set_one_shim_magnet(a_symm)
                        _, p2p_next = evaluate_shimming(halbach, pts_eval)

                    db['p2p_next'][idx_per, idx_shim] = p2p_next
                    db['sph_ampl_next'][idx_per, idx_shim] = k_next_ampl

                pbar.update(1)
    
    db.close()


def create_db_mp(n_workers=None, **kwargs):
    if n_workers is None: n_workers = cpu_count()
    intv = kwargs['n_halbach'] // n_workers
    if kwargs['n_halbach'] % n_workers > 0: intv += 1
        
    l_p = []
    for i in range(n_workers):
        end_intv = min((i+1) * intv, kwargs['n_halbach'])
        kwargs['interval'] = [i * intv, end_intv]
        p = Process(target=create_halbach, kwargs=kwargs)
        p.start()
        l_p.append(p)
        if end_intv == kwargs['n_halbach']: break

    for proc in l_p:
        proc.join()
    
    # TODO Get name of database
    if 'name' in kwargs.keys():
        db_name = kwargs['name']
    else:
        spots = kwargs['shim_segs'] * (kwargs['shim_layers'] // 2 + 1)
        db_name = f'{spots}_{kwargs["n_halbach"]}_{kwargs["n_mat"]}_{16}'
    db_merge(db_name)
    print('Database created')

#%%
