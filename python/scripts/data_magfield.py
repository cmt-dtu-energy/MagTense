#%%
import numpy as np
import ray

from magtense.magstatics import run_simulation, grid_config
from magtense.utils import plot_magfield
from multiprocessing import cpu_count
from pathlib import Path
from asyncio import Event
from typing import Optional, Tuple

from tqdm import tqdm


@ray.remote
class ProgressBarActor:
    counter: int
    delta: int
    event: Event

    def __init__(self) -> None:
        self.counter = 0
        self.delta = 0
        self.event = Event()

    def update(self, num_items_completed: int) -> None:
        """Updates the ProgressBar with the incremental
        number of items that were just completed.
        """
        self.counter += num_items_completed
        self.delta += num_items_completed
        self.event.set()

    async def wait_for_update(self) -> Tuple[int, int]:
        """Blocking call.
        Waits until somebody calls `update`, then returns a tuple of
        the number of updates since the last call to
        `wait_for_update`, and the total number of completed items.
        """
        await self.event.wait()
        self.event.clear()
        saved_delta = self.delta
        self.delta = 0
        return saved_delta, self.counter

    def get_counter(self) -> int:
        """
        Returns the total number of complete items.
        """
        return self.counter


class ProgressBar:
    progress_actor: ray.actor.ActorHandle
    total: int
    description: str
    pbar: tqdm

    def __init__(self, total: int, description: str = ""):
        # Ray actors don't seem to play nice with mypy, generating
        # a spurious warning for the following line,
        # which we need to suppress. The code is fine.
        self.progress_actor = ProgressBarActor.remote()  # type: ignore
        self.total = total
        self.description = description

    @property
    def actor(self) -> ray.actor.ActorHandle:
        """Returns a reference to the remote `ProgressBarActor`.
        When you complete tasks, call `update` on the actor.
        """
        return self.progress_actor

    def print_until_done(self) -> None:
        """Blocking call.
        Do this after starting a series of remote Ray tasks, to which you've
        passed the actor handle. Each of them calls `update` on the actor.
        When the progress meter reaches 100%, this method returns.
        """
        pbar = tqdm(desc=self.description, total=self.total)
        while True:
            delta, counter = ray.get(self.actor.wait_for_update.remote())
            pbar.update(delta)
            if counter >= self.total:
                pbar.close()
                return


def gen_magfield_sample(
    idx: int,
    res: int,
    dim: int,
    datapath: Path,
    filled_mat: np.ndarray,
    empty_mat: np.ndarray,
    check: bool = False,
    shared: bool = False,
    pba: Optional[ProgressBarActor] = None
):
    """ 
    Generate 3-D magnetic fields of experimental setup.

    Args:
        idx: Internal iterator.
        res: Resolution of magnetic field.
        dim: Dimension of demagnetizing field strength.
        datapath: Indicates where to store the data sample.
        filled_mat: Integers, which grid positions is filled.
        empty_mat: Indicate the size of area without magnetic material.
        check: Boolean for possible visual output.
        shared: True if multiprocessing is in use.
        pba: Ray actor for progress bar.
    """
    # Omit already created files for restart
    if Path(f'{datapath}/{idx}.npy').is_file(): 
        if shared: pba.update.remote(1)
        return

    spots = filled_mat.shape[1:]
    area = [1, 1, 0.5]
    tile_size = np.asarray(area) / np.asarray(spots)
    gap = 0.05
    emp_x, emp_y = {0:[4,5], 1:[4,6], 2:[3,6], 3:[3,7]}[empty_mat[idx]]
    s_x = emp_x * tile_size[0]
    s_y = emp_y * tile_size[1]
    
    filled_pos = [[i, j, k] 
                  for i in range(spots[0])
                  for j in range(spots[1])
                  for k in range(spots[2]) 
                  if filled_mat[idx][i][j][k] == 1 and (i < emp_x or i > emp_y
                    or j < emp_x or j > emp_y or k < 2 or k > 2)]
    
    tiles, _ = grid_config(spots, area, filled_pos)

    x_eval = np.linspace(s_x + gap, s_y + gap, res)
    y_eval = np.linspace(s_x + gap, s_y + gap, res)

    if dim == 2:
        xv, yv = np.meshgrid(x_eval, y_eval)
        zv = np.zeros(res * res) + area[2] / 2
    
    elif dim == 3:
        z_eval = np.linspace(-(s_y-s_x)/res, (s_y-s_x)/res, dim) + area[2] / 2
        xv, yv, zv = np.meshgrid(x_eval, y_eval, z_eval)
    
    else:
        raise ValueError('Only 2-D and 3-D magnetic field can be generated!')
    
    pts_eval = np.hstack([xv.reshape(-1,1), yv.reshape(-1,1), zv.reshape(-1,1)])
    _, H_out = run_simulation(tiles, pts_eval)

    # Tensor image with shape CxHxWxD [T]
    field = H_out * 4 * np.pi * 1e-7
    field = field.reshape((res,res,dim,3)).transpose((3,0,1,2)) \
        if dim == 3 else field.reshape((res,res,3)).transpose((2,0,1))    

    if check and idx < 10: plot_magfield(field, vmax=0.2)
    if field is not None: np.save(f'{datapath}/{idx}.npy', field)
    if shared: pba.update.remote(1)


def create_magfield_db(
    size: int,
    name: str = '',
    res: int = 256,
    dim: int = 2,
    check: bool = False,
    worker: Optional[int] = None,
    start_idx: int = 0,
    x_places: int = 10,
    z_places: int = 5, 
    seed: int = 42
) -> None:
    rng = np.random.default_rng(seed)
    filled_mat = rng.integers(2, size=(size, x_places, x_places, z_places))
    empty_mat = rng.integers(4, size=(size,))

    datapath = Path(__file__).parent.resolve() / '..' / 'data' / f'{name}_{res}'
    if not datapath.exists(): datapath.mkdir(parents=True)
    
    if worker is None: worker = cpu_count()
    if worker > 1: ray.init(num_cpus=worker, include_dashboard=False, local_mode=False)
    print(f'[INFO] #Data: {size} | #Worker: {worker} | #Path: {datapath}')

    if worker == 1:
        for idx in range(start_idx, size):
            gen_magfield_sample(idx, res, dim, datapath, filled_mat, empty_mat, check)
    else:
        pb = ProgressBar(size - start_idx)
        actor = pb.actor
        gen_magfield_ray = ray.remote(gen_magfield_sample)
        res = [gen_magfield_ray.remote(idx, res, dim, datapath, filled_mat,
               empty_mat, check, True, actor) for idx in range(start_idx, size)]
        pb.print_until_done()
        _ = [ray.get(r) for r in res]
        ray.shutdown()

# %%
