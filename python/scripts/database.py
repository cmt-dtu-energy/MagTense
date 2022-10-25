#%%
import h5py

from pathlib import Path
from typing import Optional
from multiprocessing import Process, cpu_count

from scripts.data_halbach import create_halbach, db_assign


def create_db_mp(
    target,
    assign,
    n_workers: Optional[int] = None,
    datapath: Optional[Path] = None,
    **kwargs
) -> None:
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
        kwargs['interval'] = [i*intv, end_intv]
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
            assign(db_t, Path(datapath, name))
            Path(datapath, name).unlink()
    
    print('Database created')

# %%

if __name__ == '__main__':
    db_kwargs = {
        'n_halbach': 10,
        'n_mat': 1,
        'shim_segs': 8,
        'shim_layers': 3,
        'name': 'test123'
    }
    create_db_mp(create_halbach, db_assign, n_workers=1, **db_kwargs)

