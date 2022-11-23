import platform
import itertools
import subprocess

from pathlib import Path


def main():
    mt_version = '2.1.0'
    py_versions = ['39', '310']
    cu_versions = ['cpu', 'cuda116', 'cuda117']
    lib_path = None

    for cuda, py in itertools.product(cu_versions, py_versions):                
        if platform.system() == 'Windows':
            suffix = 'pyd'
            py_lib = 'cp' + py
            arch = 'win_amd64'
            whl_arch = 'win_amd64'
        else:
            suffix = 'so'
            py_lib = 'cpython-' + py
            arch = 'x86_64-linux-gnu'
            whl_arch = 'manylinux1_x86_64'

        subprocess.run(["cp", f"magtense/compiled_libs/magtensesource.{py_lib}-{arch}_{cuda}.{suffix}", "magtense/lib"])
        if lib_path is not None: subprocess.run(["rm", lib_path])
        lib_path = f"magtense/lib/magtensesource.{py_lib}-{arch}.{suffix}"
        subprocess.run(["mv", f"MagTense/python/magtense/lib/*.{suffix}", lib_path])
        subprocess.run(["python", "setup.py", "bdist_wheel"])
        subprocess.run(["mv", f"dist/magtense-{mt_version}-py3-none-any.whl", f"dist/magtense-{mt_version}+{cuda}-py{py}-none-{whl_arch}.whl"])

        if Path('magtense.egg-info').is_dir():
            subprocess.run(["rm", "-r", "magtense.egg-info/", "build/"])

if __name__== '__main__':
    main()