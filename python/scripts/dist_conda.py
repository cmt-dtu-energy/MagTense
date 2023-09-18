import os
import platform
import subprocess


def main(mt_version = '2.2.1', py_versions = ['3.10', '3.11']):
    if platform.system() == 'Windows':
        arch = "win-64"
    else:
        arch = "linux-64"

    for py in py_versions:
        for cuda in ["cpu", "cuda12"]:
            if cuda == "cpu":
                s_make = '1,4c\\==========================\\n# Compiler names and flags\\==========================\\nUSE_CUDA = 0'
            else:
                s_make = '1,4c\\==========================\\n# Compiler names and flags\\n==========================\\nUSE_CUDA = 1'
            
            s_meta = f'1,3c\\{{% set version = "{mt_version}" %}}\\n{{% set py_version = "{py}" %}}\\n{{% set cuda_version = "{cuda}" %}}'

            subprocess.run(["sed", "-i", s_make, "src/magtense/lib/Makefile"])
            subprocess.run(["sed", "-i", s_meta, "../.conda-recipe/meta.yaml"])
            subprocess.run(["conda-build", "../.conda-recipe"])

            upload_command = [
                "anaconda",
                "upload",
                "--user",
                "cmt-dtu-energy",
                f"{os.environ.get('CONDA_PREFIX')}/conda-bld/{arch}/magtense-{mt_version}-py{py}_{cuda}.tar.bz2"
            ]
            subprocess.run(upload_command)
            subprocess.run(["conda", "build", "purge"])

if __name__== '__main__':
    main()