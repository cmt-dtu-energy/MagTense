import os
import platform
import subprocess
import shutil


def main(mt_version = '2.2.1', py_versions = ['3.11']):
    if platform.system() == 'Windows':
        arch = "win-64"
    else:
        arch = "linux-64"

    for py in py_versions:
        for cuda in ["cuda12"]:
            with open("src/magtense/lib/Makefile", 'r') as source_file:
                lines_to_append = source_file.readlines()[5:]

            with open("scripts/temp.txt", "w") as file:
                file.write("#==========================\n")
                file.write("# Compiler names and flags\n")
                file.write("#==========================\n")

                if cuda == "cpu":
                    file.write("USE_CUDA = 0\n")
                else:
                    file.write("USE_CUDA = 1\n")

                if platform.system() == 'Windows':
                    file.write("USE_NVCC = 0\n")
                else:
                    file.write("USE_NVCC = 1\n")

                file.writelines(lines_to_append)

            shutil.move("scripts/temp.txt", "src/magtense/lib/Makefile")
            
            with open("../.conda-recipe/meta.yaml", 'r') as source_file:
                lines_to_append = source_file.readlines()[3:]

            with open("scripts/temp.txt", "w") as file:
                file.write(f'{{% set version = "{mt_version}" %}}' + '\n')
                file.write(f'{{% set py_version = "{py}" %}}' + '\n')
                file.write(f'{{% set cuda_version = "{cuda}" %}}' + '\n')
                file.writelines(lines_to_append)

            shutil.move("scripts/temp.txt", "../.conda-recipe/meta.yaml")

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