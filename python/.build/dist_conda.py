from pathlib import Path
import os
import platform
import subprocess
import shutil
import tomllib


def main(py_versions=["3.12"]):
    with open(Path(__file__).parent.parent / "pyproject.toml", "rb") as f:
        mt_version = tomllib.load(f)["project"]["version"]
    if platform.system() == "Windows":
        arch = "win-64"
    else:
        arch = "linux-64"

    for py in py_versions:
        for cuda in ["cuda12"]:
            with open("src/magtense/lib/Makefile", "r") as source_file:
                lines_to_append = source_file.readlines()[4:]

            with open("scripts/temp.txt", "w") as file:
                file.write("#==========================\n")
                file.write("# Compiler names and flags\n")
                file.write("#==========================\n")

                if cuda == "cpu":
                    file.write("USE_CUDA = 0\n")
                else:
                    file.write("USE_CUDA = 1\n")

                file.writelines(lines_to_append)

            shutil.move("scripts/temp.txt", "src/magtense/lib/Makefile")

            with open("../.conda-recipe/meta.yaml", "r") as source_file:
                lines_to_append = source_file.readlines()[3:]

            with open("scripts/temp.txt", "w") as file:
                file.write(f'{{% set version = "{mt_version}" %}}' + "\n")
                file.write(f'{{% set py_version = "{py}" %}}' + "\n")
                file.write(f'{{% set cuda_version = "{cuda}" %}}' + "\n")
                file.writelines(lines_to_append)

            shutil.move("scripts/temp.txt", "../.conda-recipe/meta.yaml")

            subprocess.run(["conda-build", "../.conda-recipe"])

            upload_command = [
                "anaconda",
                "upload",
                "--user",
                "cmt-dtu-energy",
                f"{os.environ.get('CONDA_PREFIX')}/conda-bld/{arch}/magtense-{mt_version}-py{py}_{cuda}.tar.bz2",
            ]
            subprocess.run(upload_command)
            subprocess.run(["conda", "build", "purge"])


if __name__ == "__main__":
    main()
