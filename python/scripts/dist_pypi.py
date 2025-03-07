import itertools
import subprocess

from pathlib import Path


def main():
    mt_version = "2.2.1"
    py_versions = ["312"]
    cu_versions = ["cpu", "cuda12"]
    lib_path = None

    build_tag = {"cpu": 0, "cuda12": 1}

    for lib_file in Path("src/magtense/lib").glob("*.pyd"):
        subprocess.run(["rm", lib_file])

    for lib_file in Path("src/magtense/lib").glob("*.so"):
        subprocess.run(["rm", lib_file])

    for platform in ["win", "linux"]:
        if platform == "win":
            suffix = "pyd"
            arch = "win_amd64"
            whl_arch = "win_amd64"
        else:
            suffix = "so"
            arch = "x86_64-linux-gnu"
            whl_arch = "manylinux1_x86_64"

        for cuda, py in itertools.product(cu_versions, py_versions):
            if platform == "win":
                py_lib = "cp" + py
            else:
                py_lib = "cpython-" + py

            subprocess.run(
                [
                    "cp",
                    f"compiled_libs/magtensesource.{py_lib}-{arch}_{cuda}.{suffix}",
                    "src/magtense/lib",
                ]
            )
            if lib_path is not None:
                subprocess.run(["rm", lib_path])
            lib_path = f"src/magtense/lib/magtensesource.{py_lib}-{arch}.{suffix}"
            subprocess.run(
                [
                    "mv",
                    f"src/magtense/lib/magtensesource.{py_lib}-{arch}_{cuda}.{suffix}",
                    lib_path,
                ]
            )
            subprocess.run(["python3", "-m", "build", "--wheel"])
            subprocess.run(
                [
                    "mv",
                    f"dist/magtense-{mt_version}-py3-none-any.whl",
                    f"dist/magtense-{mt_version}-{build_tag[cuda]}-py{py}-none-{whl_arch}.whl",
                ]
            )

            if Path("src/magtense.egg-info").is_dir():
                subprocess.run(["rm", "-r", "src/magtense.egg-info/", "build/"])


if __name__ == "__main__":
    main()
