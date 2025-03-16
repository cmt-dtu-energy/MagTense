import itertools
import subprocess

from pathlib import Path
import argparse
import tomllib


def parse_args():
    parser = argparse.ArgumentParser(description="Build and distribute.")
    parser.add_argument(
        "--py_versions",
        type=str,
        default="312",
        help="Python versions (comma-separated)",
    )
    return parser.parse_args()


def main(py_versions):
    with open(Path(__file__).parent.parent / "pyproject.toml", "rb") as f:
        mt_version = tomllib.load(f)["project"]["version"]
    cu_versions = ["cpu", "cu12"]
    lib_path = None

    build_tag = {"cpu": 0, "cu12": 1}

    for lib_file in Path("src/magtense/lib").glob("*.pyd"):
        subprocess.run(["rm", lib_file])

    for lib_file in Path("src/magtense/lib").glob("*.so"):
        subprocess.run(["rm", lib_file])

    for platform in ["linux"]:  # ["win", "linux"]:
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
                    f"{cuda}_libs/magtensesource.{py_lib}-{arch}.{suffix}",
                    "src/magtense/lib/",
                ]
            )
            if lib_path is not None:
                subprocess.run(["rm", lib_path])
            lib_path = f"src/magtense/lib/magtensesource.{py_lib}-{arch}.{suffix}"
            if platform == "linux":
                if cuda == "cpu":
                    subprocess.run(
                        [
                            "patchelf",
                            "--force-rpath",
                            "--set-rpath",
                            "$ORIGIN/../../../../../lib",
                            lib_path,
                        ]
                    )
                elif cuda == "cu12":
                    subprocess.run(
                        [
                            "patchelf",
                            "--force-rpath",
                            "--set-rpath",
                            "$ORIGIN/../../../../../lib:$ORIGIN/../../nvidia/cublas/lib/:$ORIGIN/../../nvidia/cuda_runtime/lib/:$ORIGIN/../../nvidia/cusparse/lib/",
                            lib_path,
                        ]
                    )
            subprocess.run(
                ["cp", f".build/requirement-{py}-{cuda}.txt", "requirements.txt"]
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
    args = parse_args()
    py_versions = args.py_versions.split(",")
    main(py_versions)
