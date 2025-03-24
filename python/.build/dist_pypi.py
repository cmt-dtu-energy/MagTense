import itertools
import subprocess

from pathlib import Path
import argparse
import tomllib


def parse_args():
    parser = argparse.ArgumentParser(description="Build and distribute.")
    parser.add_argument(
        "--py_version",
        type=str,
        default="312",
        help="Python versions (comma-separated)",
    )
    parser.add_argument(
        "--platform",
        type=str,
        default="win,linux",
        help="Platforms (comma-separated)",
    )
    parser.add_argument(
        "--cvode", action=argparse.BooleanOptionalAction, help="Enable cvode"
    )
    return parser.parse_args()


def main(py_versions, platforms, cvode):
    py_folder = Path(__file__).parent.parent
    lib_folder = py_folder / "src" / "magtense" / "lib"
    with open(py_folder / "pyproject.toml", "rb") as f:
        mt_version = tomllib.load(f)["project"]["version"]
    cu_versions = ["cpu", "cu12"]
    build_tag = {"cpu": 0, "cu12": 1}

    for lib_file in lib_folder.glob("*.pyd"):
        subprocess.run(["rm", lib_file])

    for lib_file in lib_folder.glob("*.so"):
        subprocess.run(["rm", lib_file])

    subprocess.run(["rm", "-rf", f"{lib_folder}/cvode/"])
    if cvode:
        subprocess.run(["mkdir", f"{lib_folder}/cvode/"])
        subprocess.run(
            ["cp", "-r", f"{py_folder}/../cvode/lib/", f"{lib_folder}/cvode/lib/"]
        )

    for platform in platforms:
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
                ["rm", f"{lib_folder}/magtensesource.{py_lib}-{arch}.{suffix}"]
            )
            subprocess.run(
                [
                    "cp",
                    f"{py_folder}/{cuda}_libs/magtensesource.{py_lib}-{arch}.{suffix}",
                    lib_folder,
                ]
            )
            if platform == "linux":
                rpath = "$ORIGIN/../../../../../lib"
                if cvode:
                    rpath += ":$ORIGIN/../../cvode/lib/"
                if cuda == "cu12":
                    rpath += ":$ORIGIN/../../nvidia/cublas/lib/:$ORIGIN/../../nvidia/cuda_runtime/lib/:$ORIGIN/../../nvidia/cusparse/lib/"
                subprocess.run(
                    [
                        "patchelf",
                        "--force-rpath",
                        "--set-rpath",
                        f"{rpath}",
                        f"{lib_folder}/magtensesource.{py_lib}-{arch}.{suffix}",
                    ]
                )
            subprocess.run(
                [
                    "cp",
                    f"{py_folder}/.build/requirements-py{py[0]}-{cuda}.txt",
                    f"{py_folder}/requirements.txt",
                ]
            )
            subprocess.run(
                ["python3", "-m", "build", "--wheel"],
                cwd=py_folder,
            )
            subprocess.run(
                [
                    "mv",
                    f"{py_folder}/dist/magtense-{mt_version}-py{py[0]}-none-any.whl",
                    f"{py_folder}/dist/magtense-{mt_version}-{build_tag[cuda]}-py{py}-none-{whl_arch}.whl",
                ]
            )

            # if cuda == "cu12":
            #     subprocess.run(
            #         [
            #             "cp",
            #             f"{py_folder}/dist/magtense-{mt_version}-{build_tag[cuda]}-py{py}-none-{whl_arch}.whl",
            #             f"{py_folder}/dist/magtense-{mt_version}-py{py}-none-{whl_arch}.whl",
            #         ]
            #     )

            if Path(py_folder / "src" / "magtense.egg-info").is_dir():
                subprocess.run(
                    [
                        "rm",
                        "-r",
                        f"{py_folder}/src/magtense.egg-info/",
                        f"{py_folder}/build/",
                    ]
                )


if __name__ == "__main__":
    args = parse_args()
    py_versions = args.py_version.split(",")
    platforms = args.platform.split(",")
    cvode = args.cvode
    main(py_versions, platforms, cvode)
