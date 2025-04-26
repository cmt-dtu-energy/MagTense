import argparse
import itertools
import subprocess
import tomllib
from pathlib import Path


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Build and distribute.")
    parser.add_argument(
        "--py_version",
        type=str,
        default="312",
        help="Python versions (comma-separated)",
    )
    parser.add_argument(
        "--cu_version",
        type=str,
        default="cpu,cu12",
        help="Cuda / cpu versions (comma-separated)",
    )
    parser.add_argument(
        "--platform",
        type=str,
        default="win,linux",
        help="Platforms (comma-separated)",
    )
    return parser.parse_args()


def main(
    py_versions: list[str],
    cu_versions: list[str],
    platforms: list[str],
    build_tag: dict | None = None,
) -> None:
    if build_tag is None:
        build_tag = {"cpu": 0, "cu12": 1}
    py_folder = Path(__file__).parent.parent
    lib_folder = py_folder / "src" / "magtense" / "lib"
    with Path.open(py_folder / "pyproject.toml", "rb") as f:
        mt_version = tomllib.load(f)["project"]["version"]

    for platform in platforms:
        suffix = "pyd" if platform == "win" else "so"
        arch = "win_amd64" if platform == "win" else "x86_64-linux-gnu"
        whl_arch = "win_amd64" if platform == "win" else "manylinux1_x86_64"

        for lib_file in lib_folder.glob(f"*.{suffix}"):
            subprocess.run(["rm", lib_file], check=False)

        for cuda, py in itertools.product(cu_versions, py_versions):
            py_lib = "cp" + py if platform == "win" else "cpython-" + py
            subprocess.run(
                [
                    "cp",
                    f"{py_folder}/{cuda}_libs/magtensesource.{py_lib}-{arch}.{suffix}",
                    lib_folder,
                ],
                check=False,
            )
            if platform == "linux":
                rpath = "$ORIGIN/../../../../../lib/"
                if cuda == "cu12":
                    rpath += ":$ORIGIN/../../nvidia/cublas/lib/"
                    rpath += ":$ORIGIN/../../nvidia/cuda_runtime/lib/"
                    rpath += ":$ORIGIN/../../nvidia/cusparse/lib/"
                subprocess.run(
                    [
                        "patchelf",
                        "--force-rpath",
                        "--set-rpath",
                        f"{rpath}",
                        f"{lib_folder}/magtensesource.{py_lib}-{arch}.{suffix}",
                    ],
                    check=False,
                )
            subprocess.run(
                [
                    "cp",
                    f"{py_folder}/.build/requirements-py{py[0]}-{cuda}.txt",
                    f"{py_folder}/requirements.txt",
                ],
                check=False,
            )
            subprocess.run(
                ["python", "-m", "build", "--wheel"],
                cwd=py_folder,
                check=False,
            )
            subprocess.run(
                [
                    "mv",
                    f"{py_folder}/dist/magtense-{mt_version}-py{py[0]}-none-any.whl",
                    f"{py_folder}/dist/magtense-{mt_version}-{build_tag[cuda]}-py{py}-none-{whl_arch}.whl",
                ],
                check=False,
            )

            subprocess.run(
                ["rm", f"{lib_folder}/magtensesource.{py_lib}-{arch}.{suffix}"],
                check=False,
            )
            if Path(py_folder / "src" / "magtense.egg-info").is_dir():
                subprocess.run(
                    [
                        "rm",
                        "-r",
                        f"{py_folder}/src/magtense.egg-info/",
                        f"{py_folder}/build/",
                    ],
                    check=False,
                )


if __name__ == "__main__":
    args = parse_args()
    py_versions = args.py_version.split(",")
    cu_versions = args.cu_version.split(",")
    platforms = args.platform.split(",")
    main(py_versions, cu_versions, platforms)
