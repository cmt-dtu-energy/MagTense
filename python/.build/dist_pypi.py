import argparse
import itertools
import subprocess
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
    parser.add_argument(
        "--pkg_version",
        type=str,
        default="1.0.0",
        help="Package version",
    )
    return parser.parse_args()


def main(
    py_versions: list[str],
    cu_versions: list[str],
    platforms: list[str],
    pkg_version: str,
    build_tag: dict | None = None,
) -> None:
    if build_tag is None:
        build_tag = {"cpu": 0, "cu12": 1}
    py_folder = Path(__file__).parent.parent
    lib_folder = py_folder / "src" / "magtense" / "lib"

    for platform in platforms:
        suffix = "pyd" if platform == "win" else "so"
        arch = "win_amd64" if platform == "win" else "x86_64-linux-gnu"
        whl_arch = "win_amd64" if platform == "win" else "manylinux1_x86_64"

        for lib_file in lib_folder.glob(f"*.{suffix}"):
            subprocess.run(["rm", lib_file], check=False)

        for cuda, py in itertools.product(cu_versions, py_versions):
            py_lib = "cp" + py if platform == "win" else "cpython-" + py
            ext_name = f"magtensesource.{py_lib}-{arch}.{suffix}"
            src_ext = Path(f"{py_folder}/{cuda}_libs/{ext_name}")
            if not src_ext.exists():
                raise FileNotFoundError(
                    f"Expected compiled extension not found: {src_ext}\n"
                    f"Make sure the {cuda} build step has staged the extension into "
                    f"python/{cuda}_libs/ before running this script."
                )
            subprocess.run(
                [
                    "cp",
                    str(src_ext),
                    lib_folder,
                ],
                check=True,
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
                        f"{lib_folder}/{ext_name}",
                    ],
                    check=False,
                )
            subprocess.run(
                [
                    "cp",
                    f"{py_folder}/.build/requirements-py{py[0]}.txt",
                    f"{py_folder}/requirements.txt",
                ],
                check=False,
            )
            if cuda == "cpu":
                subprocess.run(
                    [
                        "sed",
                        "-i",
                        "'/^nvidia-/d'",
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
                    f"{py_folder}/dist/magtense-{pkg_version.removeprefix('v')}-py{py[0]}-none-any.whl",
                    f"{py_folder}/dist/magtense-{pkg_version.removeprefix('v')}-{build_tag[cuda]}-py{py}-none-{whl_arch}.whl",
                ],
                check=False,
            )

            subprocess.run(
                ["rm", f"{lib_folder}/{ext_name}"],
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
    main(
        args.py_version.split(","),
        args.cu_version.split(","),
        args.platform.split(","),
        args.pkg_version,
    )
