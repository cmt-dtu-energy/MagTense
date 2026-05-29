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
        build_tag = {"cpu": 0, "cu12": 1, "cu12-fmm": 3}
    py_folder = Path(__file__).parent.parent
    lib_folder = py_folder / "src" / "magtense" / "lib"

    for platform in platforms:
        suffix = "pyd" if platform == "win" else "so"
        arch = "win_amd64" if platform == "win" else "x86_64-linux-gnu"
        whl_arch = "win_amd64" if platform == "win" else "manylinux1_x86_64"

        for lib_file in lib_folder.glob(f"*.{suffix}"):
            subprocess.run(["rm", lib_file], check=False)
        for lib_file in lib_folder.glob("*.dll"):
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
            if cuda.endswith("-fmm") and platform == "linux":
                subprocess.run(
                    [
                        "cp",
                        f"{py_folder}/{cuda}_libs/libfmm3d.so",
                        lib_folder,
                    ],
                    check=False,
                )
            if cuda.endswith("-fmm") and platform == "win":
                subprocess.run(
                    [
                        "cp",
                        f"{py_folder}/{cuda}_libs/libfmm3d.dll",
                        lib_folder,
                    ],
                    check=False,
                )
            if platform == "linux":
                rpath_entries = ["$ORIGIN/../../../../../lib/"]
                if cuda.endswith("-fmm"):
                    rpath_entries.insert(0, "$ORIGIN")
                if cuda.startswith("cu12"):
                    rpath_entries += [
                        "$ORIGIN/../../nvidia/cublas/lib/",
                        "$ORIGIN/../../nvidia/cuda_runtime/lib/",
                        "$ORIGIN/../../nvidia/cusparse/lib/",
                    ]
                subprocess.run(
                    [
                        "patchelf",
                        "--force-rpath",
                        "--set-rpath",
                        ":".join(rpath_entries),
                        f"{lib_folder}/magtensesource.{py_lib}-{arch}.{suffix}",
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
            # python -m build normalizes the version in the wheel filename
            # (e.g. lowercases local version labels like fmmWorking -> fmmworking),
            # so use glob to locate the built wheel rather than constructing its exact name.
            built_wheels = sorted(
                (py_folder / "dist").glob(f"magtense-*-py{py[0]}-none-any.whl"),
                key=lambda p: p.stat().st_mtime,
            )
            if built_wheels:
                subprocess.run(
                    [
                        "mv",
                        str(built_wheels[-1]),
                        f"{py_folder}/dist/magtense-{pkg_version.removeprefix('v')}-{build_tag[cuda]}-py{py}-none-{whl_arch}.whl",
                    ],
                    check=True,
                )

            subprocess.run(
                ["rm", f"{lib_folder}/magtensesource.{py_lib}-{arch}.{suffix}"],
                check=False,
            )
            if cuda.endswith("-fmm") and platform == "linux":
                subprocess.run(
                    ["rm", f"{lib_folder}/libfmm3d.so"],
                    check=False,
                )
            if cuda.endswith("-fmm") and platform == "win":
                subprocess.run(
                    ["rm", f"{lib_folder}/libfmm3d.dll"],
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
