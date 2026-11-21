"""The compiled MagTense extension module, magtensesource.

On Windows the extension depends on DLLs that are installed outside this
folder: the Intel Fortran runtime and oneMKL in <prefix>/Library/bin and, for
the CUDA build, the CUDA libraries of the nvidia wheels. Their folders are
registered here, so that every import of magtense.lib.magtensesource finds
them, also one that does not go through magtense.micromag or
magtense.magstatics.
"""

import os
import sys
from pathlib import Path

# Windows only
if hasattr(os, "add_dll_directory"):
    package_path = Path(__file__).parent.parent

    # First entry is the installed layout
    # (<prefix>/Lib/site-packages/magtense/../../../Library/bin), the second the
    # active environment prefix, which is what applies when running from source.
    dll_paths = [
        package_path / ".." / ".." / ".." / "Library" / "bin",
        Path(sys.prefix) / "Library" / "bin",
    ]

    # The CUDA wheels changed layout between the two major versions: cu12 gave
    # every library its own nvidia/<name>/bin, cu13 puts them all together in
    # nvidia/cu13/bin/x86_64. Offer both so either generation of wheel resolves.
    nvidia_path = package_path / ".." / "nvidia"
    dll_paths.append(nvidia_path / "cu13" / "bin" / "x86_64")
    dll_paths += [
        nvidia_path / lib / "bin"
        for lib in ["cublas", "cuda_runtime", "cusparse", "nvjitlink"]
    ]

    for dll_path in dll_paths:
        if Path.is_dir(dll_path):
            os.add_dll_directory(dll_path)
            # libifcoremd.dll resolves its own dependency on libmmd.dll through
            # PATH, which add_dll_directory does not cover.
            os.environ["PATH"] = f"{dll_path.resolve()}{os.pathsep}{os.environ['PATH']}"
