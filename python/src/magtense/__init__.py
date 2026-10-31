import os
import sys
from importlib import metadata
from pathlib import Path

# Windows only. This has to live here, in the package's own __init__, rather
# than in a submodule like magstatics.py/micromag.py (which used to each carry
# their own copy of this block): Python only runs a module's top-level code
# when that specific module is imported, and the compiled extension
# (magtense.lib.magtensesource) can be reached directly - e.g.
# `from magtense.lib import magtensesource` - without ever touching those
# files. __init__.py is the one place guaranteed to run before any submodule
# of magtense does, on every import path.
if hasattr(os, "add_dll_directory"):
    # First entry is the installed layout
    # (<prefix>/Lib/site-packages/magtense/../../../Library/bin), the second the
    # active environment prefix, which is what applies when running from source.
    dll_paths = [
        Path(__file__).parent / ".." / ".." / ".." / "Library" / "bin",
        Path(sys.prefix) / "Library" / "bin",
    ]

    # The CUDA wheels changed layout between the two major versions: cu12 gave
    # every library its own nvidia/<name>/bin, cu13 puts them all together in
    # nvidia/cu13/bin/x86_64. Offer both so either generation of wheel resolves.
    nvidia_path = Path(__file__).parent / ".." / "nvidia"
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

if (Path(__file__).parent.parent.parent / "pyproject.toml").exists():
    # Set dynamically in .github/workflows/python-package-conda.yml
    # Fallback if not set
    v = "1.0.0"
    __version__ = v.removeprefix("v")
else:
    __version__ = metadata.version("magtense")
