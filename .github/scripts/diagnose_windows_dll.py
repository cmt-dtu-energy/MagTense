"""Diagnose a Windows "DLL load failed" error when importing magtensesource.

Windows' DLL loader does not report which specific dependency in the chain
failed to resolve, so this narrows it down three ways:

1. Lists what DLLs are actually on disk in every directory
   ``magstatics.py``'s own ``add_dll_directory`` logic checks, to catch a
   stale/incomplete install there.
2. Retries loading the compiled extension after adding each candidate
   directory one at a time, to see which one (if any) actually fixes it.
3. If none do, probes for the Microsoft VC++ Redistributable by name, since
   that is a plausible culprit entirely outside what our own code manages.

Uses ``importlib.util.find_spec`` rather than ``import magtense`` so the
diagnostic itself does not trip over the same failure it is debugging.
"""

import ctypes
import importlib.util
import os
import sys
from pathlib import Path


def try_load(pyd: Path, label: str) -> bool:
    try:
        ctypes.WinDLL(str(pyd))
        print(" ->", label, ": LoadLibrary SUCCEEDED")
        return True
    except OSError as e:
        print(" ->", label, ": LoadLibrary failed:", e)
        return False


def main() -> None:
    print("python executable:", sys.executable)
    print("sys.prefix:", sys.prefix)

    spec = importlib.util.find_spec("magtense")
    if spec is None or spec.origin is None:
        print("magtense package not importable even as a spec - wheel install itself may be broken")
        return
    pkg_dir = Path(spec.origin).parent
    print("magtense package dir:", pkg_dir)

    candidates = {
        "sys.prefix/Library/bin (conda-style)": Path(sys.prefix) / "Library" / "bin",
        "installed-layout Library/bin": (pkg_dir / ".." / ".." / ".." / "Library" / "bin").resolve(),
        "nvidia/cu13/bin/x86_64": (pkg_dir / ".." / "nvidia" / "cu13" / "bin" / "x86_64").resolve(),
    }
    for label, d in candidates.items():
        exists = d.is_dir()
        print()
        print("---", label, ":", d, "(exists=", exists, ") ---")
        if exists:
            dlls = sorted(p.name for p in d.glob("*.dll"))
            print(" ", len(dlls), "dll(s):", ", ".join(dlls) if dlls else "(empty)")

    ext_dir = pkg_dir / "lib"
    pyds = sorted(ext_dir.glob("*.pyd"))
    print()
    print("Extension dir", ext_dir, ":", [p.name for p in pyds])
    if not pyds:
        print("No .pyd found - wheel install is incomplete")
        return
    pyd = pyds[0]

    print()
    print("Attempting to load", pyd.name, "with no extra directories added:")
    if try_load(pyd, "baseline"):
        return

    added = []
    for label, d in candidates.items():
        if d.is_dir():
            os.add_dll_directory(str(d))
            os.environ["PATH"] = str(d) + os.pathsep + os.environ.get("PATH", "")
            added.append(label)
            print()
            print("Added:", label)
            if try_load(pyd, "after adding " + ", ".join(added)):
                return

    print()
    print("Still failing after adding every known candidate directory.")
    print("Probing for the Microsoft VC++ Redistributable by name (not managed by")
    print("our own DLL-path code at all, so a hit or miss here is informative")
    print("either way):")
    for dll_name in ["vcruntime140.dll", "vcruntime140_1.dll", "msvcp140.dll"]:
        try:
            ctypes.WinDLL(dll_name)
            print(" ", dll_name, "-> resolves OK")
        except OSError as e:
            print(" ", dll_name, "-> FAILED:", e)


if __name__ == "__main__":
    main()
