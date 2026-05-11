import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


def _is_fmm3d_wheel() -> bool:
    return (
        os.environ.get("MAGTENSE_USE_FMM3D", "").strip().lower()
        in {"1", "true", "yes", "on"}
        and os.environ.get("MAGTENSE_WHEEL_VARIANT") == "cu12-fmm"
    )


def _has_nvidia_gpu() -> bool:
    # This test validates the cu12-fmm wheel by comparing FMM output against a
    # CUDA baseline, so a CUDA runtime package alone is not enough.
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return False

    return (
        subprocess.run(
            [nvidia_smi, "-L"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        ).returncode
        == 0
    )


def _tail(text: str, max_lines: int = 200) -> str:
    lines = text.splitlines()
    if len(lines) <= max_lines:
        return text
    return "\n".join(
        [f"... showing last {max_lines} of {len(lines)} lines ...", *lines[-max_lines:]]
    )


def _read_tail(path: Path, max_lines: int = 80) -> str:
    try:
        return _tail(path.read_text(errors="replace"), max_lines=max_lines)
    except OSError as exc:
        return f"Could not read {path}: {exc}\n"


def _print_failure_diagnostics(
    test_dir: Path,
    completed: subprocess.CompletedProcess[str],
) -> None:
    print(f"\nFMM workflow failed with exit code {completed.returncode}")
    print("\n--- subprocess stdout ---")
    print(_tail(completed.stdout), end="")
    print("\n--- subprocess stderr ---", file=sys.stderr)
    print(_tail(completed.stderr), end="", file=sys.stderr)

    diagnostic_paths = [
        test_dir / "error.txt",
        *sorted((test_dir / "timer_logs").glob("*_42_timer.log")),
    ]
    for path in diagnostic_paths:
        if path.exists():
            print(f"\n--- {path.relative_to(test_dir)} ---", file=sys.stderr)
            print(_read_tail(path), end="", file=sys.stderr)


@pytest.mark.skipif(
    not _is_fmm3d_wheel(),
    reason="FMM3D quick test requires the cu12-fmm MagTense wheel.",
)
@pytest.mark.skipif(
    os.environ.get("GITHUB_ACTIONS") == "true" and not _has_nvidia_gpu(),
    reason="FMM3D quick test compares against CUDA and requires an NVIDIA GPU.",
)
def test_fmm_quick_workflow() -> None:
    # Smoke-test the installed cu12-fmm wheel through the grain workflow:
    # generate exchange/demag setup, run the CUDA reference, run quick FMM
    # variants, and validate the FMM fields against the CUDA result.
    test_dir = Path(__file__).parent
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("PYTHONFAULTHANDLER", "1")
    env.setdefault("CUDA_LAUNCH_BLOCKING", "1")

    completed = subprocess.run(
        [
            sys.executable,
            "run_fmm_test.py",
            "--quick-test",
            "--seed",
            "42",
            "--Hext-dir-type",
            "2",
        ],
        cwd=test_dir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    if completed.returncode != 0:
        _print_failure_diagnostics(test_dir, completed)
        completed.check_returncode()

    print(completed.stdout, end="")
    print(completed.stderr, end="", file=sys.stderr)
