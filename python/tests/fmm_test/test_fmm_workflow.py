import os
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


@pytest.mark.skipif(
    not _is_fmm3d_wheel(),
    reason="FMM3D quick test requires the cu12-fmm MagTense wheel.",
)
def test_fmm_quick_workflow() -> None:
    test_dir = Path(__file__).parent
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")

    subprocess.run(
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
        check=True,
    )
