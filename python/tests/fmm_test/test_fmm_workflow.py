import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.skipif(
    os.environ.get("MAGTENSE_USE_FMM3D") != "1"
    or os.environ.get("MAGTENSE_WHEEL_VARIANT") != "cu12-fmm",
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
