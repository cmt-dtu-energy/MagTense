import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.skipif(
    os.environ.get("MAGTENSE_USE_FMM3D") != "1",
    reason="FMM3D quick test requires an FMM-enabled MagTense build.",
)
def test_fmm_quick_workflow() -> None:
    test_dir = Path(__file__).parent

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
        check=True,
    )
