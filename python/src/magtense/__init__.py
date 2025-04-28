import importlib.metadata
from pathlib import Path

if (Path(__file__).parent.parent.parent / "pyproject.toml").exists():
    # Set dynamically in .github/workflows/python-package-conda.yml
    # Fallback if not set
    __version__ = "1.0.0"
else:
    __version__ = importlib.metadata.version("magtense")
