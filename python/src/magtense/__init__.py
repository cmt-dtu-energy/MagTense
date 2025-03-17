import importlib.metadata
import pathlib
import tomllib

toml_location = pathlib.Path(__file__).parent.parent.parent
if (toml_location / "pyproject.toml").exists():
    with open(toml_location / "pyproject.toml", "rb") as f:
        __version__ = tomllib.load(f)["project"]["version"]
else:
    __version__ = importlib.metadata.version("magtense")
