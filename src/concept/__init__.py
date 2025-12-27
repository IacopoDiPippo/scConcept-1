from importlib.metadata import PackageNotFoundError, version

from .model import ContrastiveModel
from .api import scConcept

__all__ = ["ContrastiveModel", "scConcept"]

try:
    __version__ = version("scConcept")
except PackageNotFoundError:
    # fallback when running from source (editable / not installed)
    __version__ = "0.0.0"
