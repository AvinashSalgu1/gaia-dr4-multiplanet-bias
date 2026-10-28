"""GEDI (Gaia Exoplanet Detection and Inference): how unmodelled companions bias Gaia DR4 planet detection, orbits, occurrence rates and yields."""

from .epochs import Epochs
from .orbits import Companion
from .simulate import Star, simulate

__all__ = ["Epochs", "Companion", "Star", "simulate"]
__version__ = "0.3.0"
