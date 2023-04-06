"""Exposure type data class."""

from dataclasses import dataclass

@dataclass
class ExpType:
    """Different kinds of exposure types."""
    bias: str = 'Bias'
    flat: str = 'Flat'
    science: str = 'Object'
