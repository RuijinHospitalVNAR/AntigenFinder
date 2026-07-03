"""预测器模块"""

from .base_predictor import BasePredictor
from .bepipred_wrapper import BepipredWrapper
from .discotope_wrapper import DiscotopeWrapper
from .graphbepi_wrapper import GraphbepiWrapper
from .epigraph_wrapper import EpigraphWrapper

__all__ = [
    'BasePredictor',
    'BepipredWrapper',
    'DiscotopeWrapper',
    'GraphbepiWrapper',
    'EpigraphWrapper',
]
