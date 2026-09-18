"""Gristmill base module.

Public names are going to be imported here.
"""

from .generate import (
    BasePrinter,
    NaiveCodePrinter,
    CPrinter,
    FortranPrinter,
    BlasFortranPrinter,
    EinsumPrinter,
    OMEinsumPrinter,
    mangle_base,
    BlockSpec,
    split_input_blocks,
)
from .optimize import optimize, verify_eval_seq, ContrStrat, RepeatedTermsStrat
from .utils import get_flop_cost

__version__ = "0.9.0"

__all__ = [
    "ContrStrat",
    "RepeatedTermsStrat",
    "optimize",
    "verify_eval_seq",
    "get_flop_cost",
    "BasePrinter",
    "mangle_base",
    "NaiveCodePrinter",
    "CPrinter",
    "FortranPrinter",
    "BlasFortranPrinter",
    "EinsumPrinter",
    "OMEinsumPrinter",
    "BlockSpec",
    "split_input_blocks",
]
