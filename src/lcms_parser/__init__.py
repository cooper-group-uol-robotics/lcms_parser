"""Toolkit for parsing UPLC-MS instrument files."""

from lcms_parser.experimental.hits import HitIdentifier
from lcms_parser.experimental.planner import ExpectedResults
from lcms_parser.fileio.waters import WatersRawFile
from lcms_parser.helpers.helpers import (
    IonTraceMode,
    normalised,
)
from lcms_parser.msdata.peak import (
    MassPeak,
    MassSpectrumExperimentalHit,
    MassSpectrumResult,
)
from lcms_parser.msdata.spectrum import MassSpectrum
from lcms_parser.traces.analog import (
    AnalogTrace,
    AnalogTracePeak,
)
from lcms_parser.traces.ion import (
    TICTrace,
    TICTracePeak,
)

__version__ = "0.3.5"
__all__ = (
    "ExpectedResults",
    "HitIdentifier",
    "WatersRawFile",
    "IonTraceMode",
    "normalised",
    "MassPeak",
    "MassSpectrumResult",
    "MassSpectrumExperimentalHit",
    "MassSpectrum",
    "TICTrace",
    "TICTracePeak",
    "AnalogTrace",
    "AnalogTracePeak",
)
