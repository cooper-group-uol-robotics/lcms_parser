"""Module to deal with identified or expected MS peaks."""

from dataclasses import dataclass

from lcms_parser.helpers.helpers import IonTraceMode


@dataclass(kw_only=True)
class MassPeak:
    """A base class containing MS peak information."""

    mz_value: float
    mode: IonTraceMode


@dataclass(kw_only=True)
class MassSpectrumResult(MassPeak):
    """A class containing mass spectrometry results."""

    formula: str
    charge: int


@dataclass(kw_only=True)
class MassSpectrumExperimentalHit(MassSpectrumResult):
    """A class for keeping track of a mass spectrometry hit."""

    mz_expected: float
    time: float

    def get_relative_error(self) -> float:
        """Get relative error.

        A helper method to get a relative error - useful when reporting
        high resolution mass spectrometry results.

        Returns
        -------
            Relative error of the mz_value result.

        """
        return self.get_absolute_error() / self.mz_expected

    def get_absolute_error(self) -> float:
        """Get absolute error.

        A helper method to get a absolute error - useful when reporting
        high resolution mass spectrometry results.

        Returns
        -------
            Absolute error of the observed result.

        """
        return abs(self.mz_value - self.mz_expected)
