"""Module to deal with mass spectra.

A mass spectrum can consist of a single or multiple scans.

"""

from typing import (
    TYPE_CHECKING,
    List,
    Optional,
)

import numpy as np
from numpy.typing import NDArray
from scipy.signal import find_peaks

from lcms_parser.helpers.helpers import (
    IonTraceMode,
    normalised,
)
from lcms_parser.msdata.peak import (
    MassPeak,
    MassSpectrumExperimentalHit,
)

if TYPE_CHECKING:
    from lcms_parser.experimental.planner import ExpectedResults


class MassSpectrum:
    """A class containing mass spectrum data."""

    def __init__(
        self,
        masses: NDArray,
        intensities: NDArray,
        mode: IonTraceMode,
    ):
        self.masses = masses
        self.intensities = intensities
        self.mode = mode
        self.experimental_hits: List[MassSpectrumExperimentalHit] = []

    def __repr__(self):
        return (
            f"MassSpectrum(masses={repr(self.masses)}, intensities="
            f"{repr(self.intensities)}, mode={self.mode})"
        )

    def __eq__(self, other):
        if isinstance(other, MassSpectrum):
            return np.allclose(self.masses, other.masses) and np.allclose(
                self.intensities, other.intensities
            )
        else:
            raise TypeError(f"Cannot compare MassSpectrum with {type(other)}.")

    def get_peaks(
        self,
        **kwargs,
    ) -> list[MassPeak]:
        """Get peaks in a MassSpectrum.

        Uses `scipy.signal.find_peaks()` to identify peaks in the MS scan.

        Parameters
        ----------
        **kwargs
            Optional arguments passed to `scipy.signal.find_peaks()`.

        Returns
        -------
            List of Peaks found in the spectrum.

        """
        peak_idx, _ = find_peaks(x=normalised(self.intensities), **kwargs)
        peaks = []

        for idx in peak_idx:
            peaks.append(MassPeak(mz_value=self.masses[idx], mode=self.mode))

        return peaks

    def identify_peaks(
        self,
        expected_results: "ExpectedResults",
        atol: float = 0.4,
        peak_params: Optional[dict[str, float]] = None,
    ) -> list[MassSpectrumExperimentalHit]:
        """Identify MS hits in a MassSpectrum.

        This function used the MassSpectrum object directly (useful when
        performing separate analysis for each peak in the chromatogram). It
        will also add identified peaks to the `experimental_hits` list.

        Parameters
        ----------
        expected_results
            ExpectedResults to look for in the MS data.
        trace, optional
            Ionisation mode, by default "ES+"
        atol, optional
                Tolerance for identification, by default 0.4 Daltons.

        Returns
        -------
            A list of MSHits identified in the experimental data.

        """
        if peak_params is None:
            peak_params = {"height": 0.5, "distance": 30}

        ms_peaks = self.get_peaks(**peak_params)

        hits = []
        for peak in ms_peaks:
            for hit in expected_results.find(peak, atol=atol):
                hits.append(
                    MassSpectrumExperimentalHit(
                        formula=hit.formula,
                        charge=hit.charge,
                        mode=hit.mode,
                        mz_value=peak.mz_value,
                        mz_expected=hit.mz_value,
                        time=0.0,
                    )
                )

        self.experimental_hits.extend(hits)

        return hits
