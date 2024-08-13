"""Module to deal with identifying hits."""

from datetime import timedelta
from typing import (
    Optional,
    Protocol,
)

from lcms_parser.experimental.planner import ExpectedResults
from lcms_parser.helpers.helpers import IonTraceMode
from lcms_parser.msdata.peak import MassSpectrumExperimentalHit
from lcms_parser.msdata.spectrum import MassSpectrum
from lcms_parser.traces.ion import TICTrace


class Experimental(Protocol):
    """Protocol class for mixins that operate on MS data."""

    def get_trace(self, mode: IonTraceMode) -> TICTrace:
        """Get trace."""
        ...

    def get_mass_spectrum(
        self,
        time: float | timedelta,
        mode: IonTraceMode,
        average: int = 0,
    ) -> MassSpectrum:
        """Get mass spectrum."""
        ...


class HitIdentifier(Experimental):
    """Mixin class to identify hits in MS .RAW file."""

    def identify_hits(
        self,
        expected_results: ExpectedResults,
        mode: IonTraceMode = "ES+",
        atol: float = 0.4,
        direct_injection: bool = True,
        time: Optional[timedelta] = None,
        tic_peak_params: Optional[dict[str, float]] = None,
        ms_peak_params: Optional[dict[str, float]] = None,
    ) -> list[MassSpectrumExperimentalHit]:
        """Identify hits in the .RAW file.

        Parameters
        ----------
        expected_results
            Expected results (m/z values).
        mode, optional
            Ionisation mode, by default "ES+".
        atol, optional
            Tolerance for m/z value to identify a hit, by default 0.4.
        direct_injection, optional
            If direct injection, the first TIC peak will be used.
        time, optional
            Time at which the MS spectrum should be analysed. If direct
            injection, then the first TIC peak will be used regardless. If no
            time is provided, the first TIC peak will be used.
        tic_peak_params, optional
            Scipy's find_peaks parameters for TIC peak identification.
        ms_peak_params, optional
            Scipy's find_peaks parameters for MS peak identification.

        Returns
        -------
            A list of experimental hits identified in the MS.

        """

        if time is None:
            direct_injection = True

        else:
            ms_time = time

        if direct_injection:
            if tic_peak_params is None:
                tic_peak_params = {"height": 0.2, "distance": 50}

            trace = self.get_trace(mode=mode)
            tic_peaks = trace.get_peaks(**tic_peak_params)
            ms_time = tic_peaks[0].time

        if ms_peak_params is None:
            ms_peak_params = {"height": 0.5, "distance": 30}

        ms = self.get_mass_spectrum(
            mode=mode,
            time=ms_time,
            average=0,
        )

        ms_peaks = ms.get_peaks(**ms_peak_params)

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
                        time=ms_time.total_seconds() / 60,
                    )
                )

        return hits
