"""Module to deal with analog traces."""

from dataclasses import dataclass
from datetime import timedelta
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.signal import (
    find_peaks,
    peak_widths,
)


@dataclass
class AnalogTracePeak:
    """A class containing trace peak information.

    Attributes
    ----------
    time
        A `timedelta` for the retention time.
    intensity
        Maximum intensity of the peak.
    integral
        Relative integral of the peak.
    relative_height
        Height from which the integral is calculated.
    lhs, rhs
        Indices

    """

    time: timedelta
    intensity: float
    integral: float
    relative_height: float
    lhs: timedelta
    rhs: timedelta


@dataclass(init=False)
class AnalogTrace:
    """A class containing an analog trace."""

    times: NDArray[np.float64]
    intensities: NDArray[np.float64]
    peaks: list[AnalogTracePeak]

    def __init__(
        self,
        times: NDArray[np.float64],
        intensities: NDArray[np.float64],
        peaks: Optional[list[AnalogTracePeak]] = None,
    ):
        self.times = times
        self.intensities = intensities / intensities.max()
        self.peaks = peaks if peaks is not None else []

    def get_peaks(
        self,
        solvent_front=0.4,
        run_end=3.0,
        **kwargs,
    ) -> list[AnalogTracePeak]:
        """Get peaks in an AnalogTrace.

        Uses `scipy.signal.find_peaks()` to identify peaks in the MS scan.

        Parameters
        ----------
        **kwargs
            Optional arguments passed to `scipy.signal.find_peaks()`.

        Returns
        -------
            List of Peaks found in the trace.

        """

        # % of maximum peak height at which the width is established
        rel_height = kwargs["rel_height"] if "rel_height" in kwargs else 0.95

        # Solvent front and end of run index
        sfidx = self.get_scan_index(solvent_front)
        eidx = self.get_scan_index(run_end)

        run_times = self.times[sfidx:eidx]
        run_int = self.intensities[sfidx:eidx]

        peak_idx, _ = find_peaks(x=run_int, **kwargs)
        _, peak_height, peak_lhs, peak_rhs = peak_widths(
            x=run_int,
            peaks=peak_idx,
            rel_height=rel_height,
        )
        peaks = []
        integrals: list[float] = []
        total_area = 0.0

        for height, lhs, rhs in zip(peak_height, peak_lhs, peak_rhs):
            lhs_idx = round(lhs)
            rhs_idx = round(rhs)
            area: float = np.trapz(
                run_int[lhs_idx:rhs_idx] - height,
                dx=1,
            )
            integrals.append(area)

            total_area += area

        for x, idx in enumerate(peak_idx):
            peaks.append(
                AnalogTracePeak(
                    time=timedelta(minutes=run_times[idx]),
                    intensity=run_int[idx],
                    integral=integrals[x] / total_area,
                    lhs=timedelta(minutes=run_times[round(peak_lhs[x])]),
                    rhs=timedelta(minutes=run_times[round(peak_rhs[x])]),
                    relative_height=peak_height[x],
                )
            )

        self.peaks.extend(peaks)
        return peaks

    def get_scan_index(
        self,
        time: float | timedelta,
    ) -> int:
        """Get index of a mass scan at a certain time.

        Parameters
        ----------
        time
            Time (in minutes or as timedelta) to find the index for.

        Returns
        -------
            Array index of the MS scan.

        """
        if isinstance(time, timedelta):
            time = time / timedelta(minutes=1)

        return (np.absolute(self.times - time)).argmin()
