"""Module to deal with ion traces."""

from dataclasses import dataclass
from datetime import timedelta
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.signal import find_peaks

from lcms_parser.helpers.helpers import (
    IonTraceMode,
    normalised,
)


@dataclass
class TICTracePeak:
    """A class containing trace peak information."""

    mode: IonTraceMode
    time: timedelta
    intensity: float


@dataclass(init=False)
class TICTrace:
    """A class containing a Total Ion Chromatogram trace."""

    mode: IonTraceMode
    times: NDArray[np.float64]
    intensities: NDArray[np.float64]
    peaks: list[TICTracePeak]

    def __init__(
        self,
        mode: IonTraceMode,
        times: NDArray[np.float64],
        intensities: NDArray[np.float64],
        peaks: Optional[list[TICTracePeak]] = None,
    ):
        self.mode = mode
        self.times = times
        self.intensities = intensities
        self.peaks = peaks if peaks is not None else []

    def get_peaks(
        self,
        **kwargs,
    ) -> list[TICTracePeak]:
        """Get peaks in a TICTrace.

        Uses `scipy.signal.find_peaks()` to identify peaks in the MS scan.

        Parameters
        ----------
        **kwargs
            Optional arguments passed to `scipy.signal.find_peaks()`.

        Returns
        -------
            List of Peaks found in the trace.

        """
        peak_idx, _ = find_peaks(x=normalised(self.intensities), **kwargs)
        peaks = []

        for idx in peak_idx:
            peaks.append(
                TICTracePeak(
                    mode=self.mode,
                    time=timedelta(minutes=self.times[idx]),
                    intensity=self.intensities[idx],
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
