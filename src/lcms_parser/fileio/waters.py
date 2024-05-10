"""Main class to interact with MassLynx RAW files.

Requires Waters Connect Python SDK.

"""

from datetime import timedelta
from os import PathLike

import numpy as np
from masslynx import (
    MassLynxRawAnalogReader,
    MassLynxRawChromatogramReader,
    MassLynxRawInfoReader,
    MassLynxRawScanReader,
)

from lcms_parser.experimental.hits import HitIdentifier
from lcms_parser.helpers.helpers import IonTraceMode
from lcms_parser.msdata.spectrum import MassSpectrum
from lcms_parser.traces.analog import AnalogTrace
from lcms_parser.traces.ion import TICTrace


class WatersRawFile(HitIdentifier):
    """Raw MassLynx file and associated properties."""

    def __init__(
        self,
        path: PathLike,
    ):
        """Initialise MSFile.

        Parameters
        ----------
        path
            A path to the .RAW file.

        """
        self._info_reader = MassLynxRawInfoReader(str(path))
        self._analog_reader = MassLynxRawAnalogReader(str(path))
        self._chromatogram_reader = MassLynxRawChromatogramReader(str(path))
        self._scan_reader = MassLynxRawScanReader(str(path))
        self._trace_function_lookup = self._get_number_functions()
        self.traces: dict[IonTraceMode, TICTrace] = {}
        self.analog: AnalogTrace = self.get_analog_trace()

        # self.TICpeaks = []
        # self.analogpeaks = []
        # self.analogpeakwidths = []
        # self.analogpeakfraction = []
        # self.rtmax = 1
        # self.MSScanmax = 1

    def _get_number_functions(self):
        """Get a lookup dictionary for function numbers of trace types.

        Waters .RAW files contain trace information under different numbers
        rather than names (e.g., ES+). This function creates a lookup table
        for implementation references. Many MassLynx SDK functions (e.g.,
        ReadTIC()) take the function number rather than trace name as the
        arguments.

        Returns
        -------
            Lookup table: value(trace_type) = function_number.

        """
        trace_function_lookup = {}
        num_functions = self._info_reader.GetNumberofFunctions()
        for func in range(num_functions):
            # Functions (electrospray TIC or DAD) are numbered in .RAW.
            # Need to identify which number corresponds to the trace.
            trace = str(
                self._info_reader.GetIonModeString(
                    self._info_reader.GetIonMode(func)
                )
            )
            trace_function_lookup[trace] = func
        return trace_function_lookup

    def __get_channel_count(self):
        """Implementation of MassLynx SDK GetChannelCount C++ function."""
        from ctypes import (
            POINTER,
            c_int,
            c_void_p,
        )

        size = c_int(0)
        ch_count = self._info_reader.massLynxDll.getChannelCount
        ch_count.argtypes = [c_void_p, POINTER(c_int)]
        self._info_reader.CheckReturnCode(
            ch_count(self._info_reader._getReader(), size)
        )
        return size.value

    def __get_number_scans(self, function):
        return self._info_reader.GetScansInFunction(function)

    def get_analog_trace(self) -> AnalogTrace:
        times, intensities = self._analog_reader.ReadChannel(0)
        data = AnalogTrace(
            times=np.array(times), intensities=np.array(intensities)
        )

        return data

    def get_trace(
        self,
        mode: IonTraceMode,
    ) -> TICTrace:
        """Get Trace from the MS data file.

        Parameters
        ----------
        mode
            Requested trace, for allowed modes see `IonTraceMode`.

        Returns
        -------
            Trace containing requested data.

        """
        if mode in self.traces:
            return self.traces[mode]

        else:
            function = self._trace_function_lookup[mode]

            times, intensities = self._chromatogram_reader.ReadTIC(function)

            data = TICTrace(
                mode=mode,
                times=np.array(times),
                intensities=np.array(intensities),
            )
            self.traces[mode] = data

            return data

    def get_mass_spectrum(
        self,
        time: float | timedelta,
        mode: IonTraceMode,
        average: int = 0,
    ) -> MassSpectrum:
        """Get MS scan at a specific time.

        Parameters
        ----------
        time
            Time (in minutes or as `timedelta`) to extract the scan for.
        mode
            Mass trace to examine.
        average, optional
            Number of scans around the time to average the MS scans for.

        Returns
        -------
            MassSpectrum at the given time.

        """

        function = self._trace_function_lookup[mode]
        trace = self.get_trace(mode=mode)
        idx = trace.get_scan_index(time=time)

        intensities = []
        for scan in range(idx - average, idx + average + 1):
            masses, intensity = self._scan_reader.ReadScan(function, scan)
            intensities.append(intensity)
        mean_intensity = np.array(intensities).mean(axis=0)

        return MassSpectrum(
            masses=masses, intensities=mean_intensity, mode=mode
        )

    def get_peak_mass_spectrum(
        self,
        mode: IonTraceMode,
        peak_idx: int = 0,
        average: int = 0,
    ) -> MassSpectrum:
        """Get MS scan for a specific trace peak.

        This is an API function for extracting MS scans at a peak in the TIC.
        By default, one scan at the apex of the peak is extracted by one can
        also specify how many scans (left and right from the apex) to average.

        Parameters
        ----------
        mode
            Mass trace to examine.Mass trace to examine.
        peak_idx, optional
            Index of the TracePeak to examine, by default it is the first peak.
        average, optional
            Number of scans around the time to average the MS scans for.

        Returns
        -------
            MassSpectrum at the given TIC peak.

        """
        trace = self.get_trace(mode)
        try:
            if len(trace.peaks) == 0:
                peak = trace.get_peaks()[peak_idx]
            else:
                peak = trace.peaks[peak_idx]

        except IndexError:
            raise IndexError("Peak does not exist.")

        return self.get_mass_spectrum(
            time=peak.time, mode=peak.mode, average=average
        )
