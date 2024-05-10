"""Module to deal with planning what is expected."""


from collections import UserDict

import numpy as np

from lcms_parser.msdata.peak import (
    MassPeak,
    MassSpectrumResult,
)


class ExpectedResults(UserDict[str, list[MassSpectrumResult]]):
    """A class containing expected MS results.

    Simple dictionary class where the keys are predicted formulae and the
    values are lists of expected `MSResult`s.

    """

    def __setitem__(
        self,
        formula: str,
        results: list[MassSpectrumResult],
    ):
        self.data[formula] = results

    def __contains__(self, ms_peak):
        if not isinstance(ms_peak, MassPeak):
            raise TypeError("Expected results must be mass peaks!")
        values = []
        for ms_results in self.values():
            for ms_result in ms_results:
                values.append(ms_result.ms_peak)
        return ms_peak in values

    def find(
        self,
        value: float | MassPeak,
        atol: float = 1e-03,
    ) -> list[MassSpectrumResult]:
        """Find an m/z value or a MassPeak in the ExpectedResults.

        Looks through the m/z values of all the predicted formulae. If
        MassPeak is provided, then modes also need to agree.

        Parameters
        ----------
        value
            The m/z value (in Daltons) to find or a required MassPeak.
        atol, optional
            Tolerance for identification, by default 1e-03. For lower accuracy
            mass spectrometers, might want to decrease (e.g., to 0.5 Da).

        Returns
        -------
            List of MSResults with the required m/z value.

        """
        match value:
            case float():
                return self.find_value(value, atol)
            case MassPeak():
                return self.find_peak(value, atol)
            case _:
                raise TypeError(
                    f"Searching of type {type(value)} in ExpectedResults "
                    "not supported."
                )

    def find_value(
        self,
        mz_value: float,
        atol: float = 1e-03,
    ) -> list[MassSpectrumResult]:
        """Find an m/z value in the ExpectedResults.

        Looks through the m/z values of all the predicted formulae.

        Parameters
        ----------
        mz_value
            The m/z value (in Daltons) to find.
        atol, optional
            Tolerance for identification, by default 1e-03. For lower accuracy
            mass spectrometers, might want to decrease (e.g., to 0.5 Da).

        Returns
        -------
            List of MSResults with the required m/z value.

        """
        results = []
        for ms_results in self.values():
            for ms_result in ms_results:
                if np.isclose(ms_result.mz_value, mz_value, atol=atol):
                    results.append(ms_result)
        return results

    def find_peak(
        self,
        ms_peak: MassPeak,
        atol: float = 1e-03,
    ) -> list[MassSpectrumResult]:
        """Find a MassPeak in the ExpectedResults.

        Looks through the m/z values of all the predicted formulae and
        matches it agains the m/z value of the required MassPeak. Unlike
        ExpectedResults.find_value(), it also requires the modes to match,
        so an m/z=200 in ES+ will not match against m/z=200 in ES-.

        Parameters
        ----------
        mz_value
            MassPeak to be identified.
        atol, optional
            Tolerance for identification, by default 1e-03. For lower accuracy
            mass spectrometers, might want to decrease (e.g., to 0.5 Da).

        Returns
        -------
            List of MSResults with the required m/z value.

        """
        results = []
        for ms_results in self.values():
            for ms_result in ms_results:
                if (
                    np.isclose(
                        ms_result.mz_value,
                        ms_peak.mz_value,
                        atol=atol,
                    )
                    and ms_result.mode == ms_peak.mode
                ):
                    results.append(ms_result)
        return results

    @classmethod
    def from_dict(
        cls,
        expected,
    ):
        """Initialise ExpectedResults from a dictionary.

        Converts a dictionary of a form:
        {
            (formula : str) = {
                (charge : int | str) = (mz_value : float | str)
            }
        }

        Parameters
        ----------
        expected
            A dictionary of expected formulae.

        """
        results = ExpectedResults()
        for formula, masses in expected.items():
            ms_results = []
            for charge, mz_value in masses.items():
                ms_results.append(
                    MassSpectrumResult(
                        mz_value=float(mz_value),
                        mode="ES+",
                        charge=charge,
                        formula=formula,
                    )
                )
            results[formula] = ms_results

        return cls(results)
