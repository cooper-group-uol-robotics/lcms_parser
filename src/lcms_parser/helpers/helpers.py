"""Module to define common functions and types."""

from typing import Literal

import numpy as np
from numpy.typing import (
    ArrayLike,
    NDArray,
)

IonTraceMode = Literal["ES+", "ES-"]


def normalised(array: ArrayLike) -> NDArray:
    """Get normalised array (min-max scaling).

    Parameters
    ----------
    array
        Array to normalise.

    Returns
    -------
        Normalised array (values between 0 and 1).

    """
    array = np.array(array)
    return (array - np.min(array)) / (np.max(array) - np.min(array))
