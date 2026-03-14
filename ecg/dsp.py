# -*- coding: utf-8 -*-
"""Pure-Python drop-in replacement for scipy.signal.butter and scipy.signal.lfilter.

Removes the scipy dependency from processor.py.

Limitations:
  - butter(): only btype='lowpass' is supported
  - lfilter(): general IIR implementation, same interface as scipy
"""

from __future__ import annotations

import math
import cmath
from collections.abc import Sequence

import numpy as np


def _poly_from_roots(roots):
    """Return coefficients of the monic polynomial with the given roots.

    Coefficients are ordered from highest to lowest degree.
    """
    p = [1.0 + 0j]
    for r in roots:
        new_p = [0j] * (len(p) + 1)
        for i, c in enumerate(p):
            new_p[i] += c
            new_p[i + 1] -= c * r
        p = new_p
    return p


def butter(order: int, Wn: float, btype: str = 'lowpass') -> tuple[list[float], list[float]]:
    """Design a digital Butterworth IIR filter.

    Drop-in replacement for scipy.signal.butter with btype='lowpass'.

    Method: bilinear transform with cutoff frequency pre-warping.
    The analog Butterworth prototype poles of order N are uniformly
    distributed on the left half of the unit circle:

        p_k = exp(j * pi * (2k + N + 1) / (2N)),  k = 0, ..., N-1

    Pre-warping applies wc = 2 * tan(pi * Wn / 2) to preserve the
    cutoff frequency exactly after the bilinear transform s = 2*(z-1)/(z+1).

    Digital poles are obtained as z_k = (1 + wc/2 * p_k) / (1 - wc/2 * p_k).
    For a lowpass filter of order N all N digital zeros are at z = -1.

    :param order: filter order (integer >= 1)
    :param Wn: normalised cutoff frequency in (0, 1), where 1 = Nyquist
    :param btype: filter type; only 'lowpass' is supported
    :returns: (b, a) — numerator and denominator coefficient lists
    """
    if btype != 'lowpass':
        raise ValueError("Only btype='lowpass' is supported.")

    # Pre-warp the cutoff frequency for the bilinear transform
    wc = 2.0 * math.tan(math.pi * Wn / 2.0)

    # Butterworth analog prototype poles, uniformly distributed on the
    # left half of the unit circle (angles in (pi/2, 3pi/2)).
    # Formula: theta_k = pi * (2k + N + 1) / (2N), k = 0, ..., N-1
    analog_poles = [
        cmath.exp(1j * math.pi * (2 * k + order + 1) / (2 * order))
        for k in range(order)
    ]

    # Bilinear transform: z = (1 + wc/2 * s) / (1 - wc/2 * s)
    half_wc = wc / 2.0
    digital_poles = [
        (1.0 + half_wc * p) / (1.0 - half_wc * p)
        for p in analog_poles
    ]

    # Digital zeros: all at z = -1 for a lowpass filter
    digital_zeros = [-1.0 + 0j] * order

    # Build numerator and denominator polynomials from their roots
    b = _poly_from_roots(digital_zeros)
    a = _poly_from_roots(digital_poles)

    # Normalise for unity DC gain (z = 1)
    gain = sum(a).real / sum(b).real
    b = [c * gain for c in b]

    return [c.real for c in b], [c.real for c in a]


def lfilter(b: Sequence[float], a: Sequence[float], x: np.ndarray) -> np.ndarray:
    """Apply a causal IIR filter to signal x.

    Drop-in replacement for scipy.signal.lfilter.

    Implements direct form I with zero initial conditions:

        y[n] = b[0]*x[n] + b[1]*x[n-1] + ... - a[1]*y[n-1] - a[2]*y[n-2] - ...

    The FIR (feedforward) part is computed in one shot with np.convolve (C speed);
    the IIR (feedback) part remains sequential due to its recursive nature.

    :param b: numerator (feedforward) coefficients
    :param a: denominator (feedback) coefficients; a[0] must not be zero
    :param x: input signal (array-like)
    :returns: numpy array of the filtered signal, same length as x
    """
    a0 = a[0]
    b = np.asarray(b) / a0
    a = np.asarray(a) / a0
    na = len(a)

    # FIR part: vectorised convolution in C, truncated to the first len(x) samples
    fir = np.convolve(b, x)[:len(x)]

    # IIR part: sequential — each sample depends on previous outputs
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = fir[i]
        for j in range(1, min(na, i + 1)):
            y[i] -= a[j] * y[i - j]

    return y
