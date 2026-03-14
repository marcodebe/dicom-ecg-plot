# -*- coding: utf-8 -*-
"""Pure-Python drop-in replacement for scipy.signal.butter and scipy.signal.lfilter.

Elimina la dipendenza da scipy per processor.py.

Limitazioni:
  - butter(): supporta solo btype='lowpass'
  - lfilter(): implementazione IIR generale, stessa interfaccia di scipy
"""

from __future__ import annotations

import math
import cmath
from collections.abc import Sequence

import numpy as np


def _poly_from_roots(roots):
    """Restituisce i coefficienti del polinomio monico con le radici date.

    I coefficienti sono ordinati dal grado più alto al più basso.
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
    """Progetta un filtro IIR Butterworth digitale.

    Drop-in replacement per scipy.signal.butter con btype='lowpass'.

    Metodo: bilinear transform con pre-warping della frequenza di taglio.
    I poli analogici del prototipo Butterworth di ordine N sono distribuiti
    uniformemente sulla semicirconferenza sinistra del cerchio unitario:

        p_k = exp(j * pi * (2k + N - 1) / (2N)),  k = 0, ..., N-1

    Il pre-warping applica wc = 2 * tan(pi * Wn / 2) per preservare
    esattamente la frequenza di taglio dopo la trasformazione bilineare
    s = 2*(z-1)/(z+1).

    I poli digitali si ottengono come z_k = (1 + wc/2 * p_k) / (1 - wc/2 * p_k).
    Per un filtro passa-basso di ordine N tutti gli N zeri digitali sono in z = -1.

    :param order: ordine del filtro (intero >= 1)
    :param Wn: frequenza di taglio normalizzata in (0, 1), dove 1 = Nyquist
    :param btype: tipo di filtro; solo 'lowpass' supportato
    :returns: (b, a) — liste di coefficienti numeratore e denominatore
    """
    if btype != 'lowpass':
        raise ValueError("Solo btype='lowpass' è supportato.")

    # Pre-warping della frequenza di taglio per la trasformazione bilineare
    wc = 2.0 * math.tan(math.pi * Wn / 2.0)

    # Poli del prototipo analogico Butterworth, uniformemente distribuiti
    # sul semipiano sinistro del cerchio unitario (angoli in (pi/2, 3pi/2)).
    # Formula: theta_k = pi * (2k + N + 1) / (2N), k = 0, ..., N-1
    analog_poles = [
        cmath.exp(1j * math.pi * (2 * k + order + 1) / (2 * order))
        for k in range(order)
    ]

    # Trasformazione bilineare: z = (1 + wc/2 * s) / (1 - wc/2 * s)
    half_wc = wc / 2.0
    digital_poles = [
        (1.0 + half_wc * p) / (1.0 - half_wc * p)
        for p in analog_poles
    ]

    # Zeri digitali: tutti in z = -1 per il filtro passa-basso
    digital_zeros = [-1.0 + 0j] * order

    # Polinomi numeratore e denominatore dalle rispettive radici
    b = _poly_from_roots(digital_zeros)
    a = _poly_from_roots(digital_poles)

    # Normalizzazione: guadagno unitario a DC (z = 1)
    gain = sum(a).real / sum(b).real
    b = [c * gain for c in b]

    return [c.real for c in b], [c.real for c in a]


def lfilter(b: Sequence[float], a: Sequence[float], x: np.ndarray) -> np.ndarray:
    """Applica un filtro IIR causale al segnale x.

    Drop-in replacement per scipy.signal.lfilter.

    Implementa la forma diretta I con condizioni iniziali nulle:

        y[n] = b[0]*x[n] + b[1]*x[n-1] + ... - a[1]*y[n-1] - a[2]*y[n-2] - ...

    La parte FIR (feedforward) è calcolata in blocco con np.convolve (C),
    la parte IIR (feedback) rimane sequenziale per natura ricorsiva.

    :param b: coefficienti del numeratore (feedforward)
    :param a: coefficienti del denominatore (feedback); a[0] != 0
    :param x: segnale in ingresso (array-like)
    :returns: numpy array del segnale filtrato, stessa lunghezza di x
    """
    a0 = a[0]
    b = np.asarray(b) / a0
    a = np.asarray(a) / a0
    na = len(a)

    # Parte FIR: convoluzione vettoriale in C, tronca ai primi len(x) campioni
    fir = np.convolve(b, x)[:len(x)]

    # Parte IIR: sequenziale (ogni campione dipende dai precedenti)
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = fir[i]
        for j in range(1, min(na, i + 1)):
            y[i] -= a[j] * y[i - j]

    return y
