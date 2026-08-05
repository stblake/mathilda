#!/usr/bin/env python3
"""Experiment 17 -- FFT, convolution, correlation (scipy.fft / numpy column).

Same six kernels as ``fft_and_signal.m``, same order and sizes.

Power-of-two AND non-power-of-two sizes are both here, because a mixed-radix
size is exactly where a weak FFT falls back to something quadratic.  numpy's
pocketfft handles all three well, so the spread across these rows is a property
of the CAS columns.

NORMALISATION.  Wolfram's ``Fourier`` divides by sqrt(n); numpy's ``fft`` does
not.  Only the CHECK values are affected -- the timing is the same transform --
so each check divides by sqrt(8) for the 8-point reference sequence.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

import numpy as np

from harness import bench, check, require, seed

require(["numpy.fft:fft", "numpy.fft:ifft", "numpy:convolve",
         "numpy:correlate", "scipy.signal:fftconvolve"])

seed()
small = np.array([1., 2., 3., 4., 5., 6., 7., 8.])
S = math.sqrt(8.0)


def r6(x):
    return int(math.floor(float(x) * 1e6 + 0.5))


d1 = np.random.random(262144)
bench("Fourier 2^18 (262144)", lambda: np.fft.fft(d1))
check("Fourier 2^18 (262144)", r6(np.fft.fft(small)[0].real / S))

d2 = np.random.random(1200000)
bench("Fourier 1200000 (mixed radix)", lambda: np.fft.fft(d2))
check("Fourier 1200000 (mixed radix)", r6(np.fft.fft(small)[-1].real / S))

d3 = np.random.random(262143)
bench("Fourier 262143 (awkward size)", lambda: np.fft.fft(d3))
# Wolfram uses e^{+2 pi i}, numpy e^{-2 pi i}: the imaginary part is conjugated,
# so compare |Im| and the row measures the transform rather than the convention.
check("Fourier 262143 (awkward size)", r6(abs(np.fft.fft(small)[1].imag) / S))

# Wolfram's InverseFourier also carries a 1/sqrt(n); the round trip is identity
# in both systems, which is what the check pins.
bench("InverseFourier 2^18", lambda: np.fft.ifft(d1))
check("InverseFourier 2^18", r6(np.fft.ifft(np.fft.fft(small))[0].real))

a = np.random.random(100000)
b = np.random.random(2048)
# Wolfram's ListConvolve[b, a] returns the VALID (no-padding) region.
bench("ListConvolve 100000 x 2048", lambda: np.convolve(a, b, mode="valid"))
# ListConvolve[{1,2},{3,4,5}] == {10, 13}: the kernel is applied reversed.
check("ListConvolve 100000 x 2048",
      r6(np.convolve(np.array([3., 4., 5.]), np.array([1., 2.]),
                     mode="valid")[0]))

bench("ListCorrelate 100000 x 2048", lambda: np.correlate(a, b, mode="valid"))
check("ListCorrelate 100000 x 2048",
      r6(np.correlate(np.array([3., 4., 5.]), np.array([1., 2.]),
                      mode="valid")[0]))
