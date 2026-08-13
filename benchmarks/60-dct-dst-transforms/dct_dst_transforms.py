#!/usr/bin/env python3
"""Experiment 60 -- Discrete cosine/sine & 2-D transforms (scipy.fft/numpy column).

Same kernels as ``dct_dst_transforms.m``, same order and labels.

Baseline is compiled (scipy.fft pocketfft), a compiled-vs-compiled O(n log n)
comparison.  NORMALIZATION is reconciled to Mathematica's convention in the
check, verified against Mathilda element-by-element:

  DCT type 1  ==  scipy.fft.dct(type=1, norm='ortho', orthogonalize=False)
  DCT type 2  ==  scipy.fft.dct(type=2, norm='ortho') with y[1:] /= sqrt(2)
  DCT type 4  ==  scipy.fft.dct(type=4, norm='ortho')
  DST type 1  ==  scipy.fft.dst(type=1, norm='ortho')
  DST type 2  ==  scipy.fft.dst(type=2, norm='ortho') with y[:-1] /= sqrt(2)
                  (the last/Nyquist term is unscaled, mirroring DCT-2's index 0)
  Fourier 2D  ==  numpy.fft.fft2 / sqrt(n_total)

Timing runs on a 10^6 vector / 512x512 matrix; checks use {1..8}/{{1,2},{3,4}}.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
import scipy.fft as fft

from harness import bench, check, require, seed

require(["scipy.fft:dct", "scipy.fft:dst", "numpy.fft:fft2"])

seed()
n = 1000000
v = np.random.random(n)
c8 = np.arange(1, 9.0)
mbig = np.random.random((512, 512))
m2 = np.array([[1.0, 2.0], [3.0, 4.0]])
S2 = math.sqrt(2.0)


def r6(x):
    return int(math.floor(float(x) * 1e6 + 0.5))


bench("FourierDCT type 1 of 10^6", lambda: fft.dct(v, type=1, norm="ortho", orthogonalize=False))
check("FourierDCT type 1 of 10^6",
      r6(fft.dct(c8, type=1, norm="ortho", orthogonalize=False).sum()))

bench("FourierDCT type 2 of 10^6", lambda: fft.dct(v, type=2, norm="ortho"))
_d2 = fft.dct(c8, type=2, norm="ortho"); _d2[1:] /= S2
check("FourierDCT type 2 of 10^6", r6(_d2.sum()))

bench("FourierDCT type 4 of 10^6", lambda: fft.dct(v, type=4, norm="ortho"))
check("FourierDCT type 4 of 10^6", r6(fft.dct(c8, type=4, norm="ortho").sum()))

bench("FourierDST type 1 of 10^6", lambda: fft.dst(v, type=1, norm="ortho"))
check("FourierDST type 1 of 10^6", r6(fft.dst(c8, type=1, norm="ortho").sum()))

bench("FourierDST type 2 of 10^6", lambda: fft.dst(v, type=2, norm="ortho"))
_s2 = fft.dst(c8, type=2, norm="ortho"); _s2[:-1] /= S2
check("FourierDST type 2 of 10^6", r6(_s2.sum()))

bench("Fourier 2D 512x512", lambda: np.fft.fft2(mbig))
check("Fourier 2D 512x512", r6(np.abs((np.fft.fft2(m2) / 2.0).flatten()).sum()))
