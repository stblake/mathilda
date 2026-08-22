"""63-image-convolve -- scipy.ndimage.convolve as the baseline for ImageConvolve.

scipy.ndimage.convolve is the right comparison rather than scipy.signal.convolve2d: ndimage
reflects the kernel (true convolution, as ImageConvolve does) and offers mode='nearest', which
replicates the edge pixel exactly like Mathilda's "Fixed" padding. convolve2d's boundary='symm'
is a mirror rather than a replication, so its edge pixels would differ and every check comparing
a border value would disagree for a reason that has nothing to do with either implementation.

The kernel is built here rather than taken from a library so that both sides use provably the
same numbers: a normalised Gaussian with sigma = r/2, which is Mathematica's GaussianMatrix
convention.
"""
import numpy as np
from scipy.ndimage import convolve


def gaussian_matrix(r):
    s = 1.0 if r == 0 else r / 2.0
    i, j = np.mgrid[-r:r + 1, -r:r + 1]
    k = np.exp(-(i ** 2 + j ** 2) / (2.0 * s * s))
    return k / k.sum()


def ramp(n):
    y, x = np.mgrid[1:n + 1, 1:n + 1]
    return ((x * 7 + y * 13) % 251) / 251.0


img = ramp(512)
small = ramp(64)
k1, k2, k4, k0 = (gaussian_matrix(r) for r in (1, 2, 4, 0))

bench("IC gaussian r=1 512x512", lambda: convolve(img, k1, mode='nearest'))
bench("IC gaussian r=2 512x512", lambda: convolve(img, k2, mode='nearest'))
bench("IC gaussian r=4 512x512", lambda: convolve(img, k4, mode='nearest'))
bench("IC marshalling floor r=0 512x512", lambda: convolve(img, k0, mode='nearest'))
bench("GaussianFilter r=2 512x512", lambda: convolve(img, k2, mode='nearest'))

g = convolve(small, k2, mode='nearest')
check("IC total r=2 64x64", round(1e6 * g.sum()))
check("IC corner r=2 64x64", round(1e9 * g[0, 0]))
check("IC centre r=2 64x64", round(1e9 * g[31, 31]))

flat = convolve(np.full((8, 8), 0.25), k2, mode='nearest')
check("IC constant preserved", round(1e9 * np.abs(flat - 0.25).max()))
