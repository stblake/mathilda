"""64-image-resize -- skimage as the baseline for ImageResize.

downscale_local_mean is the right comparison for the area path: it is exactly a block mean at an
integer factor, which is what Mathilda's area averaging reduces to when the fractional coverage
terms vanish. skimage's `resize` with anti_aliasing=True would prefilter with a Gaussian first and
so measure a different algorithm.

For the enlargement, resize(order=1, anti_aliasing=False) is plain bilinear, matching Mathilda's
default for enlargement (there are no frequencies to remove when enlarging).
"""
import numpy as np
from skimage.transform import resize, downscale_local_mean


def ramp(n):
    y, x = np.mgrid[1:n + 1, 1:n + 1]
    return ((x * 7 + y * 13) % 251) / 251.0


img = ramp(512)

bench("IR area 512->256", lambda: downscale_local_mean(img, (2, 2)))
bench("IR area 512->128", lambda: downscale_local_mean(img, (4, 4)))
bench("IR nearest 512->256", lambda: resize(img, (256, 256), order=0, anti_aliasing=False))
bench("IR bilinear 512->1024", lambda: resize(img, (1024, 1024), order=1, anti_aliasing=False))

small = downscale_local_mean(img, (2, 2))
check("IR area preserves mean", round(1e9 * (small.mean() - img.mean())))

chk = np.array([[0., 1., 0., 1.], [1., 0., 1., 0.],
                [0., 1., 0., 1.], [1., 0., 1., 0.]])
check("IR checkerboard area total", round(1e6 * downscale_local_mean(chk, (2, 2)).sum()))

q = np.array([[0., .25, .5, .75], [0., .25, .5, .75],
              [1., 1., 1., 1.], [1., 1., 1., 1.]])
check("IR block mean exact", round(1e9 * downscale_local_mean(q, (2, 2)).sum()))
