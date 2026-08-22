import numpy as np, time
from skimage.feature import corner_harris, corner_shi_tomasi, corner_peaks
a = np.fromfunction(lambda i,j: ((i+1)*7 + (j+1)*13) % 251 / 251.0, (512,512))
def t_(f, n=5):
    for _ in range(2): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
print("skimage corner_shi_tomasi  %.2f ms" % t_(lambda: corner_shi_tomasi(a, sigma=1.0)))
print("skimage corner_harris      %.2f ms" % t_(lambda: corner_harris(a, sigma=1.0)))
r = corner_shi_tomasi(a, sigma=1.0)
print("skimage corner_peaks       %.2f ms" % t_(lambda: corner_peaks(r, min_distance=1), 3))
