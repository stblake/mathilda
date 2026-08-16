import numpy as np, time
from scipy import ndimage
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (32,48,64))
def t_(f, n=3):
    for _ in range(1): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
for r in (1,2,4):
    k=2*r+1
    print("scipy uniform_filter size=%d^3  %.2f ms" % (k, t_(lambda: ndimage.uniform_filter(a, size=k, mode='nearest'), 5)))
for r in (1,2):
    k=2*r+1
    print("scipy median_filter  size=%d^3  %.2f ms" % (k, t_(lambda: ndimage.median_filter(a, size=k, mode='nearest'))))
