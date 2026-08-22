import numpy as np, time
from scipy import ndimage
z, y, x = 64, 96, 128
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (z,y,x))
def mk(n):
    return np.fromfunction(lambda m,i,j: ((m+1)*5 + (i+1)*3 + (j+1)*2) % 11 - 5.0, (n,n,n))
def t_(f, n=3):
    for _ in range(1): f()
    s = time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
for n in (3,5,7,9):
    k = mk(n)
    print("k=%d^3 scipy ndimage.convolve  %.2f ms" % (n, t_(lambda: ndimage.convolve(a, k, mode='nearest'))))
