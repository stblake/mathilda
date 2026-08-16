import numpy as np, time
from scipy import ndimage
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (32,48,64))
def t_(f, n=5):
    for _ in range(2): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
def gradmag(v):
    return np.sqrt(sum(ndimage.sobel(v, axis=ax, mode='nearest')**2 for ax in range(3)))
print("scipy 3x sobel + magnitude   %.2f ms" % t_(lambda: gradmag(a)))
print("scipy sobel one axis         %.2f ms" % t_(lambda: ndimage.sobel(a, axis=2, mode='nearest')))
