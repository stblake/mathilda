import numpy as np, time
from skimage.filters import threshold_otsu
from scipy import ndimage
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (32,48,64))
def t_(f, n=5):
    for _ in range(2): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
print("skimage threshold_otsu (3-D)      %.2f ms" % t_(lambda: a > threshold_otsu(a)))
for r in (2,4,8):
    k=2*r+1
    print("scipy uniform_filter local r=%d    %.2f ms" % (r, t_(lambda: a > ndimage.uniform_filter(a, size=k, mode='nearest'))))
def sauvola(vol, r, c2):
    k=2*r+1
    mu = ndimage.uniform_filter(vol, size=k, mode='nearest')
    msq = ndimage.uniform_filter(vol*vol, size=k, mode='nearest')
    sd = np.sqrt(np.maximum(msq - mu*mu, 0))
    return vol > mu + c2*sd
print("scipy sauvola-style r=4           %.2f ms" % t_(lambda: sauvola(a, 4, -0.2)))
