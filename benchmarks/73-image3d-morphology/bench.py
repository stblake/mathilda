import numpy as np, time
from scipy import ndimage
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (32,48,64))
def t_(f, n=5):
    for _ in range(2): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
for r in (1,2,4,8):
    k = 2*r+1
    print("scipy grey_dilation size=%d^3  %.2f ms" % (k, t_(lambda: ndimage.grey_dilation(a, size=(k,k,k), mode='nearest'))))
print("scipy grey_opening  size=9^3   %.2f ms" % t_(lambda: ndimage.grey_opening(a, size=(9,9,9), mode='nearest')))
