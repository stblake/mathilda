import numpy as np, time
from scipy import ndimage
a = np.fromfunction(lambda k,j,i: np.where(((k+1)*5 + (j+1)*3 + (i+1)*7) % 23 == 0, 0.0, 1.0), (32,48,64))
sp = np.ones((32,48,64)); sp[15,23,31] = 0.0
def t_(f, n=5):
    for _ in range(2): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
print("scipy distance_transform_edt        %.2f ms" % t_(lambda: ndimage.distance_transform_edt(a)))
print("scipy distance_transform_edt sparse %.2f ms" % t_(lambda: ndimage.distance_transform_edt(sp)))
