import numpy as np, time
from skimage.filters import threshold_local
a = np.fromfunction(lambda i,j: ((i+1)*7 + (j+1)*13) % 251 / 251.0, (512,512))
def t_(f, n=5):
    for _ in range(2): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
for r in (2,8,16,32):
    bs = 2*r+1
    print("skimage threshold_local block=%d (mean)  %.2f ms" % (bs, t_(lambda: a > threshold_local(a, bs, method='mean'))))
