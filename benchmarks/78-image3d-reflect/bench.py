import numpy as np, time
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (32,48,64))
def t_(f, n=10):
    for _ in range(3): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
for ax in (1, 2, 0):
    print("numpy flip(axis=%d) + copy  %.3f ms" % (ax, t_(lambda: np.flip(a, axis=ax).copy())))
