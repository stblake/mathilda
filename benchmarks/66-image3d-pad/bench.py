import numpy as np, time
z, y, x = 64, 96, 128
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (z,y,x))
def t(f, n=10):
    for _ in range(2): f()
    s = time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
print("np.pad constant  %.3f ms" % t(lambda: np.pad(a, 8, mode='constant')))
print("np.pad reflect   %.3f ms" % t(lambda: np.pad(a, 8, mode='reflect')))
b = np.pad(a, 8, mode='constant')
print("crop slice+copy  %.3f ms" % t(lambda: b[8:72, 8:104, 8:136].copy()))
