import numpy as np, time
a = np.array([[(i*7+j*13) % 251 / 251.0 for j in range(1,513)] for i in range(1,513)])
def t(f, n=20):
    for _ in range(3): f()
    s = time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
print("np.pad constant      %.3f ms" % t(lambda: np.pad(a, 16, mode='constant')))
print("np.pad reflect       %.3f ms" % t(lambda: np.pad(a, 16, mode='reflect')))
b = np.pad(a, 16, mode='constant')
print("crop (slice+copy)    %.3f ms" % t(lambda: b[16:528, 16:528].copy()))
