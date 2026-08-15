import numpy as np, time
from skimage.feature import structure_tensor, structure_tensor_eigenvalues
a = np.fromfunction(lambda k,j,i: ((k+1)*7 + (j+1)*13 + (i+1)*3) % 251 / 251.0, (32,48,64))
def lam_min(vol, sigma=1.0):
    S = structure_tensor(vol, sigma=sigma, order='rc')
    e = structure_tensor_eigenvalues(S)
    return e[-1]                      # eigenvalues come sorted descending
def t_(f, n=3):
    for _ in range(1): f()
    s=time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000
print("skimage structure_tensor + eigenvalues (3-D) %.2f ms" % t_(lambda: lam_min(a)))
print("skimage structure_tensor only                %.2f ms" % t_(lambda: structure_tensor(a, sigma=1.0, order='rc')))
