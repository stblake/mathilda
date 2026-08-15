import numpy as np, time
from scipy import signal
a = np.fromfunction(lambda i,j: ((i+1)*7 + (j+1)*13) % 251 / 251.0, (512,512))
t32 = a[99:131, 199:231].copy()
t8  = a[99:107, 199:207].copy()

def ncc(img, t):
    """The same identities: cross term by correlation, window stats by integral images."""
    kh, kw = t.shape; m = kh*kw
    tb = t.mean(); tn = np.sqrt(((t-tb)**2).sum())
    ci, cj = kh//2, kw//2
    p = np.pad(img, ((ci, kh-1-ci), (cj, kw-1-cj)), mode='edge')
    s1 = np.pad(p.cumsum(0).cumsum(1), ((1,0),(1,0)))
    s2 = np.pad((p*p).cumsum(0).cumsum(1), ((1,0),(1,0)))
    H, W = img.shape
    def box(s):
        return (s[kh:kh+H, kw:kw+W] - s[0:H, kw:kw+W] - s[kh:kh+H, 0:W] + s[0:H, 0:W])
    sum1, sum2 = box(s1), box(s2)
    cross = signal.correlate2d(img, t, mode='same', boundary='symm')
    var = np.maximum(sum2 - sum1*sum1/m, 0.0)
    sn = np.sqrt(var)
    out = np.zeros_like(img)
    ok = (sn > 0) & (tn > 0)
    out[ok] = (cross[ok] - sum1[ok]*tb) / (sn[ok]*tn)
    return out

def t_(f, n=5):
    for _ in range(2): f()
    s = time.perf_counter()
    for _ in range(n): f()
    return (time.perf_counter()-s)/n*1000

print("scipy NCC 32x32     %.2f ms" % t_(lambda: ncc(a, t32)))
print("scipy NCC 8x8       %.2f ms" % t_(lambda: ncc(a, t8)))
print("scipy correlate2d 32 %.2f ms" % t_(lambda: signal.correlate2d(a, t32, mode='same', boundary='symm')))
