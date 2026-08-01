#!/usr/bin/env python3
"""Experiment 10 -- Second sweep: eight real-application pipelines

Eight pipelines out of real applications rather than eight builtins.  See
``README.md``.

    python3 hpc_sweep_applications.py

NumPy is the memory-bandwidth reference for these rows: it makes the same
number of passes over the same data, so a spread against it is overhead rather
than algorithm.

WARNING: this file allocates several hundred megabytes.
"""

import numpy as np

import time


def bench(label, fn, reps=3):
    """One untimed warm-up, then the MINIMUM of `reps` timed runs.

    The minimum, not the mean: we are measuring the cost of the work, and every
    source of noise on a loaded machine can only add.
    """
    fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    print("%-52s%s ms" % (label, round(1000.0 * min(ts), 3)))


def bench1(label, fn):
    """A single timed run, no warm-up, for kernels that take seconds."""
    t0 = time.perf_counter()
    fn()
    print("%-52s%s ms  (1 run)" % (label, round(1000.0 * (time.perf_counter() - t0), 3)))


def check(label, value):
    print("%-52scheck = %s" % (label, value))


def nest(f, x, m):
    """Nest[f, x, m] -- apply f to x, m times."""
    for _ in range(m):
        x = f(x)
    return x

# ---- bsmc -- Black-Scholes Monte Carlo, 10^7 paths
# Computational finance: the whole option-pricing workload -- draw, Box-
# Muller, lognormal, payoff, reduce, over 10^7 paths.
npth=10**7
s0=100.0
strike=100.0
rr=0.05
sig=0.2
tt=1.0
def bs(m):
    u1 = np.random.rand(m); u2 = np.random.rand(m)
    z = np.sqrt(-2.0*np.log(u1))*np.cos(6.283185307179586*u2)
    st = s0*np.exp((rr-0.5*sig**2)*tt + sig*np.sqrt(tt)*z)
    pay = np.where(st-strike >= 0.0, st-strike, 0.0)
    return float(math.exp(-rr*tt)*pay.sum()/m)

# ---- tseries -- Return series: EMA + rolling vol + drawdown, 10^6
# Sequential scans: an exponential moving average, a rolling volatility and a
# maximum drawdown.
from scipy.signal import lfilter
nts=10**6
wts=250
ret=np.random.rand(nts)*0.02-0.01
def tsr(r, w):
    px = np.exp(np.cumsum(r))
    ema = lfilter([0.02],[1.0,-0.98], px[1:], zi=[0.98*px[0]])[0]
    run = np.maximum.accumulate(px)
    cs = np.cumsum(r*r)
    vol = np.sqrt((cs[w:]-cs[:-w])/w)
    return (float(ema[-1]) + float((1.0-px/run).max())
            + float(vol.max()))

# ---- logreg -- Logistic regression, 200000x32, 100 GD steps
# Machine learning: matrix-vector both ways with loop-carried state.
nlr=200000
dlr=32
Xlr=np.random.rand(nlr,dlr)*2.0-1.0
wtr=np.random.rand(dlr)*2.0-1.0
ylr=(Xlr @ wtr >= 0.0).astype(float)
def lrfit(xm, yv, m):
    nn, dd = xm.shape
    w = np.zeros(dd)
    for _ in range(m):
        pr = 1.0/(1.0+np.exp(-(xm @ w)))
        w = w - 4.0*((xm.T @ (pr-yv))/nn)
    return float(w @ w)
Xchk=np.array([[float((i*k) % 5)/5.0-0.4 for k in range(1,5)]
               for i in range(1,65)])
ychk=np.array([float(i % 2) for i in range(1,65)])

# ---- kmeans -- k-means, 100000 points, 8 dims, k = 16, 20 its
# Masked reductions, an argmin and a scatter-as-dot.
nkm=100000
dkm=8
kkm=16
Pkm=np.random.rand(nkm,dkm).T.copy()
def km(ptsT, kc, m):
    nn = ptsT.shape[1]
    cs = [ptsT[:,j].copy() for j in range(kc)]
    def dist(c):
        return np.array([((ptsT-q[:,None])**2).sum(axis=0) for q in c])
    for _ in range(m):
        d = dist(cs); mnv = d.min(axis=0)
        asg = np.zeros(nn); new = []
        for j in range(kc):
            mj = np.where(mnv-d[j] >= 0.0, 1.0, 0.0)*(1.0-asg)
            asg = asg + mj; cnt = mj.sum()
            new.append((ptsT @ mj)/cnt if cnt > 0.0 else cs[j])
        cs = new
    return float(dist(cs).min(axis=0).sum())
Pchk=np.array([[float((i*i+3*i*k) % 101)/101.0+0.001*k
                for k in range(1,5)] for i in range(1,201)]).T.copy()

# ---- nbody -- N-body all-pairs gravity, 1024 bodies, 10 steps
# Astrophysics: Outer over all pairs, with megabyte intermediates.
nbd=1024
st0=[np.random.rand(nbd) for _ in range(6)]
def grav(s, dt=0.001, eps=0.01):
    px,py,pz,vx,vy,vz = s
    dx = px[:,None]-px[None,:]
    dy = py[:,None]-py[None,:]
    dz = pz[:,None]-pz[None,:]
    r2 = dx*dx+dy*dy+dz*dz+eps
    iv = r2**-1.5
    ax = -(dx*iv).sum(axis=1); ay = -(dy*iv).sum(axis=1)
    az = -(dz*iv).sum(axis=1)
    return [px+dt*(vx+0.5*dt*ax), py+dt*(vy+0.5*dt*ay),
            pz+dt*(vz+0.5*dt*az), vx+dt*ax, vy+dt*ay, vz+dt*az]
schk=[np.array([float((7*i+3*k) % 11)/11.0 for k in range(1,9)])
      for i in range(1,7)]

# ---- heat3d -- 3D heat equation, 128^3, 50 steps
# CFD: a 3-D 7-point stencil over a 16 MB working set.
n3=128
u3=np.random.rand(n3,n3,n3)
def hstep(u):
    return u + 0.1*(np.roll(u,-1,0)+np.roll(u,1,0)
                    +np.roll(u,-1,1)+np.roll(u,1,1)
                    +np.roll(u,-1,2)+np.roll(u,1,2)-6.0*u)
uchk=np.array([[[float((i+2*j+3*k) % 5) for k in range(1,9)]
                for j in range(1,9)] for i in range(1,9)])

# ---- psd -- Welch PSD, 2^22 samples, 1024 blocks of 4096
# Signal processing: many medium FFTs, windowed.
nsg=2**22
blk=4096
nblk=1024
sgn=np.random.rand(nsg)*2.0-1.0
wnd=0.5*(1.0-np.cos(2.0*math.pi*np.arange(blk)/(blk-1)))
def psd(s, w, nb, bl):
    acc = np.zeros(bl)
    rt = math.sqrt(bl)
    for i in range(nb):
        acc = acc + np.abs(np.fft.fft(s[i*bl:i*bl+bl]*w)/rt)**2
    return acc
schk=np.sin(0.1*np.arange(1,4097))
wchk=0.5*(1.0-np.cos(2.0*math.pi*np.arange(0,256)/255.0))

# ---- imgpipe -- Gaussian blur + Sobel edges, 1024^2
# Vision: a separable Gaussian blur and a Sobel edge filter, rank-2
# correlation throughout.
imr=np.random.rand(1024,1024)
gk=np.array([[1.,4.,6.,4.,1.]])/16.0
gkT=gk.T.copy()
sxk=np.array([[-1.,0.,1.],[-2.,0.,2.],[-1.,0.,1.]])
syk=np.array([[-1.,-2.,-1.],[0.,0.,0.],[1.,2.,1.]])
def edge(im):
    b = spsig.correlate(spsig.correlate(im, gk, mode='valid'),
                      gkT, mode='valid')
    gx = spsig.correlate(b, sxk, mode='valid')
    gy = spsig.correlate(b, syk, mode='valid')
    return np.sqrt(gx*gx+gy*gy)
ichk=np.array([[float((i*j) % 17) for j in range(1,65)]
               for i in range(1,65)])


def main():
    print("Experiment 10 -- Second sweep: eight real-application pipelines")
    print("")
    bench1("Black-Scholes Monte Carlo, 10^7 paths", lambda: bs(npth))
    check("Black-Scholes Monte Carlo, 10^7 paths", 100.0*0.5*(1.0+math.erf(0.35/math.sqrt(2.0))) - 100.0*math.exp(-0.05)*0.5*(1.0+math.erf(0.15/math.sqrt(2.0))))
    bench1("Return series: EMA + rolling vol + drawdown, 10^6", lambda: tsr(ret,wts))
    check("Return series: EMA + rolling vol + drawdown, 10^6", tsr(np.array([float((3*i) % 7)/700.0-0.004 for i in range(1,2001)]), 50))
    bench1("Logistic regression, 200000x32, 100 GD steps", lambda: lrfit(Xlr,ylr,100))
    check("Logistic regression, 200000x32, 100 GD steps", lrfit(Xchk,ychk,50))
    bench1("k-means, 100000 points, 8 dims, k = 16, 20 its", lambda: km(Pkm,kkm,20))
    check("k-means, 100000 points, 8 dims, k = 16, 20 its", km(Pchk,4,10))
    bench1("N-body all-pairs gravity, 1024 bodies, 10 steps", lambda: nest(grav,st0,10))
    check("N-body all-pairs gravity, 1024 bodies, 10 steps", float(sum(a.sum() for a in nest(grav,schk,3))))
    bench1("3D heat equation, 128^3, 50 steps", lambda: nest(hstep,u3,50))
    check("3D heat equation, 128^3, 50 steps", float(nest(hstep,uchk,5).sum()))
    bench1("Welch PSD, 2^22 samples, 1024 blocks of 4096", lambda: psd(sgn,wnd,nblk,blk))
    check("Welch PSD, 2^22 samples, 1024 blocks of 4096", float(psd(schk,wchk,16,256).sum()))
    bench1("Gaussian blur + Sobel edges, 1024^2", lambda: edge(imr))
    check("Gaussian blur + Sobel edges, 1024^2", float(edge(ichk).sum()))


if __name__ == "__main__":
    main()
