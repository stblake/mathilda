"""Branch scan for d11 (six masses): which (theta, T) solve r(T) = (4,1), and
which has the smallest T.  Vectorised fixed-step RK4 (h = 2e-4) over 40001
launch angles in [-pi, pi], t <= 2.3; every local minimum (in theta) of the
closest approach to (4,1) below 0.02 is refined by fsolve on (theta, T) with
solve_ivp DOP853 rtol=atol=1e-12.  Result (T, theta), sorted:
  1.9160969183 -0.274881254   <- shortest
  1.9176286707  0.5512858114
  1.9999300548 -0.196971957
  2.0236768226  0.3285693099
  2.065343358   0.7322047148
  2.098161426   0.3548510553
  2.1526863822 -0.4612274201  ... (and later arrivals)
~30 min single core.  Agrees with D. Smith's scan (AppendixD/Smith13.pdf)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
M = np.array([[2,-1],[2,0],[2,1],[3,-1],[3,0],[3,1]],float)
def acc(P):  # P: (n,2)
    d = M[None,:,:]-P[:,None,:]
    r3 = (d**2).sum(-1)**1.5
    return (d/r3[...,None]).sum(1)
th = np.linspace(-np.pi, np.pi, 40001)
P = np.zeros((th.size,2)); V = np.stack([np.cos(th),np.sin(th)],1)
h=2e-4; best=np.full(th.size,9.); tb=np.zeros(th.size)
t=0
for k in range(int(2.3/h)):
    # RK4 on (P,V)
    a1=acc(P); k1p=V; k1v=a1
    a2=acc(P+h/2*k1p); k2p=V+h/2*k1v; k2v=a2
    a3=acc(P+h/2*k2p); k3p=V+h/2*k2v; k3v=a3
    a4=acc(P+h*k3p); k4p=V+h*k3v; k4v=a4
    P=P+h/6*(k1p+2*k2p+2*k3p+k4p); V=V+h/6*(k1v+2*k2v+2*k3v+k4v)
    t+=h
    dd=np.hypot(P[:,0]-4,P[:,1]-1)
    m=dd<best; best[m]=dd[m]; tb[m]=t
# local minima of best below 0.02
idx=[i for i in range(1,th.size-1) if best[i]<best[i-1] and best[i]<=best[i+1] and best[i]<0.02]
def F(z):
    s=solve_ivp(lambda t,u:[u[2],u[3],*acc(np.array([u[:2]]))[0]],[0,z[1]],[0,0,np.cos(z[0]),np.sin(z[0])],rtol=1e-12,atol=1e-12,method='DOP853')
    return [s.y[0,-1]-4,s.y[1,-1]-1]
res=[]
for i in idx:
    z,info,ier,msg=fsolve(F,[th[i],tb[i]],full_output=True)
    if ier==1 and max(abs(np.array(F(z))))<1e-9: res.append((z[1],z[0]))
for T,a in sorted(set((round(T,10),round(a,10)) for T,a in res)): print(T,a)
print("nmin", len(idx), "nres", len(res))
for i in idx: print("cand", th[i], tb[i], best[i])
