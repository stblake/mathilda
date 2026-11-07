# D13 landscape: equilibria of the ODE (saddle at ||u||=0.2152) and bisection for the threshold; trajectories near it.
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
A=np.array([[-1,10,0],[0,-1,10],[0,0,-1.]]); B=np.array([[0,-0.7,-0.03],[0,0.7,-0.1],[0.03,0.1,0]])
det=lambda n: np.linalg.det(A+n*B)
ns=np.linspace(0,50,5001); d=[det(n) for n in ns]
for i in range(len(ns)-1):
    if d[i]*d[i+1]<0:
        n=brentq(det,ns[i],ns[i+1]); w,v=np.linalg.eig(A+n*B); k=np.argmin(abs(w)); u=np.real(v[:,k]); u=u/np.linalg.norm(u)*n
        print("equilibrium norm",n,u, "jac eigs", np.linalg.eigvals(A+n*B+np.outer(B@u,u)/n))
def f(t,u):
    n=np.linalg.norm(u); return A@u+n*(B@u)
def M(z0,T=50):
    s=solve_ivp(f,[0,T],[0,0,z0],method='DOP853',rtol=1e-13,atol=1e-30)
    return np.linalg.norm(s.y[:,-1]) if s.status==0 else np.inf
lo,hi=6.3e-3,1e-2
for _ in range(60):
    m=(lo+hi)/2
    if M(m)<1: lo=m
    else: hi=m
print(lo,hi,M(lo),M(hi))
for z in [lo*(1-1e-3),lo*(1-1e-6), lo*(1-1e-9), lo, hi]:
    s=solve_ivp(f,[0,50],[0,0,z],method='DOP853',rtol=1e-13,atol=1e-30,dense_output=True)
    print(z, [ "%.3g"%np.linalg.norm(s.sol(t)) for t in [1,2,5,10,20,30,40,45,50]] if s.status==0 else s.t[-1])
