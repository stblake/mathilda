"""Reference for SIAM100 Appendix D problem 15 (experiment 190): blow-up time of
u_t = Lap u + e^u on [0,3]^2, u = 0 on the boundary and at t = 0.

Chebyshev collocation on [-3/2, 3/2]^2, folded by the even symmetry, with the
centre-clustering map x = 1.5 (s - beta sin(pi s)/pi), Radau (rtol 1e-13).
Two formulations:
  v-form (main):  v = e^-u, v_t = Lap v - |grad v|^2/v - 1, v = 1 on the
     boundary; blow-up <=> v(centre) = 0, which v reaches LINEARLY and smoothly
     (v ~ tau + |x|^2/(8|log tau|)).  Stop at v(centre) = eps0, add v/|v_t|.
     beta=0.95, eps0=1e-7: N=48 2.1197288160, 56 2.1197288115, 64 2.1197288099
     beta=0.95, eps0=1e-8: N=48 2.1197288159, 56 2.1197288113, 64 2.1197288098
     beta=0.98, eps0=1e-7: N=48 2.1197288184
  u-form (check): stop at u(centre) = M, add e^-M (u ~ -log(t_inf - t)):
     beta=0.95, M=24, N=48: 2.1197287976; beta=0.9, M=20, N=48: 2.1197286166
Consensus: 2.1197288 (8 digits).  WL 4th-order FD on the v-form (d15_ref.wl)
corroborates 2.11973 (6 digits; n=61 gives 2.1197213, h^4 convergence).
Usage: d15_ref.py v|u beta eps0|M N...
"""
import numpy as np, time, sys
from scipy.integrate import solve_ivp
def cheb(N):
    x=np.cos(np.pi*np.arange(N+1)/N); c=np.r_[2,np.ones(N-1),2]*(-1)**np.arange(N+1)
    X=np.tile(x,(N+1,1)).T; dX=X-X.T; D=np.outer(c,1/c)/(dX+np.eye(N+1)); D-=np.diag(D.sum(1)); return D,x
def fold(A,N,h):
    F=np.zeros((h,h))
    for i in range(1,h+1):
        for j in range(1,N):
            F[i-1,(j if j<=h else N-j)-1]+=A[i,j]
    return F
def ops(N,beta):
    D,s=cheb(N); gp=1.5*(1-beta*np.cos(np.pi*s)); Dx=D/gp[:,None]; D2=Dx@Dx
    h=N//2; F1=fold(Dx,N,h); F2=fold(D2,N,h); I=np.eye(h)
    return h, np.kron(F2,I)+np.kron(I,F2), np.kron(F1,I), np.kron(I,F1)
def run_v(N,beta,eps0):
    h,K,Gx,Gy=ops(N,beta)
    def f(t,w):
        v=1+w; gx=Gx@w; gy=Gy@w
        return K@w-(gx*gx+gy*gy)/v-1
    def jac(t,w):
        v=1+w; gx=Gx@w; gy=Gy@w; g2=gx*gx+gy*gy
        return K-(2*gx[:,None]*Gx+2*gy[:,None]*Gy)/v[:,None]+np.diag(g2/v**2)
    c=h*h-1; ev=lambda t,w:1+w[c]-eps0; ev.terminal=True
    sol=solve_ivp(f,[0,50],np.zeros(h*h),method='Radau',jac=jac,rtol=1e-13,atol=1e-14,events=ev)
    assert sol.status==1
    w=sol.y[:,-1]; return sol.t[-1]+(1+w[c])/(-f(0,w)[c])
def run_u(N,beta,M):
    h,K,_,_=ops(N,beta); c=h*h-1
    ev=lambda t,u:u[c]-M; ev.terminal=True
    sol=solve_ivp(lambda t,u:K@u+np.exp(u),[0,50],np.zeros(h*h),method='Radau',
                  jac=lambda t,u:K+np.diag(np.exp(u)),rtol=1e-13,atol=1e-13,events=ev)
    assert sol.status==1
    return sol.t[-1]+np.exp(-M)
if __name__=="__main__":
    form=sys.argv[1]; beta=float(sys.argv[2]); p=float(sys.argv[3])
    for N in map(int,sys.argv[4:]):
        t0=time.time(); T=(run_v if form=='v' else run_u)(N,beta,p)
        print(form,N,"%.15f"%T,"%.1fs"%(time.time()-t0),flush=True)
