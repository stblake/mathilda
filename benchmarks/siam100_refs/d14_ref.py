"""Reference for SIAM100 Appendix D problem 14 (experiment 189).

-Laplace u = exp(alpha |x|^2) on the regular pentagon inscribed in the unit
circle (vertices exp(i(pi/2 + 2 pi k/5))), u = 0 on four sides, du/dn = 0 on
the bottom side (W2, W3).  Find alpha with Int_{Neumann side} u ds = e^alpha.
(Rotation invariance: any choice of the Neumann side gives the same answer.)

Split u = u_p + h with the radial particular solution
    u_p(r) = -sum_{k>=1} alpha^(k-1) r^(2k) / (4 k k!)     (-Lap u_p = e^{alpha r^2})
and h = Re F(z) harmonic with h = -u_p on the Dirichlet sides and
dh/dn = -du_p/dn on the Neumann side, by boundary least squares in two
DIFFERENT bases sharing only the boundary data:
  Method A ("MPS"): corner (Fourier-Bessel type) expansions zeta^nu at every
    vertex -- nu = k*5/3 at the three Dirichlet-Dirichlet corners, (k-1/2)*5/3
    at the two Dirichlet-Neumann corners (the r^(5/6) singularity) -- plus
    z^m, m <= deg.  Residual 1.3e-12 at nk=22, deg=60.
  Method B ("lightning", Gopal-Trefethen): poles clustered exponentially at
    each vertex along the exterior bisector, plus z^m.
Int u ds along the Neumann side by graded Gauss-Legendre panels; brentq on
I(alpha) - e^alpha over [-3, 0] (I(alpha) e^{-alpha} is strictly decreasing,
so the root is unique).

Run (double precision):
  MPS 18/50 -2.4628347560427652, 22/60 -2.4628347560428643, 26/70 -2.462834756042903
  lightning 100 -2.462834756045318, 140 -2.462834756042388
=> consensus 12 digits -2.46283475604 (13th: 3 by MPS, 2 by lightning).
Method C (d14_ref.wl): WL FEM at MaxCellMeasure 2e-5 / 5e-6 with Richardson
(rate h^(5/3)) agrees to ~5-6 digits; it corroborates, not more.
"""
import numpy as np, math, sys
from scipy.optimize import brentq
W = np.exp(1j*(np.pi/2 + 2*np.pi*np.arange(5)/5))
sides = [(W[k], W[(k+1)%5]) for k in range(5)]   # side 2 = (W2,W3) bottom: Neumann
NEU = 2
def up(r2, a):
    # -sum_{k>=1} a^{k-1} r^{2k}/(4 k k!)
    s = np.zeros_like(r2); term_pow = r2.copy(); fact=1.0; ak=1.0
    for k in range(1,60):
        fact *= k
        s += ak*term_pow/(4*k*fact)
        term_pow = term_pow*r2; ak *= a
    return -s
def dup_dr_over_r(r2, a):   # u_p'(r)/r = -(e^{a r^2}-1)/(2 a r^2)
    x = a*r2
    return -np.where(np.abs(x)>1e-8, np.expm1(x)/np.where(x==0,1,x), 1+x/2)/2
def samples(M, sig=None):
    j = np.arange(1, M+1)
    sig = 29.0/(np.sqrt(M)-1)
    s = np.exp(-sig*(np.sqrt(M)-np.sqrt(j)))     # clustered near 0, in (0,1]
    t = np.concatenate([s/2, 1 - s/2])
    return np.sort(t)
def build(Np, deg, Ms, sig=4.0):
    poles=[]
    for k in range(5):
        d = np.exp(-min(sig, 26.0/(np.sqrt(Np)-1))*(np.sqrt(Np)-np.sqrt(np.arange(1,Np+1))))
        poles += list(W[k]*(1 + d))
    poles=np.array(poles)
    pts=[]; typ=[]; nrm=[]
    for k,(a,b) in enumerate(sides):
        t = samples(Ms)
        z = a + (b-a)*t
        n = -1j*(b-a)/abs(b-a)       # outward normal for CCW orientation
        pts.append(z); typ.append(np.full(len(z), k==NEU)); nrm.append(np.full(len(z), n))
    return poles, np.concatenate(pts), np.concatenate(typ), np.concatenate(nrm)
def basis(z, poles, deg, deriv=False):
    cols=[]
    for p in poles:
        cols.append(-1/(z-p)**2 if deriv else 1/(z-p))
    for m in range(deg+1):
        if deriv: cols.append(m*z**(m-1) if m>0 else 0*z)
        else: cols.append(z**m)
    return np.array(cols).T
def solve(a, Np, deg, Ms):
    poles, z, isN, n = build(Np, deg, Ms)
    B = basis(z, poles, deg); Bd = basis(z, poles, deg, True)
    # h = Re(B c), c complex -> real columns [Re B, -Im B]
    A_D = np.hstack([B.real, -B.imag])
    G = Bd*n[:,None]    # Re(F' n)
    A_N = np.hstack([G.real, -G.imag])
    A = np.where(isN[:,None], A_N, A_D)
    r2 = np.abs(z)**2
    rhs_D = -up(r2, a)
    xn = (z*np.conj(n)).real           # x . n
    rhs_N = -dup_dr_over_r(r2, a)*xn
    rhs = np.where(isN, rhs_N, rhs_D)
    sc = np.linalg.norm(A, axis=0); sc[sc==0]=1
    c, *_ = np.linalg.lstsq(A/sc, rhs, rcond=None); c/=sc
    res = np.max(np.abs(A@c - rhs))
    nc = len(poles)+deg+1
    cc = c[:nc] + 1j*c[nc:]
    # integral along Neumann side: graded Gauss-Legendre panels toward both ends
    a0,b0 = sides[NEU]; L = abs(b0-a0)
    xg, wg = np.polynomial.legendre.leggauss(30)
    br = [0.0] + list(0.5*0.15**np.arange(25,0,-1)) + [0.5]
    br = np.array(br); br = np.concatenate([br, 1-br[::-1][1:]])
    tot = 0.0
    for lo,hi in zip(br[:-1],br[1:]):
        t = lo + (hi-lo)*(xg+1)/2; w = wg*(hi-lo)/2
        zz = a0 + (b0-a0)*t
        h = (basis(zz, poles, deg) @ cc).real
        u = up(np.abs(zz)**2, a) + h
        tot += np.sum(w*u)*L
    return tot, res
TH = 3*np.pi/5
def cornerinfo(k):
    # interior wedge at W[k]: between side (k-1 -> k) reversed and side (k -> k+1)
    w = W[k]; d1 = W[(k+1)%5]-w; d0 = W[(k-1)%5]-w
    psi = np.angle(d1)          # phi=0 along side k (towards W[k+1]); CCW interior goes to d0 at phi=TH
    # side index of phi=0 edge is k, phi=TH edge is k-1
    return w, psi, k, (k-1)%5
def cbasis(z, nk, deg, deriv=False):
    cols=[]
    for k in range(5):
        w, psi, s0, s1 = cornerinfo(k)
        zz = (z - w)*np.exp(-1j*psi)
        rho = np.abs(zz)
        phi = np.angle(zz*np.exp(-1j*(TH/2 - np.pi)))  # cut along exterior bisector
        phi = np.mod(phi, 2*np.pi) + (TH/2 - np.pi)
        neu0 = (s0==NEU); neu1=(s1==NEU)
        for j in range(1, nk+1):
            if not neu0 and not neu1: nu = j*np.pi/TH; f = lambda ph: np.sin(nu*ph)   # D on both
            elif neu1: nu = (j-0.5)*np.pi/TH                                            # D at phi=0, N at TH
            else: nu = (j-0.5)*np.pi/TH                                                 # N at 0, D at TH
            # complex form g = zeta^nu, Re/Im give cos/sin
            # use analytic function zeta^nu with this branch: zeta = rho e^{i phi}
            if not deriv:
                g = rho**nu*np.exp(1j*nu*phi)
            else:
                g = nu*rho**(nu-1)*np.exp(1j*(nu-1)*phi)*np.exp(-1j*psi)  # d/dz
            if neu0 and not neu1:
                # N at phi=0, D at phi=TH: sin(nu(TH-phi)) = Im(e^{i nu TH} conj?) -> use Re/Im general combination: keep both
                cols.append(g)
            else:
                cols.append(g)
    for m in range(deg+1):
        if deriv: cols.append(m*z**(m-1) if m>0 else 0*z)
        else: cols.append(z**m)
    return np.array(cols).T
def solve2(a, nk, deg, Ms):
    _, z, isN, n = build(10, deg, Ms)
    B = cbasis(z, nk, deg); Bd = cbasis(z, nk, deg, True)
    A_D = np.hstack([B.real, -B.imag]); G=Bd*n[:,None]; A_N=np.hstack([G.real,-G.imag])
    A = np.where(isN[:,None], A_N, A_D)
    r2=np.abs(z)**2; xn=(z*np.conj(n)).real
    rhs=np.where(isN, -dup_dr_over_r(r2,a)*xn, -up(r2,a))
    sc=np.linalg.norm(A,axis=0); sc[sc==0]=1
    c,*_=np.linalg.lstsq(A/sc,rhs,rcond=None); c/=sc
    res=np.max(np.abs(A@c-rhs)); nc=B.shape[1]; cc=c[:nc]+1j*c[nc:]
    a0,b0=sides[NEU]; L=abs(b0-a0); xg,wg=np.polynomial.legendre.leggauss(30)
    br=np.array([0.0]+list(0.5*0.15**np.arange(25,0,-1))+[0.5]); br=np.concatenate([br,1-br[::-1][1:]])
    tot=0.0
    for lo,hi in zip(br[:-1],br[1:]):
        t=lo+(hi-lo)*(xg+1)/2; w=wg*(hi-lo)/2; zz=a0+(b0-a0)*t
        tot+=np.sum(w*(up(np.abs(zz)**2,a)+(cbasis(zz,nk,deg)@cc).real))*L
    return tot,res
if __name__ == "__main__":
    from scipy.optimize import brentq
    for nk, deg, Ms in [(18, 50, 300), (22, 60, 300), (26, 70, 400)]:
        a = brentq(lambda a: solve2(a, nk, deg, Ms)[0] - math.exp(a), -3, 0, xtol=1e-16, rtol=1e-15)
        print("MPS", nk, deg, repr(a), flush=True)
    for Np in [100, 140]:
        a = brentq(lambda a: solve(a, Np, 2*Np//3 + 10, 3*Np + 40)[0] - math.exp(a), -3, 0, xtol=1e-16, rtol=1e-15)
        print("lightning", Np, repr(a), flush=True)
