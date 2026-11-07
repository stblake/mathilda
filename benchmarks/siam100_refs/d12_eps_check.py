# D12: limiting ODE (with the adiabatic-invariant correction) vs the full eps-system
# (eps = 1e-2, 5e-3, 2.5e-3) and vs the naive constrained pendulum; first hitting times of y=0.
import numpy as np
from scipy.integrate import solve_ivp
def full(alpha, eps):
    def f(t,u):
        x,y,vx,vy=u
        g=x*x+y*y-1; c=1+alpha*x*x
        # U = y + eps^-2/2 c g^2
        Ux = eps**-2*(0.5*2*alpha*x*g*g + c*g*2*x)
        Uy = 1 + eps**-2*(c*g*2*y)
        return [vx,vy,-Ux,-Uy]
    ev=lambda t,u:u[1]; ev.terminal=True; ev.direction=-1
    s=solve_ivp(f,[0,50],[0,1,1,1],method='DOP853',rtol=1e-11,atol=1e-12,events=ev)
    return s.t_events[0][0] if len(s.t_events[0]) else None
def lim(alpha, corr=True):
    def f(t,u):
        th,w=u
        acc=-np.cos(th)
        if corr: acc+= alpha*np.sin(th)*np.cos(th)/(2*np.sqrt(1+alpha*np.cos(th)**2))
        return [w,acc]
    ev=lambda t,u:np.sin(u[0]); ev.terminal=True; ev.direction=-1
    s=solve_ivp(f,[0,50],[np.pi/2,-1],method='DOP853',rtol=1e-12,atol=1e-13,events=ev)
    return s.t_events[0][0] if len(s.t_events[0]) else None
for a in [0.5,1,3]:
    print(a, lim(a), lim(a,False), [full(a,e) for e in [1e-2,5e-3,2.5e-3]])
