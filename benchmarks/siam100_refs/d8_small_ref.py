"""Sanity/reduced rows of d8, method 1 (mpmath, dps 80):
  c_3(R^1)^2 = (2 pi)^-1 Int_R dx/(1+x^2+x^4+x^6)                      (= 1/4)
  c_2(R^2)^2 = (2 pi)^-2 Int_{R^2} dx dy/(1+x^2+y^2+x^4+x^2 y^2+y^4)
             = (2 pi)^-2 (1/2) Int_0^{2 pi} Int_0^inf dt/(1 + t + q t^2) d theta,
  q = 1 - cos^2 sin^2 (polar, t = r^2), radial integral in closed form
  Int_0^inf dt/(1+t+q t^2) = 2/s (pi/2 - atan(1/s)), s = sqrt(4q - 1)."""
from mpmath import mp, mpf, quad, sqrt, pi, cos, sin, atan, inf, nstr
mp.dps = 80
c3 = sqrt(quad(lambda x: 1/(1+x**2+x**4+x**6), [-inf, 0, inf]) / (2*pi))
print("c3R1", nstr(c3, 75))
def rad(th):
    q = 1 - cos(th)**2*sin(th)**2
    s = sqrt(4*q - 1)
    return 2/s*(pi/2 - atan(1/s))
c2 = sqrt(quad(rad, [0, pi/2]) * 4 / 2 / (2*pi)**2)
print("c2R2", nstr(c2, 75))
