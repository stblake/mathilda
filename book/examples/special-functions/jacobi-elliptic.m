# Jacobi elliptic functions, parameter m = k^2
{JacobiSN[u, 0], JacobiCN[u, 0], JacobiDN[u, 0]}
{JacobiSN[u, 1], JacobiCN[u, 1], JacobiDN[u, 1]}
JacobiSN[0.7, 0.3]
N[JacobiSN[7/10, 3/10], 30]
# The two Pythagorean identities, at 30 digits
s = JacobiSN[7/10`30, 3/10]; c = JacobiCN[7/10`30, 3/10]; d = JacobiDN[7/10`30, 3/10];
{s^2 + c^2, 3/10 s^2 + d^2}
# sn reaches 1 at the quarter period K(m)
JacobiSN[EllipticK[0.3], 0.3]
D[JacobiSN[u, m], u]
D[JacobiAmplitude[u, m], u]
# sn is the inverse of the incomplete elliptic integral F
InverseJacobiSN[1, m]
InverseJacobiSN[JacobiSN[0.7, 0.3], 0.3]
# The amplitude winds continuously, 2 Pi per period 4K
JacobiAmplitude[10., 1.7]
JacobiSN[0.7 + 0.2 I, 0.3]
