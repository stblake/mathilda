# Section 4.6.9: A complete textbook RSA, in one continuous session.
p = 61
q = 53
n = p q
phi = EulerPhi[n]
phi == (p - 1) (q - 1)
e = 17
GCD[e, phi]
ExtendedGCD[e, phi]
d = PowerMod[e, -1, phi]
Mod[e d, phi]
msg = 42
c = PowerMod[msg, e, n]
PowerMod[c, d, n]
PowerMod[c, d, n] == msg
