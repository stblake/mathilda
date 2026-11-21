# Certificates Mathilda already makes for itself.
s = Reduce`SOSCertificate[x^4 + y^4 - x^2 y^2 + 1, {}, {x, y}];
Expand[Total[(#1 (#2 . #3 . #2)) & @@@ s] - (x^4 + y^4 - x^2 y^2 + 1)]
Reduce`SOSCertificate[x^2 - 1, {}, {x}]
PossibleZeroQ[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]
Minimize[x^4 - 2 x^2, x]
