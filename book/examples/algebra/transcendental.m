# 4.2 Transcendental equations and algebraic collapse
Reduce[Sin[x] == 1/2 && 0 <= x < 2 Pi, x]
Reduce[Cos[x] == 0, x]
Reduce[E^x == 5, x]
Reduce[Log[x] == 2, x]
Simplify[Root[#^2 - 2 &, 2]^2]
Simplify[Root[#^2 - 2 &, 1] + Root[#^2 - 2 &, 2]]
