# Section 4.5.10: NSolve -- domains, polynomial systems, and a transcendental fallback.
NSolve[x^4 - 1 == 0, x]
NSolve[x^4 - 1 == 0, x, Reals]
NSolve[{x^2 + y^2 == 1, x^3 - y^3 == 2}, {x, y}]
NSolve[E^x - x == 7, x, Reals]
