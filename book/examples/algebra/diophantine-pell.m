# Pell's equation x^2 - 61 y^2 == 1 — three views, transcript display for §4.2.11.
Solve[x^2 - 61 y^2 == 1 && x > 0 && y > 0 && x < 10^10, {x, y}, Integers]
Solve[x^2 - 61 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers]
FindInstance[x^2 - 61 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers]
