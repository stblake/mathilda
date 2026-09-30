# DIRECT is deterministic. On a symmetric box its very first sample is the centre.
rast = 20 + x^2 - 10 Cos[2 Pi x] + y^2 - 10 Cos[2 Pi y];
NMinimize[{rast, -5.12 <= x <= 5.12 && -5.12 <= y <= 5.12}, {x, y}, Method -> {"DIRECT", "PostProcess" -> False}]
NMinimize[{rast, -3 <= x <= 5.12 && -4 <= y <= 5.12}, {x, y}, Method -> {"DIRECT", "PostProcess" -> False}]
NMinimize[{rast, -3 <= x <= 5.12 && -4 <= y <= 5.12}, {x, y}, Method -> "DIRECT"]
