# 4.3.1 Vector calculus -- Grad, Div, Curl, Laplacian
Grad[x^2 y + y^2 z, {x, y, z}]
Grad[{x y, y z, z x}, {x, y, z}]
Div[{x^2, y^2, z^2}, {x, y, z}]
Curl[{-y, x, 0}, {x, y, z}]
Laplacian[x^2 + y^2 + z^2, {x, y, z}]
Grad[k q/r, {r, t, p}, "Spherical"]
Div[{r^2, 0, 0}, {r, t, p}, "Spherical"]
