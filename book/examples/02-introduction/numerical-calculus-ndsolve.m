# Simple harmonic motion y'' = -y, y(0)=1, y'(0)=0; exact solution is cos(x).
sol = NDSolve[{y''[x] == -y[x], y[0] == 1, y'[0] == 0}, y, {x, 0, 10}]
First[y[N[Pi]] /. sol]
First[y[N[2 Pi]] /. sol]
