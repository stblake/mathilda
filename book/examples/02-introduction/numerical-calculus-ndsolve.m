# Van der Pol oscillator y'' + mu (y^2 - 1) y' + y = 0 (mu = 3): a nonlinear ODE
# with no elementary closed form, solved as an InterpolatingFunction on [0, 20].
sol = NDSolve[{y''[x] + 3 (y[x]^2 - 1) y'[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 20}]
First[y[5.0] /. sol]
First[y[20.0] /. sol]
