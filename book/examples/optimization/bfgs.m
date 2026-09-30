# Several variables: BFGS on Rosenbrock's function, with monitors and a supplied gradient.
rosen = (1 - x)^2 + 100 (y - x^2)^2;
FindMinimum[rosen, {{x, -1.2}, {y, 1}}]
D[rosen, {{x, y}}]
steps = 0; FindMinimum[rosen, {{x, -1.2}, {y, 1}}, StepMonitor :> steps++]; steps
evals = 0; FindMinimum[rosen, {{x, -1.2}, {y, 1}}, EvaluationMonitor :> evals++]; evals
FindMinimum[rosen, {{x, -1.2}, {y, 1}}, MaxIterations -> 10]
FindMinimum[rosen, {{x, -1.2}, {y, 1}}, Gradient -> {-2 (1 - x) - 400 x (y - x^2), 200 (y - x^2)}]
FindMinimum[Sin[x] Sin[2 y], {x, y}]
