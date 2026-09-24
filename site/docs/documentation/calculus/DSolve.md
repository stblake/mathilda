# DSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DSolve[eqn, y, x] solves a differential equation for the function y with`**

**`DSolve[eqn, y[x], x] returns {{y[x] -> expr}} with the solution as an`**

**`DSolve[{eqn1, ...}, {y1, ...}, x] solves a system; DSolve[eqn, u, {x, y}]`**

<details>
<summary>Notes</summary>

independent variable x, returning {{y -\> Function\[{x}, ...\]}}. expression in x. solves a partial differential equation. Initial/boundary conditions are given as equations at points, e.g. y\[0\] == a, y'\[x0\] == b; DSolve fits the generated constants C\[k\] to them. Every returned branch is verified by back-substitution before it is kept. DSolve is a cascade polyalgorithm: DSolve\[eqn, y, x, Method -\> "\<Name\>"\] dispatches directly to one method (strict, no fallback), and every method is also callable on its own as DSolve\`\<Name\>\[...\] with its own ?-documentation. Methods (tried roughly in this order): First order: Quadrature              y^(n)\[x\] == f(x): integrate n times. LinearFirstOrder        y'\[x\] + p y == q: integrating factor. Separable               y'\[x\] == g(x)h(y). Bernoulli               y'\[x\] == A y\[x\] + B y^n\[x\]. Homogeneous             y'\[x\] == F(y\[x\]/x). Exact                   M + N y'\[x\] == 0, M\_y == N\_x (+ mu(x)/mu(y) factor search). Clairaut                y\[x\] == x y'\[x\] + f(y'\[x\]) (+ singular envelope). Riccati                 y'\[x\] == q0 + q1 y\[x\] + q2 y\[x\]^2: linearize to 2nd-order. Lagrange                y == x phi(y'\[x\]) + psi(y'\[x\]) (parametric solution). Chini                   y'\[x\] == f y\[x\]^n + g y\[x\] + h (reducible-to-autonomous). Abel                    y'\[x\] == f3 y\[x\]^3 + f2 y\[x\]^2 + f1 y\[x\] + f0. FirstOrderSubstitution  y'\[x\] == F(a x + b y\[x\] + c). Factorable              factor in y'\[x\], solve each factor. NthAlgebraic            algebraic in the top derivative. AlmostLinear            f(x)g(y\[x\])y'\[x\] + k(x)l(y\[x\]) + m(x) == 0. LinearCoefficients      y'\[x\] == (a1 x + b1 y\[x\] + c1)/(a2 x + b2 y\[x\] + c2). SeparableReduced        x y'\[x\]/y == G(x^n y\[x\]). LieSymmetry (LieGroup)  heuristic Lie point-symmetry backstop. Linear constant-coefficient: LinearConstantCoefficients   any order, homogeneous + inhomogeneous. UndeterminedCoefficients     tidy particular for poly/exp/sinusoid forcing. Linear variable-coefficient / 2nd order: EulerCauchy             a\_n x^n y\[x\]^(n) + ... == g. ExactODE                higher-order exact linear equation. NormalForm              y''\[x\] + P y'\[x\] + Q y\[x\] -\> z'' == r z (Kovacic prerequisite). SpecialFunctionForm     Airy/Bessel/Kummer/Gauss recognizers. Kovacic                 Liouvillian solutions of z''==r z. FrobeniusSeries         series about a regular singular point. PowerSeries             series about an ordinary point. FirstOrderPowerSeries   order-1 power series (pinned-only). OperatorFactor (DFactor)  factor a linear operator into first-order factors. Nonlinear higher-order: ReductionOfOrder        y''\[x\] == F(x,y'\[x\]) (missing y). AutonomousReduction     y''\[x\] == f(y\[x\],y'\[x\]) (missing x). Liouville               y''\[x\] + g(y\[x\]) y'\[x\]^2 + h(x) y'\[x\] == 0. Systems (nfun\>1): DecoupleSystem          each equation involves one function. TriangularSystem        dependency graph is a DAG. LinearFirstOrderSystem  Y'==A Y+b(x), any constant A (Jordan e^{Ax}). First-order PDE: PDELinearFirstOrder     a u\_{v1}+b u\_{v2}+c u==f (characteristics). Options: GeneratedParameters (constant head, default C), Assumptions, Method, IncludeSingularSolutions.  Attributes: Protected.

</details>

## Examples (22)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (10)

```mathematica
In[1]:= DSolve[y''[x] + 4 y'[x] + 5 y[x] == 0, y[x], x]
Out[1]= {{y[x] -> C[1] Cos[x] E^(-2 x) - C[2] Sin[x] E^(-2 x)}}
```

BVP

```mathematica
In[2]:= DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi/2] == 1}, y[x], x]
Out[2]= {{y[x] -> Sin[x]}}
```

Kovacic Case 1

```mathematica
In[3]:= DSolve[y''[x] - (x^2 + 3) y[x] == 0, y[x], x]
Out[3]= {{y[x] -> (C[1] (x^2)^(3/4) E^(1/2 x^2))/Sqrt[x] + (C[2] (x^2)^(1/4) HypergeometricPFQ[{1}, {1/2}, x^2] E^(-1/2 x^2))/Sqrt[x]}}
```

No closed form → series

```mathematica
In[4]:= DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x]
Out[4]= {{y[x] -> C[1] + C[2] x + -1/6 C[1] x^3 + -1/12 C[2] x^4 + 1/120 C[1] x^5 + 1/180 (C[1] + C[2]) x^6 + O[x]^7}}
```

Airy

```mathematica
In[5]:= DSolve[y''[x] - x y[x] == 0, y[x], x]
Out[5]= {{y[x] -> C[1] AiryAi[x] + C[2] AiryBi[x]}}
```

Bessel

```mathematica
In[6]:= DSolve[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y[x], x]
Out[6]= {{y[x] -> C[1] BesselJ[2, x] + C[2] BesselY[2, x]}}
```

First-order substitution

```mathematica
In[7]:= DSolve[y'[x] == (x + y[x])^2, y[x], x]
Out[7]= {{y[x] -> -x - Tan[-C[1] - x]}}
```

2nd-order autonomous, nonlinear

```mathematica
In[8]:= DSolve[y[x] y''[x] == y'[x]^2, y[x], x]
Out[8]= {{y[x] -> C[1] E^(C[2] x)}}
```

Defective + singular

```mathematica
In[9]:= DSolve[{y'[t] == 0, x'[t] + y[t] == 0}, {y[t], x[t]}, t]
Out[9]= {{y[t] -> C[1], x[t] -> C[2] - C[1] t}}
```

Inconsistent BVP

```mathematica
In[10]:= DSolve[{y''[x] + y[x] == 0, y[0] == 1, y[Pi] == 1}, y, x]
Out[10]= {}
```

### Scope (12)

Transport

```mathematica
In[11]:= DSolve[D[u[t,x],t] + c D[u[t,x],x] == 0, u, {t,x}]
Out[11]= {{u -> Function[{t, x}, C[1][-c t + x]]}}
```

```mathematica
In[12]:= DSolve[D[u[x,y],x] + 3 D[u[x,y],y] + u[x,y] == 1, u, {x,y}]
Out[12]= {{u -> Function[{x, y}, E^(-x) (E^x + C[1][-3 x + y])]}}
```

Wave/d'Alembert

```mathematica
In[13]:= DSolve[D[u[t,x],{t,2}] == c^2 D[u[t,x],{x,2}], u, {t,x}]
Out[13]= {{u -> Function[{t, x}, C[1][-Sqrt[c^2] t + x] + C[2][Sqrt[c^2] t + x]]}}
```

Laplace/elliptic

```mathematica
In[14]:= DSolve[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, u, {x,y}]
Out[14]= {{u -> Function[{x, y}, C[1][-I x + y] + C[2][I x + y]]}}
```

```mathematica
In[15]:= DSolve`SeparationOfVariables[D[u[x,t],t] == D[u[x,t],{x,2}], u, {x,t}]
Out[15]= {{u -> Function[{x, t}, E^(-C[3] t) (C[1] E^(-1/2 Sqrt[-4 C[3]] x) + C[2] E^(1/2 Sqrt[-4 C[3]] x))]}}

In[16]:= PDEClassify[D[u[t,x],{t,2}] == c^2 D[u[t,x],{x,2}] /. c -> 2, u, {t,x}]
Out[16]= "Hyperbolic"

In[17]:= {PDEClassify[D[u[x,t],t] == D[u[x,t],{x,2}], u, {x,t}], PDEClassify[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, u, {x,y}]}
Out[17]= {"Parabolic", "Elliptic"}

In[18]:= DSolve[{D[u[x,t],{t,2}] == 4 D[u[x,t],{x,2}], u[x,0] == Sin[x], Derivative[0,1][u][x,0] == 0}, u, {x,t}]
Out[18]= {{u -> Function[{x, t}, 1/2 (Sin[-2 t + x] + Sin[2 t + x])]}}

In[19]:= DSolve[y'[x] + y[x] == a Sin[x], y[x], x]
Out[19]= {{y[x] -> E^(-x) (C[1] + (-1/4 - 1/4*I) a E^((1 + I) x) + (-1/4 + 1/4*I) a E^((1 - I) x))}}

In[20]:= DSolve[{y'[x] == 3 y[x], y[0] == 5}, y[x], x]
Out[20]= {{y[x] -> 5 E^(3 x)}}

In[21]:= DSolve[{y'[x] == -3 y[x]^2, y[1] == 2}, y[x], x]
Out[21]= {{y[x] -> 1/(-5/2 + 3 x)}}

In[22]:= DSolve[y''[x] == 7, y, x]
Out[22]= {{y -> Function[{x}, C[1] + C[2] x + 7/2 x^2]}}
```

## Algorithm

dsolve.c — DSolve dispatcher (cascade polyalgorithm).

Mirrors src/calculus/integrate.c: a Method-option enum selects either the automatic cascade (try each method until one returns a non-NULL result) or a

```text
single pinned method (strict, no fallback).  A per-command fail-memo keyed on
```

eval_toplevel_id() suppresses the fixed-point loop's redundant re-entry on a problem the deterministic cascade already declined, and g_dsolve_depth distinguishes the outermost user call from internal recursions.

The shared problem substrate (parse / verify / fit / assemble) is in dsolve_common.c; each method is one file src/calculus/dsolve_<method>.c.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Integrate](../../calculus/Integrate/), [HoldAll](../../expression-information/HoldAll/), [PolynomialQ](../../algebra/PolynomialQ/), [Solve](../../solutions-of-equations/Solve/), [TimeConstrained](../../time-and-date/TimeConstrained/), [Piecewise](../../control-flow/Piecewise/), [Log](../../elementary-functions/Log/), [Hypergeometric1F1](../../special-functions/Hypergeometric1F1/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)
- Tests: [`tests/test_dsolve_m12_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m12_stress.c)
- Tests: [`tests/test_dsolve_m14_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m14_stress.c)
- Tests: [`tests/test_dsolve_m17_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m17_stress.c)
