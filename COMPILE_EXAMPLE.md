# A worked example: `Compile[]` and the 2-D wave equation

This is a tutorial on Mathilda's `Compile[]`, built around one problem that is
big enough to be worth compiling: an explicit finite-difference solver for the
wave equation on a square, marched in time.

It covers what the compiler can and cannot express, how to tell which of those
you are in, and what the speed actually is — measured against the same program
in the interpreter, against Wolfram Language 14.0 (`wolframscript`, both its
bytecode and its native-C compilation targets), and against `NDSolve`.

Every number below was measured on one machine (Apple M-series, macOS 15.6,
Mathilda at `-O3`) on 2026-07-27. The scripts are reproduced in full so they can
be re-run.

**Contents**

1. [The problem, and how we know the answer is right](#the-problem)
2. [The scheme](#the-scheme)
3. [Writing it: the interpreted version](#the-interpreted-version)
4. [Compiling it](#compiling-it)
5. [Is it actually compiled? `CompileDiagnostics`](#is-it-actually-compiled)
6. [What `Compile[]` does to the numbers](#results)
7. [Against Wolfram Language](#against-wolfram-language)
8. [Against `NDSolve`](#against-ndsolve)
9. [What is and is not in the compilable subset](#the-compilable-subset)
10. [Measurement traps found while writing this](#measurement-traps)

---

## The problem

Solve

$$u_{tt} = u_{xx} + u_{yy} \qquad \text{on } [0,1]^2,$$

with $u = 0$ on the boundary, $u(0,x,y) = \sin \pi x \, \sin \pi y$, and
$u_t(0,x,y) = 0$. The exact solution is a single standing mode,

$$u(t,x,y) = \sin \pi x \, \sin \pi y \, \cos\!\left(\sqrt 2\, \pi t\right).$$

We march to $T = 0.5$.

Having a closed-form solution matters for more than tidiness. A solver can be
wrong in two independent ways — the *discretisation* can be a poor approximation
of the PDE, and the *implementation* can be a poor realisation of the
discretisation — and a comparison against $u(t,x,y)$ alone cannot tell them
apart. So we check both, separately:

- **Physical error** — distance from the continuum solution. Should fall like
  $O(h^2)$. This validates the scheme.
- **Discrete error** — distance from the *exact solution of the difference
  equations*, derived in [The scheme](#the-scheme). Should be at roundoff, and stay there for hundreds
  of steps. This validates the implementation, and it is the one that catches
  an off-by-one in a stencil.

The second check is the reason this tutorial can make claims about correctness
at all. It is worth setting up whenever you write a scheme by hand.

---

## The scheme

Uniform grid $x_i = (i-1)h$, $y_j = (j-1)h$, $h = 1/(n-1)$; the standard
second-order explicit stencil,

$$u^{k+1}_{i,j} = 2u^k_{i,j} - u^{k-1}_{i,j} + \lambda^2\left(u^k_{i+1,j} + u^k_{i-1,j} + u^k_{i,j+1} + u^k_{i,j-1} - 4u^k_{i,j}\right),$$

with Courant number $\lambda = c\,\Delta t/h$. In 2-D this is stable for
$\lambda \le 1/\sqrt2$; we use $\lambda = 1/2$ throughout, so
$\Delta t = h/2$ and reaching $T = 0.5$ takes $n - 2$ steps from the two
starting levels.

**The exact discrete solution.** Write $\varphi_{ij} = \sin \pi x_i \sin \pi y_j$.
The discrete 5-point Laplacian has $\varphi$ as an eigenvector:

$$\varphi_{i+1,j} + \varphi_{i-1,j} + \varphi_{i,j+1} + \varphi_{i,j-1} - 4\varphi_{ij} = \left(4\cos \pi h - 4\right)\varphi_{ij}.$$

Substituting $u^k = \varphi \cos(k\theta)$ into the scheme and using
$\cos((k{+}1)\theta) + \cos((k{-}1)\theta) = 2\cos\theta\cos(k\theta)$ gives a
single condition,

$$\boxed{\;\cos\theta = 1 - 2\lambda^2\left(1 - \cos \pi h\right).\;}$$

So if we *start* the march from $u^0 = \varphi$ and $u^1 = \varphi\cos\theta$,
then $u^k = \varphi\cos(k\theta)$ holds **exactly**, for every $k$, in exact
arithmetic. Any deviation is floating-point roundoff or a bug. That is the
discrete reference.

```mathematica
h    = 1./(n - 1);
lam  = 0.5;
phi  = Table[Sin[Pi (i - 1) h] Sin[Pi (j - 1) h], {i, 1, n}, {j, 1, n}];
ct   = 1 - 2 lam^2 (1 - Cos[Pi h]);     (* = Cos[theta] *)
u0   = phi;
u1   = ct phi;
steps = n - 2;                          (* (steps+1) dt == 0.5 *)
```

---

## The interpreted version

One time step, written the way you would write it without a compiler in mind —
build the new grid with `Table`:

```mathematica
istep[up_, uc_, lam_, n_] :=
  Table[
    If[i == 1 || i == n || j == 1 || j == n, 0.,
      2 uc[[i, j]] - up[[i, j]] +
        lam^2 (uc[[i + 1, j]] + uc[[i - 1, j]] +
               uc[[i, j + 1]] + uc[[i, j - 1]] - 4 uc[[i, j]])],
    {i, 1, n}, {j, 1, n}];

up = u0; uc = u1;
Do[tmp = istep[up, uc, lam, n]; up = uc; uc = tmp, {steps}];
Max[Abs[Flatten[uc - Cos[(steps + 1) ArcCos[ct]] phi]]]
```

At `n = 41` this runs in **4.12 s** and reports a discrete error of
`2.78e-14` — the scheme is right.

It is also unusably slow. $41^2 \times 39 \approx 65{,}000$ stencil evaluations took four
seconds, or about 63 µs each, essentially all of it spent building and
destroying `Expr` trees.

---

## Compiling it

The compiled step is the same arithmetic, but written as an in-place update of a
grid rather than as a `Table` that constructs one:

```mathematica
step = Compile[{{up, _Real, 2}, {uc, _Real, 2}, {lam, _Real}},
  Module[{n = Length[uc], un = uc},
    Do[
      un[[i, j]] = 2 uc[[i, j]] - up[[i, j]] +
        lam^2 (uc[[i + 1, j]] + uc[[i - 1, j]] +
               uc[[i, j + 1]] + uc[[i, j - 1]] - 4 uc[[i, j]]),
      {i, 2, n - 1}, {j, 2, n - 1}];
    un]];
```

Four things in that fragment are worth pointing at.

**`{up, _Real, 2}` declares a rank-2 array argument.** The third element of an
argument spec is the rank. A `List` passed in is packed into a flat machine
buffer at the boundary and unpacked on the way out; an `NDArray` is used
directly, with no copy.

**`un = uc` inside `Module` makes an array LOCAL, and it is a copy.** This is
what lets the body write into `un` at all. Argument arrays are *borrowed* — the
caller still owns them — so writing through one is not in the compilable subset,
deliberately: for a `List` argument, packed into a temporary at the boundary, the
write would vanish without a trace. Copying into a local is what the
interpreter's value semantics do anyway.

Because `un` starts as a copy of `uc`, the boundary ring is already zero and the
loop only has to touch the interior. That is why the compiled body has no `If`
and the interpreted one does.

**`un[[i, j]] = ...` writes the buffer in place.** No new grid is allocated per
element. The subscripts are resolved and range-checked one axis at a time, which
matters: `un[[1, n + 5]]` is inside the buffer as a linear offset and would
quietly write into the next row, so the check cannot be on the flattened index.

**`Do[..., {i, ...}, {j, ...}]` takes several iterators**, nesting them with the
last varying fastest, exactly as the interpreter does.

Driving it:

```mathematica
up = NDArray[u0]; uc = NDArray[u1];
Do[tmp = step[up, uc, lam]; up = uc; uc = tmp, {steps}];
```

`NDArray[...]` is worth the keystrokes. Passing plain `List`s works and gives the
same answer, but each call then packs two grids in and unpacks one out; at
`n = 41` that is the difference between 7.2 ms and 18.3 ms. Keeping the state
packed across steps costs one conversion at the start rather than three per step.

---

## Is it actually compiled?

**The single most useful habit when working with `Compile[]`.** A body outside
the compilable subset does not fail loudly — the object is still built, and
calling it silently runs the interpreter. You get the right answer, 10–40× more
slowly, with no diagnostic. And the subset is a *cliff, not a slope*: one
unsupported head anywhere costs the entire body, not just that node.

`CompileDiagnostics` answers the question directly:

```mathematica
In[1]:= CompileDiagnostics[{{up, _Real, 2}, {uc, _Real, 2}, {lam, _Real}},
          Module[{n = Length[uc], un = uc},
            Do[un[[i, j]] = 2 uc[[i, j]] - up[[i, j]] +
                 lam^2 (uc[[i + 1, j]] + uc[[i - 1, j]] +
                        uc[[i, j + 1]] + uc[[i, j - 1]] - 4 uc[[i, j]]),
               {i, 2, n - 1}, {j, 2, n - 1}];
            un]]

Out[1]= {"Compiled" -> True, "ResultType" -> "Array", "Instructions" -> 72,
         "CommonSubexpressions" -> 1, "InstructionsUnoptimized" -> 92}
```

72 bytecode instructions per interior point, one repeated subexpression hoisted
(`uc[[i, j]]`, which appears twice), and the optimiser removed 20 of the 92
instructions the emitter first produced.

When a body does *not* compile, the report names the **innermost** subexpression
that stopped it — not the enclosing expression, which is rarely the culprit:

```mathematica
In[2]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + BarnesG[x]]

Out[2]= {"Compiled" -> False,
         "Reason" -> "no machine lowering for this head at these argument types",
         "Subexpression" -> "BarnesG[x]"}
```

For the internal auto-compiled builtins (`Plot`, `NIntegrate`, `Table`, …), the
environment variable `MATHILDA_COMPILE_DIAG=1` prints the same report to stderr
whenever one of them falls back.

---

## Results

Same problem, same scheme, same starting levels, both arms at top level, timed
in-process. `steps = n - 2` so every row reaches $T = 0.5$.

### Compiled vs interpreted, `n = 41`

| | time | discrete error | speedup |
|---|---:|---:|---:|
| interpreted (`Table`) | 4.119 s | 2.78e-14 | 1× |
| compiled, `List` arguments | 0.0183 s | 2.78e-14 | 225× |
| **compiled, `NDArray` state** | **0.00724 s** | 2.78e-14 | **569×** |

Both arms produce a bit-identical discrete error, which is the check that the
compiled program is running the same scheme and not a subtly different one.

### Scaling of the compiled solver

| $n$ | steps | time | discrete error | physical error |
|---:|---:|---:|---:|---:|
| 41 | 39 | 0.0069 s | 2.8e-14 | 2.27e-4 |
| 101 | 99 | 0.0987 s | 7.8e-14 | 3.63e-5 |
| 153 | 151 | 0.357 s | 4.2e-13 | 1.57e-5 |
| 201 | 199 | 0.818 s | 7.9e-14 | 9.09e-6 |
| 401 | 399 | 6.62 s | 1.8e-12 | 2.27e-6 |
| **641** | **639** | **27.3 s** | 9.7e-12 | 8.87e-7 |

The physical error falls by a factor of ~4 whenever $n$ doubles — second order,
as the scheme promises. The discrete error stays at roundoff across 639 steps on
a 641 × 641 grid ($2.6 \times 10^8$ stencil updates), which is the implementation check.

The last row is the headline: **a 641 × 641 grid marched 639 steps in 27 s** —
$2.6 \times 10^8$ stencil updates, about 104 ns each. Extrapolating the interpreted
per-update cost measured at `n = 41` (62.8 µs), the same march interpreted would
take on the order of four and a half hours. That one is an extrapolation, not a
measurement: it was not run.

---

## Against Wolfram Language

The identical source, run under `wolframscript` (Wolfram Language 14.0). WL's
`Compile` has two backends: `"WVM"` (its bytecode virtual machine, the default)
and `"C"` (generate C, compile it, load it as a library). A C compiler *was*
available and the C target genuinely produced a `LibraryFunction` — this was
checked rather than assumed, with `FreeQ[stepC, LibraryFunction]`.

Every configuration reports the same discrete error to the last digit, which is
what licenses comparing the times at all.

| $n$ | Mathilda `Compile` | WL `Compile` (WVM) | WL `Compile` (C) | Mathilda advantage |
|---:|---:|---:|---:|---:|
| 41 | 0.0069 s | 0.0147 s | 0.0137 s | 2.1× |
| 101 | 0.0987 s | 0.1877 s | 0.2248 s | 1.9× |
| 153 | 0.357 s | 0.655 s | — | 1.8× |
| 201 | 0.818 s | 1.501 s | 1.823 s | 1.8× |
| 401 | 6.62 s | 12.17 s | 14.87 s | 1.8× |
| **641** | **27.3 s** | **52.4 s** | 65.7 s | **1.9×** |

Two things to say honestly about this.

**Mathilda's compiled code is about 1.8–2.1× faster here, and the ratio is
stable across two and a half orders of magnitude of problem size.** For a
stencil written this way — indexed reads and writes into rank-2 machine arrays,
inside nested counted loops — Mathilda's register VM is doing meaningfully less
work per element.

**WL's native-C target is not faster than its own bytecode VM on this body**, and
past `n = 101` it is consistently *slower* (14.87 s vs 12.17 s at `n = 401`).
That is a useful data point for Mathilda's own roadmap, which lists a native
backend as the next codegen step: it suggests the cost in a tensor-heavy kernel
sits in array element access rather than in bytecode dispatch, so replacing
dispatch alone buys little. Whatever a native backend is worth, this measurement
says it is not automatic.

### And the interpreters

Compiled-vs-interpreted ratios only mean something relative to the interpreter
they are measured against, so here are both, at `n = 41`, top level:

| | time |
|---|---:|
| Mathilda, interpreted | 4.119 s |
| Wolfram, `Table` auto-compilation **off** | 0.326 s |
| Wolfram, default (`Table` auto-compiled) | 0.056 s |
| Wolfram, explicit `Compile` (WVM) | 0.0147 s |
| **Mathilda, explicit `Compile`** | **0.0072 s** |

Read fairly, this says two different things at once. **Mathilda's *interpreter*
is about 12× slower than Wolfram's on this array-heavy code** — that is a real
gap and this tutorial is not going to pretend otherwise. And **Mathilda's
*compiled* code is about 2× faster than Wolfram's.** The 569× headline figure in
[Results](#results) is therefore partly a compliment to `Compile[]` and partly a comment on the
interpreter it is being compared with; the cross-system compiled numbers are the
ones that isolate the compiler.

A detail worth knowing if you benchmark WL yourself: by default WL
auto-compiles a `Table` whose body is numeric once it has 250 or more elements
(`SystemOptions["CompileOptions"]`, `TableCompileLength`). The 0.326 s row above
had that switched off; the 0.056 s row is the default, and it is *already
compiled code* despite looking like an interpreter baseline.

---

## Against `NDSolve`

The obvious question about all of this is why write a scheme by hand when
`NDSolve` will solve the same PDE.

```mathematica
NDSolve[{D[u[t, x, y], {t, 2}] == D[u[t, x, y], {x, 2}] + D[u[t, x, y], {y, 2}],
         u[0, x, y] == Sin[Pi x] Sin[Pi y],
         Derivative[1, 0, 0][u][0, x, y] == 0,
         u[t, 0, y] == 0, u[t, 1, y] == 0,
         u[t, x, 0] == 0, u[t, x, 1] == 0},
        u, {t, 0, 0.5}, {x, 0, 1}, {y, 0, 1}]
```

| method | time | physical error at $T=0.5$ |
|---|---:|---:|
| Mathilda `NDSolve` | 0.031 s | 5.4e-5 |
| Wolfram `NDSolve` | 0.081 s | 3.1e-4 |
| Mathilda compiled FD, `n = 41` | 0.0069 s | 2.27e-4 |
| Mathilda compiled FD, `n = 153` | 0.357 s | 1.57e-5 |

**At matched accuracy, `NDSolve` wins comfortably.** Mathilda's `NDSolve` reaches
5.4e-5 in 0.031 s; getting the hand-written scheme below that needs roughly
`n = 153`, which costs 0.357 s — about 12× more. `NDSolve` uses a higher-order
spatial discretisation and an adaptive time integrator, and on a smooth problem
like this one those are simply better tools than a fixed second-order explicit
stencil with a Courant-limited step.

That is the honest framing for `Compile[]`. It does not make your scheme
competitive with a good library solver on a problem the library solver was built
for. What it does is make the schemes that *aren't* in the library — a custom
flux limiter, an unusual boundary condition, a coupled constraint, a stochastic
term, a bespoke stencil you are researching — run at machine speed instead of
`Expr`-tree speed. The 569× is what you get back when you have to write it
yourself.

(Both `NDSolve` rows use default settings, and the two systems pick different
default spatial grids, so this is a defaults-vs-defaults comparison, not an
algorithmic one.)

---

## The compilable subset

Indexed arrays are new (M3c). What the compiler now accepts:

**Reading `Part` — the full spec vocabulary.** Two different lowerings, chosen by
the shape of the subscript list, but between them they cover everything `Part`
accepts on a dense array:

| form | example | lowering |
|---|---|---|
| one scalar subscript per axis | `v[[i]]`, `m[[i, j]]`, `t[[i, j, k]]` | inline; no allocation |
| counting from the end | `v[[-1]]`, `m[[2, -2]]` | inline |
| a computed index | `v[[k + 1]]`, `v[[2 k]]` | inline |
| span | `v[[2 ;; 5]]`, `v[[1 ;; 7 ;; 2]]`, `v[[-3 ;; -1]]` | delegated |
| `All` | `m[[All, 2]]` | delegated |
| list of positions | `v[[{1, 4, 2}]]` | delegated |
| partial indexing | `m[[2]]` (a row), `t[[2, 3]]` | delegated |
| any mixture | `m[[k, 2 ;; 4]]` | delegated |

The inline path is the fast one — a stencil lives there. The delegated path calls
the same `ndarray_part` the interpreter calls, so the compiled answer *is* the
interpreted answer by construction rather than by agreement, at the cost of
allocating the result.

**Writing `Part`**, with the same spec vocabulary: `u[[i, j]] = x`,
`u[[2 ;; 4]] = 0.`, `u[[All, 2]] = 0.`, `u[[{1, 3}]] = 7.`, and
`u[[1 ;; 3]] = w` from a matching array. The compound forms `+=`, `-=`, `*=`,
`/=` work on a scalar position.

**Creating arrays:** `ConstantArray[v, n]` and `ConstantArray[v, {n1, n2, ...}]`,
at any rank. The rank must be evident from the source (a bare dimension, or the
length of a literal dimension list); the dimensions themselves are ordinary
expressions evaluated per call.

**Array locals:** `Module[{u = ...}, ...]` and `With`, including returning one as
the result.

What is deliberately **not** in the subset, and why:

- **Writing through an argument array.** Borrowed, not owned — see [Compiling it](#compiling-it).
- **`ConstantArray[0, n]` with an *integer* fill.** The interpreter's
  `ConstantArray[0, n]` holds exact integer zeros and an `NDArray` has no integer
  dtype, so compiling it to a float64 buffer would answer *differently*, not just
  faster. `ConstantArray[0., n]` is the compilable spelling.
- **Array-valued `If` branches, `Sum` accumulators, or `Nest` state.** These would
  duplicate a handle without duplicating ownership.
- **`Table` as an array constructor inside a compiled body.** Use
  `ConstantArray` plus a `Do` loop.

The rule behind all of these, and the one to keep in mind generally: *the
compiled path must never answer where the interpreter declines, nor differently.*
When the compiled path cannot honour that, it bails, and you get the interpreter
— correct, and slow, and now visible via `CompileDiagnostics`.

---

## Measurement traps

Four things went wrong while measuring this. All four produced plausible numbers.

**A `Hold`-ed body silently compiled to nothing.** Trying to share one body
between WL's two compilation targets via `Evaluate@ReleaseHold@body` evaluated
`Length[uc]` on the free symbol `uc` first, getting `0`, so `Do[..., {i, 2, -1}]`
never ran. The compiled function returned its input unchanged and looked
*extremely* fast — 0.03 s for what should take 52 s. The error column caught it:
1.6 instead of 1e-11. **Report an error next to every timing.** A wrong answer is
usually a fast answer.

**Mathilda's `Module` costs 3× on this loop.** The identical interpreted march
takes 4.09 s at top level and 12.28 s wrapped in a `Module`. Quoting the second
would have inflated the `Compile[]` speedup from 569× to 1700× on the strength of
an unrelated interpreter overhead. Both arms of a comparison have to be in the
same scoping context. (The overhead itself is a real Mathilda issue, recorded
separately.)

**WL's "interpreted" baseline was already compiled.** `TableCompileLength` is 250
by default, and the grid has 1681 elements. The honest interpreter baseline needs
`SetSystemOptions` to switch that off; it is 6× slower than the default.

**`Timing` is CPU time, not wall clock.** Harmless here — the compiled solver is
single-threaded, and the two agreed to 0.3% over a 33 s sweep — but a threaded
path reports as *N times slower* under `Timing`. It is worth confirming they
agree before trusting either.

---

## Appendix: complete runnable script

```mathematica
step = Compile[{{up, _Real, 2}, {uc, _Real, 2}, {lam, _Real}},
  Module[{n = Length[uc], un = uc},
    Do[un[[i, j]] = 2 uc[[i, j]] - up[[i, j]] +
         lam^2 (uc[[i + 1, j]] + uc[[i - 1, j]] +
                uc[[i, j + 1]] + uc[[i, j - 1]] - 4 uc[[i, j]]),
       {i, 2, n - 1}, {j, 2, n - 1}];
    un]];

run[n_] := Module[{h, lam, phi, ct, steps, up, uc, tmp, t},
  h = 1./(n - 1); lam = 0.5; steps = n - 2;
  phi = Table[Sin[Pi (i - 1) h] Sin[Pi (j - 1) h], {i, 1, n}, {j, 1, n}];
  ct = 1 - 2 lam^2 (1 - Cos[Pi h]);
  up = NDArray[phi]; uc = NDArray[ct phi];
  t = Timing[Do[tmp = step[up, uc, lam]; up = uc; uc = tmp, {steps}];][[1]];
  {"time" -> t,
   "discrete error" -> Max[Abs[Flatten[Normal[uc] - Cos[(steps + 1) ArcCos[ct]] phi]]],
   "physical error" -> Max[Abs[Flatten[Normal[uc] -
      Table[Sin[Pi (i - 1) h] Sin[Pi (j - 1) h] Cos[Sqrt[2] Pi 0.5],
            {i, 1, n}, {j, 1, n}]]]]}];

run[41]
run[641]   (* ~27 s *)
```

---

## See also

- [`docs/design/compile.md`](docs/design/compile.md) — the compiler's design.
- [`docs/design/compile_state.md`](docs/design/compile_state.md) — current state,
  milestones, and the measurement traps found so far.
- [`docs/spec/builtins/control-flow.md`](docs/spec/builtins/control-flow.md) —
  `Compile`, `CompiledFunction`, `CompileDiagnostics`.
