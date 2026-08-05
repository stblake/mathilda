# Compilation and auto-compilation

Mathilda evaluates by rewriting expression trees. That is what makes it a
*symbolic* system — `x + x` becomes `2 x`, `D[Sin[x], x]` becomes `Cos[x]` — but
it is also why a tight numeric loop is slow: every `s = s + i` builds and frees a
small tree. When the work is purely numeric, Mathilda can leave the tree
rewriter behind and run the same computation as **machine bytecode over machine
numbers**. There are two ways in:

- **`Compile[]`** — you ask for it. It turns a typed body into a
  `CompiledFunction` you call like any other function.
- **Automatic compilation** — Mathilda does it behind your back. Numeric
  `Do`/`For`/`While`/`Map`/`Nest`/`Fold` loops, and the sample-point bodies of
  `Plot`, `Table`, `NIntegrate` and friends, are compiled and run on the fast
  path without a single keystroke from you. The switch is `$AutoCompilation`,
  on by default.

Both are bound by one contract, and it is the thing to hold onto through this
whole tutorial: **the fast path must give the interpreter's answer, or step
aside.** When a body falls outside the compilable subset, compilation does not
fail — it silently runs the interpreter instead. Correct, slower, and (as we'll
see) visible if you ask.

A third mechanism sits alongside these and is easy to confuse with them:
**automatic packed arrays** (`$AutoArrayPacking`) store a large list of machine
numbers as one dense buffer instead of thousands of boxed values. It is not
compilation — no bytecode is involved — but it is the other half of "run numeric
work at machine speed", and for whole-array math it is the more important half.
We'll meet it in [§5](#5-whole-array-math-and-autoarraypacking).

Every transcript below was produced by the actual Mathilda binary.

!!! note "The numbers in this tutorial"
    Timings were measured on an **Intel Core i9-9880H @ 2.30 GHz** (8 cores,
    macOS), Mathilda 0.032 at `-O3`, against **Python 3.11.15 + NumPy 2.4.4** on
    the same machine. Wall-clock time (`AbsoluteTiming` in Mathilda,
    `time.perf_counter` in Python), minimum of several repetitions after a
    warm-up. Absolute times are hardware-specific; the **ratios** are the
    portable part. Reproduce the Mathilda side with the `bench` helper defined in
    [§8](#8-measuring-honestly).

## 1. Your first CompiledFunction

`Compile[{args}, body]` builds a `CompiledFunction`. Each argument is a
`{name, type}` pair — `_Integer`, `_Real`, `_Complex`, or an array type we'll
reach later.

```mathematica
In[1]:= sq = Compile[{{x, _Integer}}, x^2 + 2 x + 1];

In[2]:= sq
Out[2]= CompiledFunction[{x}, x^2 + 2 x + 1]

In[3]:= sq[5]
Out[3]= 36
```

It behaves like any other function, but it is not a rewrite rule — it is a
compiled program that takes a machine integer in and gives a machine integer
back. The results are exactly the interpreter's:

```mathematica
In[4]:= Table[sq[k], {k, 0, 5}]
Out[4]= {1, 4, 9, 16, 25, 36}

In[5]:= Table[k^2 + 2 k + 1, {k, 0, 5}]
Out[5]= {1, 4, 9, 16, 25, 36}
```

### Did it actually compile? `CompileDiagnostics`

This is the single most useful habit when working with `Compile[]`, and the
reason is worth stating plainly: **a body outside the compilable subset does not
tell you.** The object is still built; calling it silently runs the interpreter.
You get the right answer, more slowly, with no diagnostic. `CompileDiagnostics`
asks the question directly, without building the object:

```mathematica
In[6]:= CompileDiagnostics[{{x, _Integer}}, x^2 + 2 x + 1]
Out[6]= {"Compiled" -> True, "ResultType" -> "Integer", "Instructions" -> 5, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 7}
```

`"Compiled" -> True`, a result type of `Integer`, and five bytecode
instructions — the optimiser folded the seven the emitter first produced down to
five.

### Reading the bytecode: `CompilePrint`

`CompileDiagnostics` says *how many* instructions; `CompilePrint` says *which*.
It takes the compiled object and prints the program:

```
In[7]:= CompilePrint[sq]

Signature   CompiledFunction[{x : Integer}, x^2 + 2 x + 1]
Arguments   1
              R0   : Integer      x
Result      R1 : Integer
Registers   3 scalar, 0 array, 0 tile   (frame 3 slots)
Program     5 instructions, 0 CSE

    0  POWI_I     R1, R0, 2                       R1 = R0^2
    1  MUL_IK     R2, R0, 2                       R2 = R0 * 2
    2  ADD_I      R1, R1, R2                      R1 = R1 + R2
    3  ADD_IK     R1, R1, 1                       R1 = R1 + 1
    4  RET        R1                              return R1
```

The `_I` suffixes are the point: this is **integer** register arithmetic.
`POWI_I` is integer power by repeated multiplication, `MUL_IK` folds the constant
`2` into the instruction, `ADD_IK` adds the constant `1`. No `Expr` nodes, no
`pow()` call, no allocation. That is the whole difference from the interpreter,
and it is why the loops below run tens of times faster.

## 2. Machine reals, and where `Compile` earns its keep

Change the argument type to `_Real` and the same machinery emits floating-point
bytecode:

```mathematica
In[8]:= poly = Compile[{{x, _Real}}, x^3 - 2 x + 1];

In[9]:= poly
Out[9]= CompiledFunction[{x}, x^3 - 2 x + 1]

In[10]:= poly[1.5]
Out[10]= 1.375

In[11]:= CompileDiagnostics[{{x, _Real}}, x^3 - 2 x + 1]
Out[11]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 6, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 12}
```

A single call is too cheap to be worth compiling — the win is in a *loop*. The
canonical case is a **scalar recurrence**, where each step depends on the last so
there is nothing to vectorise. The logistic map is the textbook example:

```mathematica
In[12]:= logistic = Compile[{{x0, _Real}, {n, _Integer}}, Module[{x = x0}, Do[x = 3.7 x (1 - x), {n}]; x]];

In[13]:= logistic[0.5, 20]
Out[13]= 0.264924
```

A `Module` with a `Do` loop compiles as one program: the accumulator `x` lives in
a register and the loop body is straight-line bytecode. Ten million iterations:

| path | time (n = 10⁷) | vs interpreter |
|---|---:|---:|
| Mathilda interpreter (`$AutoCompilation = False`) | 9.44 s | 1× |
| Mathilda auto-compiled (default) | 0.172 s | **55×** |
| Mathilda explicit `Compile[]` | 0.125 s | **75×** |
| Python 3.11 scalar loop | 0.445 s | — |
| NumPy | — | *no vectorisation for a sequential recurrence* |

Two things to read off this. First, `Compile[]` is **75× faster than the
interpreter** and **3.5× faster than a Python 3.11 loop**. Second, the row that
matters most is the empty NumPy cell: a recurrence `xₙ₊₁ = f(xₙ)` cannot be
expressed as a whole-array operation, because each element needs the one before
it. This is exactly the shape where an array library has nothing to offer and a
compiler is the only way to reach machine speed. **Compile's home turf is the
code that will not vectorise.**

## 3. Machine integers, and the overflow contract

Integer loops compile the same way, into the integer bytecode of [§1](#1-your-first-compiledfunction):

```mathematica
In[14]:= isum = Compile[{{n, _Integer}}, Module[{s = 0}, Do[s = s + i, {i, 1, n}]; s]];

In[15]:= isum[100]
Out[15]= 5050
```

There is one subtlety that the parity contract forces, and it is worth
understanding because it is where a naive compiler would go wrong. A machine
`int64` accumulator can **overflow**, and the interpreter answers an overflowing
integer sum by promoting to an arbitrary-precision bignum — it never wraps. So
the compiled path must not wrap either. It doesn't: an operation that overflows
`int64` abandons the compiled run and hands the whole loop back to the
interpreter, which produces the exact bignum.

```mathematica
In[16]:= Module[{s = 0}, Do[s = s + i^2, {i, 1, 10000000}]; s]
Out[16]= 333333383333335000000

In[17]:= Sum[i^2, {i, 1, 10000000}]
Out[17]= 333333383333335000000
```

That sum is `3.3 × 10²⁰`, far past the `int64` ceiling of `9.2 × 10¹⁸` — and the
answer is the exact bignum, identical to `Sum`, not a silently wrapped machine
integer. You pay for the overflow by re-running in the interpreter, but only when
one actually happens; the common no-overflow case runs at machine speed.

For the sum `Σᵢ₌₁ⁿ i` at `n = 10⁷` (which stays inside `int64`):

| path | time (n = 10⁷) | vs interpreter |
|---|---:|---:|
| Mathilda interpreter (`$AutoCompilation = False`) | 6.21 s | 1× |
| Mathilda auto-compiled (default) | 0.070 s | **89×** |
| Mathilda explicit `Compile[]` | 0.072 s | **86×** |
| Python 3.11 `for`-loop | 0.395 s | — |
| NumPy `arange(...).sum()` | 0.019 s | — |

Here NumPy's SIMD reduction (`0.019 s`) beats the compiled *scalar loop*
(`0.072 s`) — and that is the honest framing, not a defeat. A running sum **is** a
reduction, and a reduction vectorises; written as one it is faster in every
system, including Mathilda's own `Total[Range[10^7]]`. The loop here exists to
show integer register code and the overflow contract, not to argue a scalar loop
should beat a vector reduction. The lesson that keeps recurring: **vectorise what
vectorises; compile what won't.**

## 4. Auto-compilation: the win you already have

You do not have to write `Compile[]` to get most of this. When `$AutoCompilation`
is on — the default — Mathilda compiles the body of a numeric `Do`, `For`,
`While`, `Map`, `Nest`, `Fold`, and the per-sample bodies inside `Plot`, `Table`,
`NIntegrate`, `NSum`, `FindRoot` and the plot samplers, behind your back.

```mathematica
In[18]:= $AutoCompilation
Out[18]= True
```

The two rows labelled "auto-compiled" in the tables above are the whole story: an
ordinary interpreted-looking loop, written with no thought of compilation, runs
**55× (logistic) and 89× (integer sum)** faster than it would with the switch
off — reaching within a small factor of an explicit `Compile[]`. Here is the
same logistic recurrence written as a plain function, with no `Compile`:

```mathematica
In[19]:= logisticInterp[x0_, n_] := Module[{x = x0}, Do[x = 3.7 x (1 - x), {n}]; x];
```

`$AutoCompilation` exists so you can see the two paths side by side. Turn it off,
and the *same* function runs through the interpreter:

```mathematica
In[20]:= $AutoCompilation = False
Out[20]= False

In[21]:= logisticInterp[0.5, 20]
Out[21]= 0.264924

In[22]:= $AutoCompilation = True
Out[22]= True

In[23]:= logisticInterp[0.5, 20]
Out[23]= 0.264924
```

Same answer, both ways — which is the contract. The switch changes speed and
nothing else. And it agrees with the explicitly-compiled version too:

```mathematica
In[24]:= logisticInterp[0.5, 20] == logistic[0.5, 20]
Out[24]= True
```

!!! tip "When you want to *check* a compiled result"
    Because a compiled body is contracted to give the interpreter's answer,
    `$AutoCompilation = False` is a one-line way to confirm it. If a numeric loop
    ever gives you a result you doubt, flip the switch and compare; a difference
    is a bug, and a match is reassurance. Set the environment variable
    `MATHILDA_COMPILE_DIAG=1` to have the internal auto-compilers print to stderr
    whenever one of them falls back to the interpreter.

## 5. Whole-array math and `$AutoArrayPacking`

Scalar recurrences are Compile's territory. **Whole-array** math — apply a
function to every element, sum a vector, multiply matrices — is a different
regime, and there the decisive lever is not compilation but **packing**.

A large list of machine numbers is stored, by default, as one dense buffer (a
*packed array*) rather than thousands of individually-boxed values. It is still
an ordinary `List` — same head, same printed form, same elements — and only
`NDArrayQ` tells the two apart:

```mathematica
In[25]:= NDArrayQ[Range[1000000]]
Out[25]= True
```

`$AutoArrayPacking` controls it, and like `$AutoCompilation` it changes storage
and speed, never answers:

```mathematica
In[26]:= $AutoArrayPacking = False
Out[26]= False

In[27]:= NDArrayQ[Range[1000000]]
Out[27]= False

In[28]:= $AutoArrayPacking = True
Out[28]= True

In[29]:= NDArrayQ[Range[1000000]]
Out[29]= True
```

The difference this makes is enormous, because on a packed buffer the operation
runs as one tight machine loop (or a BLAS call), while on an unpacked list every
element is a boxed `Real` reached through a pointer. On `N = 10⁶`:

| kernel | unpacked (interp) | packed (interp) | packing speedup | `Compile[]` | NumPy 2.4.4 |
|---|---:|---:|---:|---:|---:|
| `Total[v]` | 53.3 ms | 0.24 ms | **222×** | — | 0.35 ms |
| `v . v` (dot) | 422 ms | 0.21 ms | **2010×** | — | 0.22 ms |
| `Total[Sin[v]^2 + Cos[v]^2]` | 1767 ms | 8.8 ms | **201×** | 10.3 ms | 11.2 ms |

Three things this table says, each worth pausing on.

**Packing alone buys 200×–2000×.** No `Compile`, no rewriting — the same
expression, on data that happens to be a buffer instead of a list.

**Packed Mathilda is right on NumPy's heels.** `v . v` is `0.21 ms` against
NumPy's `0.22 ms` — both dispatch to the same Accelerate `ddot`. The transcendental
sum is `8.8 ms` against `11.2 ms`. Running the same BLAS and the same
vectorised kernels, a packed Mathilda array is a machine-bandwidth reference, not
a symbolic system limping along.

**`Compile[]` barely beats the packed interpreter here** — `10.3 ms` compiled
against `8.8 ms` interpreted, i.e. compiling is *marginally slower*. That is not a
bug; it is the whole point. The packed interpreter already runs one machine loop
per array operation, so there is almost nothing for a compiler to add. Compiling
a whole-array expression fuses the passes into one, which saves the intermediate
buffers — but when the kernels are already vectorised and multi-threaded, the
saving is in the noise. **For code that vectorises, reach for packing, not
`Compile`.** Save `Compile` for the loops that don't — the recurrences of
[§2](#2-machine-reals-and-where-compile-earns-its-keep).

You can still compile array code when you want a single fused pass, and the array
argument type is `{name, _Real, rank}`:

```mathematica
In[30]:= vk = Compile[{{v, _Real, 1}}, Total[Sin[v]^2 + Cos[v]^2]];

In[31]:= vk
Out[31]= CompiledFunction[{v}, Total[Sin[v]^2 + Cos[v]^2]]

In[32]:= vk[{0.0, 1.0, 2.0}]
Out[32]= 3.0

In[33]:= CompileDiagnostics[{{v, _Real, 1}}, Total[Sin[v]^2 + Cos[v]^2]]
Out[33]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 15, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 15}
```

For a much deeper worked example of compiled indexed-array code — a finite-difference
PDE stencil, where indexing *doesn't* vectorise and `Compile` wins by 50×+, with a
full bytecode reading and a comparison against Wolfram Language — see
[the wave-equation study](https://github.com/stblake/mathilda/blob/main/docs/compile_example/COMPILE_EXAMPLE.md).

## 6. Compiling with associations

`Compile[]` accepts an **association** as an argument — a feature Mathematica's
`Compile` does not have — declared with the `_Association` type. Inside the body,
association operations run as native kernels with **no evaluator calls at
runtime**.

A read-only parameter bag: look keys up and combine them.

```mathematica
In[34]:= mean = Compile[{{a, _Association}}, (Lookup[a, "x"] + Lookup[a, "y"])/2];

In[35]:= mean[<|"x" -> 3.0, "y" -> 5.0|>]
Out[35]= 4.0

In[36]:= CompileDiagnostics[{{a, _Association}}, (Lookup[a, "x"] + Lookup[a, "y"])/2]
Out[36]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 5, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 8}
```

The set-algebra and higher-order operations compile too — `KeyDrop`, `KeyTake`,
`Counts`, `Map`, `Select`, `Append` — and compose to any depth:

```mathematica
In[37]:= kd = Compile[{{a, _Association}}, Total[Values[KeyDrop[a, "b"]]]];

In[38]:= kd[<|"a" -> 1.0, "b" -> 2.0, "c" -> 3.0|>]
Out[38]= 4.0

In[39]:= sqmap = Compile[{{a, _Association}}, Total[Values[Map[Function[u, u^2], a]]]];

In[40]:= sqmap[<|"a" -> 1.0, "b" -> 2.0, "c" -> 3.0|>]
Out[40]= 14.0
```

The interesting case for performance is a loop of **lookups by a runtime key**.
Each probe is an `O(1)` hash lookup into the association's index, done natively:

```mathematica
In[41]:= asum = Compile[{{a, _Association}}, Module[{s = 0.}, Do[s = s + Lookup[a, i], {i, 1, 5}]; s]];

In[42]:= asum[AssociationThread[Range[5], N[Range[5]^2]]]
Out[42]= 55.0

In[43]:= CompileDiagnostics[{{a, _Association}}, Module[{s = 0.}, Do[s = s + Lookup[a, i], {i, 1, 5}]; s]]
Out[43]= {"Compiled" -> True, "ResultType" -> "Real", "Instructions" -> 14, "CommonSubexpressions" -> 0, "InstructionsUnoptimized" -> 16}
```

Scaling that to a 10⁴-entry association, looked up 10⁷ times:

| path | time (10⁷ lookups) | vs interpreter |
|---|---:|---:|
| Mathilda interpreter (`$AutoCompilation = False`) | 8.48 s | 1× |
| Mathilda explicit `Compile[]` | 0.356 s | **24×** |
| Python 3.11 `dict` | 0.548 s | — |
| NumPy | — | *no associative array* |

Compiled association lookups run **24× faster than the interpreter and 1.5×
faster than a Python `dict`**, with no `evaluate()` in the inner loop — and, as
the empty NumPy cell notes, this is a data structure array libraries do not have.

## 7. The cliff, and how to see it

Everything above works because the body was inside the compilable subset. What
happens when it isn't? Not an error — a fall back to the interpreter. And the
subset is a **cliff, not a slope**: *one* unsupported head anywhere in the body
sends the *entire* body to the interpreter, not just that node. So the habit from
[§1](#1-your-first-compiledfunction) is not optional. `CompileDiagnostics` names
the innermost subexpression that stopped it:

```mathematica
In[44]:= CompileDiagnostics[{{x, _Real}}, Sin[x] + BarnesG[x]]
Out[44]= {"Compiled" -> False, "Reason" -> "no machine lowering for this head at these argument types", "Subexpression" -> "BarnesG[x]"}
```

`Sin[x]` compiles fine; `BarnesG[x]` has no machine kernel, and it takes the
whole `Sin[x] + BarnesG[x]` down with it. The report points at `BarnesG[x]`, the
actual culprit, not the enclosing `Plus`. The same is true inside control flow:

```mathematica
In[45]:= CompileDiagnostics[{{x, _Real}}, If[x > 0, Sqrt[x], BarnesG[x]]]
Out[45]= {"Compiled" -> False, "Reason" -> "no machine lowering for this head at these argument types", "Subexpression" -> "If[x > 0, Sqrt[x], BarnesG[x]]"}
```

When you hit a cliff, the options are: rewrite the offending piece into the
subset, split the body so the compilable part still compiles, or accept the
interpreter for that case. What you must not do is *assume* it compiled — which
is the entire reason `CompileDiagnostics` exists.

## 8. Measuring honestly

A few traps bite everyone who benchmarks a compiler, and every table above was
built to avoid them.

- **Use wall-clock time, not `Timing[]`.** `Timing` reports CPU time summed over
  threads; on a threaded path (BLAS, a parallel map) it reads as *N times slower*
  than the truth. Use `AbsoluteTiming`, whose first element is elapsed seconds:

  ```mathematica
  SetAttributes[bench, HoldAll];
  bench[e_] := Min[Table[First[AbsoluteTiming[e]], {5}]];
  ```

- **Force the result.** A benchmark that discards its output can end up measuring
  the discarding. Sum the array, take its last element — make the answer
  observable.

- **Put an error column next to every timing.** A wrong answer is usually a
  *fast* answer — an empty loop, a body that compiled to nothing. Every table
  here was checked for value parity (Mathilda `==` NumPy `==` Python) before its
  times were trusted; a timing row means nothing until the answers agree.

- **Remember which number moved.** A compiled-vs-interpreted speedup is a
  statement about *two* numbers. Automatic packed arrays have made Mathilda's
  interpreter dramatically faster on array code, so `Compile`'s advantage *over
  it* has shrunk — which is good news about the interpreter, not bad news about
  the compiler.

## Where to go next

- [`docs/compile_example/COMPILE_EXAMPLE.md`](https://github.com/stblake/mathilda/blob/main/docs/compile_example/COMPILE_EXAMPLE.md)
  — one problem in depth: a compiled 2-D wave-equation solver, the full bytecode,
  and a comparison against Wolfram Language's two `Compile` backends.
- [`docs/design/performance.md`](https://github.com/stblake/mathilda/blob/main/docs/design/performance.md)
  — the broad sweep: Mathilda vs Mathematica 14.0 vs NumPy across 70+ HPC
  kernels.
- The [documentation centre](../documentation/index.md) — the reference pages
  for `Compile`, `CompiledFunction`, `CompileDiagnostics`, `CompilePrint`,
  `$AutoCompilation`, and `$AutoArrayPacking`.
