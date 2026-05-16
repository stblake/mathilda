# Series Expansions in Mathilda

In [previous](./Mathilda_patterns_blog.md) [notes](./Mathilda_functional_blog.md) I toured pattern matching and functional programming in Mathilda. This time I want to look at one of the more mathematically ambitious pieces of the system: the `Series` function, which computes Taylor, Laurent and Puiseux expansions of a wide class of symbolic inputs. As usual, every example below was typed straight into the REPL; the `Out[]` lines are the interpreter's actual output.

Mathilda aims to be a faithful recreation of Mathematica's core evaluator and pattern matcher in a small C99 codebase, with the calculus capabilities of a Derive-sized CAS [1]. `Series` is where several of its subsystems — arithmetic, pattern matching, rule rewriting, derivatives, `Apart`, partial fractions — come together.

## The basics

`Series[f, {x, x0, n}]` expands `f` as a power series in `(x - x0)` up to and including order $n$. The result is displayed with an explicit `O[x - x0]^(n+1)` term so you always know where the truncation is.

```mathematica
In[1]:= Series[Exp[x], {x, 0, 6}]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + 1/720 x^6 + O[x]^7

In[2]:= Series[Sin[x], {x, 0, 9}]
Out[2]= x - 1/6 x^3 + 1/120 x^5 - 1/5040 x^7 + 1/362880 x^9 + O[x]^10

In[3]:= Series[Cos[x], {x, 0, 8}]
Out[3]= 1 - 1/2 x^2 + 1/24 x^4 - 1/720 x^6 + 1/40320 x^8 + O[x]^9
```

Coefficients are exact rationals — no floating-point artifacts. The generating function $1/(1-x)$ just gives $\sum x^n$:

```mathematica
In[4]:= Series[1/(1-x), {x, 0, 6}]
Out[4]= 1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + O[x]^7
```

Less obvious generating functions come out cleanly too. The Fibonacci numbers are the coefficients of $1/(1 - x - x^2)$:

```mathematica
In[5]:= Series[1/(1-x-x^2), {x, 0, 8}]
Out[5]= 1 + x + 2 x^2 + 3 x^3 + 5 x^4 + 8 x^5 + 13 x^6 + 21 x^7 + 34 x^8 + O[x]^9
```

and the Bernoulli numbers are (up to $n!$) the coefficients of $x/(e^x - 1)$. The odd-index ones after $B_1$ are zero, so we see a string of missing powers:

```mathematica
In[6]:= Series[x/(Exp[x]-1), {x, 0, 6}]
Out[6]= 1 - 1/2 x + 1/12 x^2 - 1/720 x^4 + 1/30240 x^6 + O[x]^7
```

## Log, inverse trig, and other transcendentals

The usual suspects are all there:

```mathematica
In[7]:= Series[Log[1+x], {x, 0, 6}]
Out[7]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + 1/5 x^5 - 1/6 x^6 + O[x]^7

In[8]:= Series[ArcTan[x], {x, 0, 7}]
Out[8]= x - 1/3 x^3 + 1/5 x^5 - 1/7 x^7 + O[x]^8

In[9]:= Series[ArcSin[x], {x, 0, 7}]
Out[9]= x + 1/6 x^3 + 3/40 x^5 + 5/112 x^7 + O[x]^8

In[10]:= Series[Sinh[x], {x, 0, 7}]
Out[10]= x + 1/6 x^3 + 1/120 x^5 + 1/5040 x^7 + O[x]^8

In[11]:= Series[Tanh[x], {x, 0, 7}]
Out[11]= x - 1/3 x^3 + 2/15 x^5 - 17/315 x^7 + O[x]^8

In[12]:= Series[ArcTanh[x], {x, 0, 7}]
Out[12]= x + 1/3 x^3 + 1/5 x^5 + 1/7 x^7 + O[x]^8
```

`Tan[x]` is secretly harder than `Sin[x]` or `Cos[x]` — its coefficients involve the Bernoulli numbers — but it comes out the same way:

```mathematica
In[13]:= Series[Tan[x], {x, 0, 7}]
Out[13]= x + 1/3 x^3 + 2/15 x^5 + 17/315 x^7 + O[x]^8
```

## Compositions, products, and algebraic combinations

Series form an algebra: you can add, multiply, and compose them, and `Series` applies those operations for you automatically.

Composition:

```mathematica
In[14]:= Series[Exp[Sin[x]], {x, 0, 6}]
Out[14]= 1 + x + 1/2 x^2 - 1/8 x^4 - 1/15 x^5 - 1/240 x^6 + O[x]^7

In[15]:= Series[Log[Cos[x]], {x, 0, 8}]
Out[15]= -1/2 x^2 - 1/12 x^4 - 1/45 x^6 - 17/2520 x^8 + O[x]^9

In[16]:= Series[Log[1+Sin[x]], {x, 0, 5}]
Out[16]= x - 1/2 x^2 + 1/6 x^3 - 1/12 x^4 + 1/24 x^5 + O[x]^6

In[17]:= Series[Exp[Exp[x]-1], {x, 0, 5}]
Out[17]= 1 + x + x^2 + 5/6 x^3 + 5/8 x^4 + 13/30 x^5 + O[x]^6
```

The last example computes the exponential generating function of the Bell numbers $1, 1, 2, 5, 15, 52, \dots$; multiplying coefficients by $n!$ recovers them.

Products distribute through the series algebra without any hand-holding:

```mathematica
In[18]:= Series[Exp[x] Cos[x], {x, 0, 6}]
Out[18]= 1 + x - 1/3 x^3 - 1/6 x^4 - 1/30 x^5 + O[x]^7

In[19]:= Series[(Sin[x]+Cos[x])^2, {x, 0, 6}]
Out[19]= 1 + 2 x - 4/3 x^3 + 4/15 x^5 + O[x]^7
```

A classic identity check: $\sin^2 + \cos^2 \equiv 1$, and the series machinery confirms it to any order you like. The `O[x]^{11}` term survives because that is the accuracy we asked for:

```mathematica
In[20]:= Series[Sin[x]^2 + Cos[x]^2, {x, 0, 10}]
Out[20]= 1 + O[x]^11
```

Another classical identity is the cancellation in $\tan(\sin x) - \sin(\tan x)$, whose first seven Taylor coefficients vanish. Mathilda finds the leading $x^7/30$ and the next correction:

```mathematica
In[21]:= Series[Tan[Sin[x]] - Sin[Tan[x]], {x, 0, 9}]
Out[21]= 1/30 x^7 + 29/756 x^9 + O[x]^10
```

## Laurent series: negative powers

When the input has a pole at the expansion point, `Series` produces a Laurent expansion automatically. The `1/x` term is part of the truncation window, not an error:

```mathematica
In[22]:= Series[1/Sin[x], {x, 0, 5}]
Out[22]= 1/x + 1/6 x + 7/360 x^3 + 31/15120 x^5 + O[x]^6

In[23]:= Series[Cot[x], {x, 0, 5}]
Out[23]= 1/x - 1/3 x - 1/45 x^3 - 2/945 x^5 + O[x]^6

In[24]:= Series[1/(Exp[x]-1), {x, 0, 5}]
Out[24]= 1/x - 1/2 + 1/12 x - 1/720 x^3 + 1/30240 x^5 + O[x]^6
```

Compound Laurent inputs — products, sums, compositions — work the same way:

```mathematica
In[25]:= Series[Sin[x]/(1-Cos[x]), {x, 0, 5}]
Out[25]= 2/x - 1/6 x - 1/360 x^3 - 1/15120 x^5 + O[x]^6
```

## Puiseux series: fractional powers

Branch cuts and radicals are treated as Puiseux expansions — power series in $(x - x_0)^{1/q}$ for some integer $q$. The simplest case is a square root:

```mathematica
In[26]:= Series[Sqrt[1+x], {x, 0, 5}]
Out[26]= 1 + 1/2 x - 1/8 x^2 + 1/16 x^3 - 5/128 x^4 + 7/256 x^5 + O[x]^6

In[27]:= Series[(1+x)^(1/3), {x, 0, 4}]
Out[27]= 1 + 1/3 x - 1/9 x^2 + 5/81 x^3 - 10/243 x^4 + O[x]^5
```

When we expand `ArcSin[x]` *at a branch point* $x = 1$ the answer is genuinely a Puiseux series in $\sqrt{x - 1}$, and Mathilda finds the correct $\pi/2$ constant term and the imaginary-unit prefactor from analytic continuation:

```mathematica
In[28]:= Series[ArcSin[x], {x, 1, 3}]
Out[28]= 1/2 Pi + (-I) Sqrt[2] Sqrt[x - 1]
       + (1/12*I) Sqrt[2] (x - 1)^(3/2)
       + (-3/160*I) Sqrt[2] (x - 1)^(5/2) + O[x - 1]^(7/2)
```

Symbolic fractional powers in the input also work — the binomial series in $\alpha$ is tabulated lazily:

```mathematica
In[29]:= Series[(1+x)^a, {x, 0, 4}]
Out[29]= 1 + a x + 1/2 a (-1 + a) x^2 + 1/6 a (-2 + a) (-1 + a) x^3
       + 1/24 a (-3 + a) (-2 + a) (-1 + a) x^4 + O[x]^5
```

## Expansion away from zero

Any point works, not just $x = 0$. Around $x = 1$:

```mathematica
In[30]:= Series[Log[x], {x, 1, 5}]
Out[30]= (x - 1) - 1/2 (x - 1)^2 + 1/3 (x - 1)^3
       - 1/4 (x - 1)^4 + 1/5 (x - 1)^5 + O[x - 1]^6

In[31]:= Series[Exp[x], {x, 1, 3}]
Out[31]= E + E (x - 1) + 1/2 E (x - 1)^2 + 1/6 E (x - 1)^3 + O[x - 1]^4
```

The expansion point can itself be symbolic. With a completely generic function $f$ at a generic point $a$ the answer is literally Taylor's formula with explicit `Derivative` nodes:

```mathematica
In[32]:= Series[f[x], {x, a, 3}]
Out[32]= f[a] + Derivative[1][f][a] (x - a)
       + 1/2 Derivative[2][f][a] (x - a)^2
       + 1/6 Derivative[3][f][a] (x - a)^3 + O[x - a]^4
```

## Expansions at infinity

`Series[f, {x, Infinity, n}]` expands in powers of $1/x$:

```mathematica
In[33]:= Series[(x+1)/(x-1), {x, Infinity, 4}]
Out[33]= 1 + 2/x + 2 (1/x)^2 + 2 (1/x)^3 + 2 (1/x)^4 + O[1/x]^5
```

## Leading-term form

Sometimes you do not want to pick an order in advance — you just want to know the leading non-zero term. The two-argument rule form `Series[f, x -> x0]` returns the first non-zero coefficient with the O-term landing at the next non-trivial power:

```mathematica
In[34]:= Series[Sin[x] - x, x -> 0]
Out[34]= -1/6 x^3 + O[x]^5

In[35]:= Series[(1-Cos[x])/x^2, {x, 0, 6}]
Out[35]= 1/2 - 1/24 x^2 + 1/720 x^4 - 1/40320 x^6 + O[x]^7
```

## Multivariate series

You can expand in several variables at once. Mathilda treats a two-spec call as an iterated expansion: expand in $x$ first, then expand each coefficient in $y$.

```mathematica
In[36]:= Series[Exp[x+y], {x, 0, 3}, {y, 0, 3}]
Out[36]= 1 + y + 1/2 y^2 + 1/6 y^3 + O[y]^4
       + (1 + y + 1/2 y^2 + 1/6 y^3 + O[y]^4) x
       + (1/2 + 1/2 y + 1/4 y^2 + 1/12 y^3 + O[y]^4) x^2
       + (1/6 + 1/6 y + 1/12 y^2 + 1/36 y^3 + O[y]^4) x^3 + O[x]^4

In[37]:= Series[Sin[x]*Exp[y], {x, 0, 3}, {y, 0, 3}]
Out[37]= (1 + y + 1/2 y^2 + 1/6 y^3 + O[y]^4) x
       + (-1/6 - 1/6 y - 1/12 y^2 - 1/36 y^3 + O[y]^4) x^3 + O[x]^4
```

## `Normal` and working with truncated series

`Normal[s]` drops the `O[]` term and returns the honest polynomial approximation. That polynomial can then be used in any of Mathilda's ordinary algebraic machinery:

```mathematica
In[38]:= Normal[Series[Cos[x], {x, 0, 8}]]
Out[38]= 1 - 1/2 x^2 + 1/24 x^4 - 1/720 x^6 + 1/40320 x^8

In[39]:= CoefficientList[Normal[Series[Sin[x], {x, 0, 15}]], x]
Out[39]= {0, 1, 0, -1/6, 0, 1/120, 0, -1/5040, 0, 1/362880, 0,
         -1/39916800, 0, 1/6227020800, 0, -1/1307674368000}
```

A truncated series is stored internally as a structured `SeriesData` expression — a variable, expansion point, coefficient list, two exponent bounds, and a common denominator for fractional exponents. `FullForm` reveals it:

```mathematica
In[40]:= FullForm[Series[Exp[x], {x, 0, 3}]]
Out[40]= SeriesData[x, 0, List[1, 1, Rational[1, 2], Rational[1, 6]], 0, 4, 1]
```

That is the internal representation the implementation manipulates; the pretty-printer is a view on top of it.

## How it works

Here is the inside of `Series`, roughly in the order the evaluator executes.

**1. A single data structure for Taylor, Laurent, and Puiseux.** Internally all three live in one struct (`src/series.c:263`):

```c
typedef struct {
    Expr*   x;          /* expansion variable */
    Expr*   x0;         /* expansion point */
    Expr**  coefs;      /* symbolic Expr* coefficients */
    size_t  coef_count;
    int64_t nmin;       /* leading exponent numerator */
    int64_t order;      /* O-term exponent numerator */
    int64_t den;        /* common denominator of all exponents */
} SeriesObj;
```

The object represents $\sum_i c_i\,(x-x_0)^{(n_{\min}+i)/d} + O\bigl((x-x_0)^{\text{order}/d}\bigr)$. A plain Taylor series has `nmin = 0` and `den = 1`; a Laurent series allows `nmin < 0`; a Puiseux series has `den > 1`. Everything else is the same code path.

**2. Eager truncation with a padded working order.** Unlike FriCAS, which uses lazy infinite streams, Mathilda commits to a finite working order. To protect against loss of precision when a series gets multiplied by a Laurent factor (e.g. `1/Sin[x]^10`) or inverted with a leading root, the working order is the requested order *plus* a small pad (12 when the expansion point is numeric, 2 when it is symbolic, which keeps symbolic coefficient expressions from blowing up).

**3. Dispatch on the head.** The driver `series_expand` recurses on the expression tree (`src/series.c:1417`):

- Inputs free of the expansion variable become a constant series.
- `Plus` and `Times` recurse and combine via `so_add` / `so_mul` (truncated convolution).
- `Power[base, exp]` handles three cases: neither side contains $x$ (constant), only base contains $x$ (binomial / integer-power path), or the exponent contains $x$ (rewritten as `Exp[exp * Log[base]]` and recursed).
- Reciprocal heads (`Sec`, `Csc`, `Cot`, `Sech`, `Csch`, `Coth`) get rewritten to their reciprocal form before expansion.
- Elementary functions (`Exp`, `Log`, `Sin`, `Cos`, `Sinh`, `Cosh`, `Tan`, `Tanh`, inverse trig and inverse hyperbolic) hit dedicated *kernels* — each kernel tabulates the first $N$ scalar Taylor coefficients of that function and composes with the inner series via a truncated Faà di Bruno convolution (`so_compose_scalar_kernel`).
- Anything else falls back to naive Taylor via `D[]`, using the system's own differentiator.

**4. Algebraic fast paths.** Some shapes are so common that a direct rewrite is dramatically cheaper than the generic kernel. The `Log` branch detects inputs of the form $a + b\,x^{p/q}$ (via `series_split_two_term`, a Maxima-inspired structural probe) and rewrites them as

$$\log(a + b\,x^{p/q}) = \log a + \log\!\bigl(1 + (b/a)\,x^{p/q}\bigr),$$

so the inner `Log[1 + ·]` kernel sees a pure monomial and a cheap closed-form binomial series is used. Similar peeling handles vanishing-argument logs by extracting the leading power.

**5. `Apart` preprocessing for rationals.** When the input is a rational function in $x$, `do_series_single` asks `Apart` to partial-fraction it first (`src/series.c:1914`). Each summand is then a simple pole and expands in closed form, which avoids an $O(N \cdot \deg q)$ Newton inversion of the composite denominator. The check is cheap, so non-rational inputs pay nothing.

**6. Prefactor extraction.** Inputs like $x^a\,\exp(x)$ get the symbolic $x^a$ factored out before the inner series is computed, keeping the coefficient arithmetic symbol-free:

```mathematica
In[41]:= Series[x^a Exp[x], {x, 0, 4}]
Out[41]= x^a (1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + O[x]^5)
```

**7. Expansion at Infinity via substitution.** `Series[f, {x, Infinity, n}]` substitutes $x \mapsto 1/u$, expands at $u = 0$, then relabels the result in terms of `Power[x, -1]`. All the difficult work lives in the ordinary at-zero expander.

**8. Branch-point Puiseux expansions.** At regular points of an inverse trig kernel, the at-zero composition works. At branch points (e.g. `ArcSin` at $x = 1$), the kernel returns `NULL` and `so_apply_arc_branch_point` takes over — it builds the Puiseux expansion directly in $\sqrt{x - x_0}$, respecting the correct sign and imaginary-unit prefactor for the chosen principal branch.

**9. Fallback via `D[]`.** For heads that aren't tabulated — or for tabulated kernels whose expansion point isn't special enough to hit — Mathilda falls back to the Taylor-via-derivatives definition, $f(x) = \sum_k f^{(k)}(x_0)\,(x-x_0)^k/k!$, relying on the system's own differentiator. This is why `Series[f[x], {x, a, 3}]` returns a symbolic Taylor polynomial with `Derivative[k][f][a]` coefficients.

## How it compares: FriCAS and Maxima

The `Series` implementation in Mathilda sits between two much larger designs — AXIOM/FriCAS at one end and Maxima at the other. Each made interesting, different choices.

### FriCAS: layered lazy streams

FriCAS is built as a tower of domains:

- `UnivariateTaylorSeries(Coef, var, cen)` — represents series as a lazy infinite `Stream Coef`. Coefficients are demanded on request, and all arithmetic delegates to `StreamTaylorSeriesOperations`.
- `UnivariateLaurentSeries` — a record `{expon: Int, ps: UTS}` meaning $x^{\text{expon}} \cdot \text{ps}(x)$. Negative powers live in `expon`; the Taylor tail handles everything else. Leading zeros get absorbed by `removeZeroes`, which increments `expon` and divides the Taylor part by $x$.
- `UnivariatePuiseuxSeries` — a record `{expon: Rat, lSeries: ULS}` meaning $\text{lSeries}(x^{\text{expon}})$. Binary operations use `ratGcd` to align ramification indices, then `multiplyExponents` to rescale before delegating to the Laurent layer.

Crucially, FriCAS's arithmetic and transcendental functions never truncate. Reciprocals, compositions, even `revert` (series reversion via Lagrange inversion) and Lambert sums are lazy recurrences on streams. If you later want more coefficients, you just ask.

Mathilda collapses the three layers into one `SeriesObj { nmin, order, den }` and commits up front to a finite order. That sacrifices laziness but keeps the C implementation small (2200 lines, one file), keeps memory layout simple, and lets every coefficient be an ordinary Mathilda symbolic `Expr*` without dragging in an abstract ring tower. It also side-steps some classical pitfalls — an eager implementation doesn't have to worry about leading-zero removal or about whether a requested coefficient depends on something not yet computed; the contract is "expand to this order, now."

The two designs also differ on where transcendentals live. In FriCAS, `exp`, `log`, `sin`, … are implemented as lazy stream operations in `StreamTaylorSeriesOperations`, and they apply to any series over any supported ring. In Mathilda, each transcendental has a hand-written Taylor-coefficient table (`kernel_coefs`) and is composed with the inner series via Faà di Bruno; unknown heads drop through to `D[]`-based Taylor. The result is a smaller library of transcendentals but one that integrates naturally with the rest of the system: a coefficient can be `Log[a]` or `Derivative[k][f][a]` or a piecewise-constant without any special machinery.

### Maxima: recursive series-of-series with algebraic fast paths

Maxima has two series packages. `taylor.lisp` is the Taylor-polynomial engine used for most exact expansions. `series.lisp` (the one `powerseries[f, x, x0]` dispatches to) is different: it treats an expansion as a symbolic infinite sum

```
((%sum) term k 0 $inf)
```

and applies rewriting passes until each piece is in closed form. The two passes that matter most:

- `sp2expand` — the fast-path dispatch. `sp2log`, `sp2expt`, tabulated-by-`defprop` handlers for `%sin`, `%cos`, `%exp`, and so on.
- `sratexpnd` — handles rational functions, and if the denominator factors with distinct non-zero roots it uses the distinct-roots theorem directly rather than going through full partial fractions.

The single most influential piece of Maxima's design that shows up in Mathilda is `split-two-term-poly`: the idea that you can probe an expression structurally for the shape $a + b\,x^{p/q}$ without doing a full series expansion, and then rewrite

$$\log(a + b\,x^{p/q}) \;=\; \log a + \log\!\bigl(1 + (b/a)\,x^{p/q}\bigr)$$

so the generic log kernel sees a pure binomial. Mathilda keeps that fast path more or less verbatim — `series_split_two_term` in `src/series.c` is the probe, and the `Log` branch of `series_expand` is the rewrite. The same idea is what makes `Series[Log[Cos[x]], ...]` or `Series[Log[1+Sin[x]], ...]` fall out cleanly.

Where Maxima and Mathilda diverge is in the answer's *shape*. Maxima's `series.lisp` is willing to return an infinite-summation form — symbolic coefficients parametrised by an index $k$. Mathilda always returns a truncated `SeriesData` object with a concrete O-term; that is the price for matching Mathematica-style output.

### A one-paragraph summary of the trade-offs

FriCAS chose generality and laziness: a tower of categories, streams everywhere, transcendentals expressed as lazy recurrences that work over any ring. Maxima chose algebraic rewriting: turn a series request into a symbolic sum and attack it with fast-path lemmas like `split-two-term-poly`. Mathilda sits between them — one flat Puiseux struct, eager truncation, a Mathematica-style `SeriesData` output, but with Maxima's fast paths for `Log`, a Faà di Bruno kernel for each elementary head, and `Apart` + prefactor peeling gluing the whole thing together. It fits in a single C file and never leaves the Mathilda expression tree, which is the design constraint the system cares about most.

[1] https://en.wikipedia.org/wiki/Derive_(computer_algebra_system)
