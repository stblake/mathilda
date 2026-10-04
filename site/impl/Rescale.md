---
source: src/list/rescale.c
---
**Algorithm.** `builtin_rescale` maps a value from one range to another:
`Rescale[x, {min, max}]` is `(x - min)/(max - min)`, the three-argument form adds
the affine target range `y0 + (y1 - y0)(x - min)/(max - min)`, and the
one-argument `Rescale[list]` reduces to `Rescale[list, {Min[list], Max[list]}]`
and re-evaluates (for a lone scalar `Min == Max`, giving `Indeterminate`, as in
Mathematica). An argument count outside 1–3 raises `Rescale::argb`.

**No per-element threading.** The affine formula is built *once* with the whole
`x` substituted in, as a `Plus`/`Times`/`Power` tree, and handed to the
evaluator. Because `Plus`, `Times` and `Power` are `Listable`, they thread at
every level — including nested lists — so no explicit per-element recursion is
needed and the arithmetic is identical to it. This is the difference between an
interpreted loop and a buffer op: the old per-element `Rescale[el, range]`
evaluation cost 2.17 s over 10⁶ packed reals; handing the whole array to the
Listable heads reaches the threaded, vectorised ND kernels instead.

**Data structures / limits.** Pure `Expr`-tree construction; nothing in this
file reads an element. `Rescale` is on `pack.c`'s `AWARE` list (it only rewrites
the call) but deliberately not `int64_ok` — `Rescale[Range[10]]` is a list of
exact rationals. Exact input stays exact and symbolic elements pass through. The
range arguments must be two-element `List`s. `ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.
