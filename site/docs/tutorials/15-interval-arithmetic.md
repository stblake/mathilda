# Interval arithmetic

Ordinary floating-point answers a question you did not ask. When Mathilda prints
`Sin[N[Pi]]` as `1.22465*10^-16`, that number is not wrong — it is the exact sine
of the *machine* value of π, which differs from true π in the sixteenth digit. But
it is silent about its own uncertainty: nothing on the screen tells you the true
answer is `0`.

**Interval arithmetic** replaces each number by an enclosing range and each
operation by an operation on ranges, carrying a guarantee through every step: the
computed interval always **contains** the true result. It is the arithmetic of
*validated* computing — proving that a root exists, bounding a function over a
region, or watching rounding error grow instead of trusting it away.

Mathilda's `Interval` is unusual in being both **exact** and **rigorous** at once.
Exact endpoints stay exact and symbolic — `Sin[Interval[{2, 7}]]` keeps the bound
`Sin[2]` in closed form — while inexact endpoints are rounded *outward* so the
enclosure never lies. This tutorial works through the classic examples, from the
dependency problem to chaos, showing off both sides.

Every transcript below was produced by the actual Mathilda binary.

## Constructing intervals

`Interval[{min, max}]` is the closed range between its endpoints. Mathilda
normalises on the way in: endpoints are sorted, a bare number becomes a point, and
overlapping or touching pieces are merged.

```mathematica
In[1]:= Interval[{2, 5}]
Out[1]= Interval[{2, 5}]

In[2]:= Interval[{5, 2}]
Out[2]= Interval[{2, 5}]

In[3]:= Interval[7]
Out[3]= Interval[{7, 7}]

In[4]:= Interval[{1, 3}, {2, 5}]
Out[4]= Interval[{1, 5}]
```

The multi-argument form is a **union** of ranges. When the pieces are disjoint they
are kept separate and sorted — an `Interval` is really a finite union of intervals:

```mathematica
In[5]:= Interval[{0, 1}, {2, 3}]
Out[5]= Interval[{0, 1}, {2, 3}]
```

`Min` and `Max` read off the extreme endpoints:

```mathematica
In[6]:= Min[Interval[{2, 5}]]
Out[6]= 2

In[7]:= Max[Interval[{0, 1}, {2, 3}]]
Out[7]= 3
```

## Arithmetic, kept exact

The four operations thread through intervals. With exact endpoints the results stay
exact — no floating point is introduced:

```mathematica
In[1]:= Interval[{1, 6}] + Interval[{0, 2}]
Out[1]= Interval[{1, 8}]

In[2]:= Interval[{2, 3}] - Interval[{5, 6}]
Out[2]= Interval[{-4, -2}]

In[3]:= 3 Interval[{1, 4}]
Out[3]= Interval[{3, 12}]

In[4]:= Interval[{1/2, 2}] Interval[{-3, 4}]
Out[4]= Interval[{-6, 8}]
```

Squaring an interval that straddles zero pins the lower bound to the exact `0`,
because a square is never negative:

```mathematica
In[5]:= Interval[{-2, 5}]^2
Out[5]= Interval[{0, 25}]
```

Dividing by an interval that contains zero produces an **unbounded, disjoint**
result — the reciprocal blows up on either side of 0, and Mathilda represents that
exactly, with symbolic infinities and exact rational endpoints:

```mathematica
In[6]:= 1/Interval[{-2, 5}]
Out[6]= Interval[{-Infinity, -1/2}, {1/5, Infinity}]
```

## The dependency problem

Interval arithmetic has one famous subtlety, and it is worth meeting head on. Each
operation treats its operands as **independent** ranges, even when they are the same
quantity. So subtracting an interval from itself does **not** give zero:

```mathematica
In[1]:= Interval[{-1, 1}] - Interval[{-1, 1}]
Out[1]= Interval[{-2, 2}]
```

The result `[-2, 2]` is correct as a set operation — it encloses every difference
`x - y` with `x, y` in `[-1, 1]` — but it is pessimistic if you meant a *single* `x`.
This is the **dependency problem**, and it means the form in which you write an
expression matters. Compare the true square with the "independent" product:

```mathematica
In[2]:= Interval[{-1, 1}]^2
Out[2]= Interval[{0, 1}]

In[3]:= Interval[{-1, 1}] Interval[{-1, 1}]
Out[3]= Interval[{-1, 1}]
```

`^2` knows both factors are the same `x`, so it returns the *exact* range `[0, 1]`.
Written as a product, the two copies are treated independently and the answer widens
to `[-1, 1]`. The consequence for evaluating a whole expression is that you get a
**guaranteed enclosure**, but not necessarily a tight one:

```mathematica
In[4]:= Interval[{-1, 2}]^3 - 2 Interval[{-1, 2}] + 1
Out[4]= Interval[{-4, 11}]
```

Every value of `x^3 - 2x + 1` for `x` in `[-1, 2]` lies in `[-4, 11]` — that is
certified. (The tight range is about `[-0.09, 5]`; recovering it needs subdivision,
splitting `[-1, 2]` into many small intervals and taking the union.)

## Functions: symbolic where it can, numeric where it must

Elementary and special functions thread through intervals too, and here the exact /
rigorous split is on full display. `Sin` over `[2, 7]` reaches its minimum `-1` at
`3π/2 ≈ 4.71`, which lies inside the range — so the lower bound is the *exact* `-1`,
while the upper bound stays the closed-form `Sin[2]`:

```mathematica
In[1]:= Sin[Interval[{2, 7}]]
Out[1]= Interval[{-1, Sin[2]}]
```

Widen the range far enough to contain a crest as well as a trough and both bounds
snap to the exact `±1`:

```mathematica
In[2]:= Sin[Interval[{2, 10}]]
Out[2]= Interval[{-1, 1}]
```

Give it inexact endpoints and the bounds come out numeric, rounded outward:

```mathematica
In[3]:= Sin[Interval[{2.5, 5.5}]]
Out[3]= Interval[{-1, 0.598472}]
```

Monotone functions keep their endpoints symbolic where possible — `Exp[1]` is `E`,
`Sqrt[9]` is `3`:

```mathematica
In[4]:= Cos[Interval[{0, 1}]]
Out[4]= Interval[{Cos[1], 1}]

In[5]:= Exp[Interval[{0, 1}]]
Out[5]= Interval[{1, E}]

In[6]:= Sqrt[Interval[{4, 9}]]
Out[6]= Interval[{2, 3}]
```

Functions with **poles** or **discontinuities** produce unions automatically.
`Sec = 1/Cos` blows up where `Cos` crosses zero, so `Sec` over an interval spanning
`π/2` splits into two unbounded pieces:

```mathematica
In[7]:= Sec[Interval[{1, 2}]]
Out[7]= Interval[{-Infinity, Sec[2]}, {Sec[1], Infinity}]
```

`ArcCot` jumps from `-π/2` to `π/2` at the origin, and threading it over a range that
crosses zero reproduces that jump as a gap:

```mathematica
In[8]:= ArcCot[Interval[{-1, 1}]]
Out[8]= Interval[{-1/2 Pi, -1/4 Pi}, {1/4 Pi, 1/2 Pi}]
```

Special functions thread on the sub-domain where each is monotone. `Gamma` is
increasing above its minimum near `1.4616`, and `Γ(2) = 1`, `Γ(3) = 2`:

```mathematica
In[9]:= Gamma[Interval[{2, 3}]]
Out[9]= Interval[{1, 2}]
```

## Membership, unions, and comparisons

`IntervalMemberQ` tests containment, and — because comparisons on intervals are
decided exactly — it is trustworthy right down to the last bit:

```mathematica
In[1]:= IntervalMemberQ[Interval[{0, 2}], 1]
Out[1]= True

In[2]:= IntervalMemberQ[Interval[{0, 2}], Pi]
Out[2]= False
```

`IntervalUnion` and `IntervalIntersection` are the set operations on intervals:

```mathematica
In[3]:= IntervalUnion[Interval[{1, 3}], Interval[{5, 7}]]
Out[3]= Interval[{1, 3}, {5, 7}]

In[4]:= IntervalIntersection[Interval[{1, 5}], Interval[{3, 8}]]
Out[4]= Interval[{3, 5}]
```

Relational operators return a definite `True` or `False` whenever the interval is
**disjoint** from the other operand, and stay symbolic when they overlap — so a
comparison you can settle, settles, and one you cannot is left honest:

```mathematica
In[5]:= Interval[{5, 8}] > Pi
Out[5]= True

In[6]:= Interval[{1, 4}] > Pi
Out[6]= Interval[{1, 4}] > Pi
```

## Rigorous enclosures — the payoff

Now the reason interval arithmetic exists. Return to the opening example. In plain
floating point, the sine of the machine value of π is a small non-zero number, and
nothing warns you it should be zero:

```mathematica
In[1]:= Sin[N[Pi]]
Out[1]= 1.22465e-16
```

`N[Pi]` is π only to machine precision — accurate to about `10^-16`. Feed `Sin` an
interval that **admits that uncertainty**, and the answer rigorously encloses the
truth, `0`:

```mathematica
In[2]:= Sin[Interval[{N[Pi] - 10^-15, N[Pi] + 10^-15}]]
Out[2]= Interval[{-7.65714e-16, 1.01064e-15}]
```

The enclosure straddles zero. The same discipline unmasks the textbook rounding
error: `0.1 + 0.2 - 0.3` is not zero in binary floating point —

```mathematica
In[3]:= 0.1 + 0.2 - 0.3
Out[3]= 5.55112e-17
```

— but the interval computation, rounding each step outward, produces a range that is
**certified to contain** the true value `0`:

```mathematica
In[4]:= Interval[0.1] + Interval[0.2] - Interval[0.3]
Out[4]= Interval[{-1.11022e-16, 1.66533e-16}]

In[5]:= IntervalMemberQ[Interval[0.1] + Interval[0.2] - Interval[0.3], 0]
Out[5]= True
```

**Proving a root exists.** `Sin` is decreasing across `[3, 7/2]`, so its interval
image is the exact range there; the image straddling zero means `Sin` genuinely takes
the value `0` somewhere in that range (it does, at π):

```mathematica
In[6]:= IntervalMemberQ[Sin[Interval[{3, 7/2}]], 0]
Out[6]= True
```

**Sensitive dependence on initial conditions.** The logistic map `x → 4x(1-x)` at
`r = 4` is chaotic. Starting from the point `0.5`, iterate it as an interval and watch
the enclosure widen — each step roughly quadruples the width (the Lyapunov exponent is
`ln 4`), making the loss of information to rounding visible and quantified:

```mathematica
In[7]:= NestList[4 # (1 - #) &, Interval[0.5], 5]
Out[7]= {Interval[{0.5, 0.5}], Interval[{1.0, 1.0}], Interval[{-3.55271e-15, 2.66454e-15}], Interval[{-1.42109e-14, 1.06581e-14}], Interval[{-5.68434e-14, 4.26326e-14}], Interval[{-2.27374e-13, 1.7053e-13}]}
```

Crucially, every one of those widening intervals still **contains** the true orbit
value — the enclosure is loosening, but it never lies:

```mathematica
In[8]:= IntervalMemberQ[Nest[4 # (1 - #) &, Interval[0.5], 5], 0]
Out[8]= True
```

## How it works, and where it stops

Under the hood, an `Interval` is a canonical union of endpoint pairs. Exact endpoints
are computed by the ordinary evaluator, so they stay exact and symbolic; an inexact
(machine or arbitrary-precision) endpoint is nudged one ULP **outward** — lower bounds
down, upper bounds up — so the enclosure is always rigorous. Every result you have
seen encloses the truth, and that is checked continuously by a randomised
containment stress test (`make check-interval`).

Two honest limitations are worth stating. First, Mathilda treats a machine number as an
*exact point*: `Interval[N[Pi]]` is the point `{N[Pi]}`, not a band around true π — so,
as in the opening example, *you* supply the uncertainty you want propagated. Second, a
special function whose interval straddles an extremum or a pole (`Gamma` across its
minimum, say) is left symbolic rather than risk a bound it cannot justify. In both
cases the rule is the same one that makes interval arithmetic worth having: never
print a range that might not contain the answer.
