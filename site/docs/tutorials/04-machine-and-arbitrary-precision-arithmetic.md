# Machine and arbitrary-precision arithmetic

So far every number you have typed has just worked. Behind that simplicity
Mathilda is constantly making a choice: should this number be kept *exactly*, or
approximated as a *decimal*? And if it is approximated, with how many digits?
This tutorial pulls back the curtain on Mathilda's three kinds of numbers —
**exact** integers and rationals, fast **machine-precision** reals, and
**arbitrary-precision** reals — so you know which one you are holding and how to
move between them.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself (without the prompt) and you will see the same
`Out[...]`.

## Exact integers never overflow

In most programming languages an integer is a fixed-size machine word — 64 bits —
and arithmetic silently wraps around once a result grows too large. Mathilda has
no such ceiling. Integers are **arbitrary precision**: they grow to whatever size
the answer needs, using the GMP big-integer library under the hood.

```mathematica
In[1]:= 20!
Out[1]= 2432902008176640000

In[2]:= 2^100
Out[2]= 1267650600228229401496703205376
```

`20!` already exceeds what a 64-bit integer can hold, and `2^100` is a 31-digit
number — both are computed exactly, every digit correct. There is no special
"big integer" type to opt into; an integer is just an integer, and Mathilda
promotes its internal storage automatically when a value gets large:

```mathematica
In[1]:= Head[2^100]
Out[1]= Integer
```

The head is still plain `Integer`. Whether a value fits in a machine word or
spans hundreds of digits is an implementation detail you never have to think
about.

## Exact rationals stay fractions

Divide two integers and Mathilda keeps the result as an exact **rational**,
reduced to lowest terms, rather than collapsing it to a decimal:

```mathematica
In[1]:= 1/3 + 1/6
Out[1]= 1/2

In[2]:= 3/6 + 2/4
Out[2]= 1
```

`1/3 + 1/6` is added over a common denominator and reduced to `1/2`; `3/6 + 2/4`
is `1/2 + 1/2 = 1`, an integer, so Mathilda simplifies it all the way. No
rounding happens at any step. We can ask Mathilda how much precision an exact
number carries:

```mathematica
In[1]:= Precision[1/3]
Out[1]= Infinity
```

`Precision` reports the number of reliable significant digits. For an exact
quantity the answer is `Infinity` — there is no rounding error, so *every* digit
is reliable. Exact integers and rationals are Mathilda's default and its most
trustworthy numbers; reach for anything else only when you have a reason to.

## Machine-precision reals: fast but approximate

Write a number with a decimal point and you opt into a different world. `1.0`,
`3.5`, `0.2` are **machine-precision reals** — IEEE double-precision floating
point, the same fast hardware arithmetic your CPU does natively. A single decimal
point anywhere is contagious: it turns the whole calculation approximate.

```mathematica
In[1]:= 1.0 / 3.0
Out[1]= 0.333333

In[2]:= Precision[1.0]
Out[2]= MachinePrecision
```

Compare `In[1]` here with the exact `1/3` above: the decimal points forced a
floating-point division, and Mathilda printed a rounded decimal. `Precision`
reports the special value `MachinePrecision` rather than a number, because the
width is fixed by the hardware. That width is about 16 decimal digits:

```mathematica
In[1]:= $MachinePrecision
Out[1]= 15.9546
```

(Mathilda displays only the first six significant digits by default; the full
≈16 digits are still stored internally.)

Machine reals are fast, but those ~16 digits are *all* you get, and rounding
error accumulates. The classic demonstration:

```mathematica
In[1]:= 0.1 + 0.2 - 0.3
Out[1]= 2.77556e-17
```

The answer should be exactly zero. It is not, because `0.1`, `0.2`, and `0.3`
cannot be represented exactly in binary floating point — each is rounded as it is
entered, and the tiny errors survive into the result. Do the same sum *exactly*
and the error vanishes:

```mathematica
In[1]:= 1/10 + 2/10 - 3/10
Out[1]= 0
```

This is the central trade-off: machine reals are quick and fine for most
numerical work, but they are approximate from the moment you type the decimal
point. When the correctness of every digit matters, stay exact.

## N: from exact to approximate, on demand

You rarely type decimal points yourself for symbolic constants — you keep things
exact and ask for a numeric value only at the end. That is what `N` is for. `N[expr]`
evaluates `expr` to a machine-precision real:

```mathematica
In[1]:= N[Pi]
Out[1]= 3.14159

In[2]:= N[1/7]
Out[2]= 0.142857
```

`Pi` on its own stays the exact symbol `Pi` (try it); wrapping it in `N` produces
its decimal value. Likewise `N[1/7]` turns the exact fraction into a decimal.
This is the usual workflow: compute symbolically for as long as you can, then
`N` the final result.

## Arbitrary-precision reals: as many digits as you ask for

Machine precision caps you at ~16 digits. When you need more, give `N` a second
argument: `N[expr, d]` computes `expr` to `d` significant decimal digits using
arbitrary-precision arithmetic (the MPFR library under the hood).

```mathematica
In[1]:= N[Pi, 50]
Out[1]= 3.1415926535897932384626433832795028841971693993751

In[2]:= N[Sqrt[2], 40]
Out[2]= 1.4142135623730950488016887242096980785697

In[3]:= N[E, 60]
Out[3]= 2.718281828459045235360287471352662497757247093699959574966968
```

Fifty digits of π, forty of √2, sixty of *e* — all correct. Ask for a thousand
and Mathilda will compute a thousand. Unlike a machine real, this is genuine
extra precision, and `Precision` confirms it:

```mathematica
In[1]:= Precision[N[Pi, 50]]
Out[1]= 50.272
```

The result carries about 50 reliable digits (the few extra are internal guard
digits that protect the requested precision from rounding). Contrast that with
`Precision[N[Pi]]`, which would report `MachinePrecision`.

## Choosing the right kind of number

Putting it together, Mathilda gives you a ladder of three numeric kinds:

| Kind | How you get it | Precision | Use when |
|------|----------------|-----------|----------|
| Exact integer / rational | type `2^100`, `1/3` | `Infinity` (every digit exact) | correctness matters; the default |
| Machine real | a decimal point, or `N[expr]` | `MachinePrecision` (~16 digits) | fast numerics, plotting, "good enough" |
| Arbitrary-precision real | `N[expr, d]` | `d` digits, your choice | many correct digits required |

A cautionary note on going *down* the ladder: applying `N` to a large exact
integer throws away its exact value:

```mathematica
In[1]:= N[2^100]
Out[1]= 1.26765e+30
```

That machine real has lost all but the first ~16 of the original 31 digits. The
lesson that runs through this whole tutorial: **keep numbers exact for as long as
you can, and approximate only at the last step — and even then, ask for as many
digits as you actually need.**

## Where to next

You now know the three kinds of numbers Mathilda works with, how to tell them
apart with `Precision` and `Head`, how to drop to machine precision or climb to
arbitrary precision with `N`, and — most importantly — when each is the right
tool.

- **[5. Algebra](05-algebra.md)** — put exact arithmetic to work: expand and
  factor polynomials, reshape rational expressions, simplify, and solve
  equations and systems.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions, including the complete details of `N`, `Precision`,
  `Accuracy`, and the constants `Pi`, `E`, and `$MachinePrecision` used above.
