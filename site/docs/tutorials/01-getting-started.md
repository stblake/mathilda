# Getting started

Welcome to Mathilda, a tiny Mathematica-like computer algebra system written in
C. This first tutorial gets you from a fresh checkout to a working REPL session,
and walks through the handful of things you need to know to start computing:
the `In[]`/`Out[]` loop, the surface syntax, getting help, and assigning values.

Every transcript below was produced by the actual Mathilda binary — type the
`In[...]` lines yourself and you'll see the same `Out[...]`.

## Building Mathilda

Mathilda builds with a plain makefile. From the project root:

```bash
make -j
```

This compiles the C99 sources into a single executable called `Mathilda` in the
project root. The build links three libraries:

- **GMP** — arbitrary-precision integer arithmetic, so big numbers never
  overflow.
- **MPFR** — arbitrary-precision floating point for high-precision numerics.
- **GNU Readline** — line editing, history (up/down arrows), and multiline
  input at the prompt.

Once the build finishes, launch the interactive REPL:

```bash
./Mathilda
```

You'll be greeted with a short banner and the first input prompt:

```
Mathilda - A tiny, LLM-generated, Mathematica-like computer algebra system.

This program is free, open source software and comes with ABSOLUTELY NO WARRANTY.

End a line with '\' to enter a multiline expression. Press Return to evaluate.
Exit by evaluating Quit[] or CONTROL-C.

In[1]:=
```

## The In[]/Out[] loop

Mathilda works just like a calculator with a memory. You type an expression at
the `In[n]:=` prompt, press Return, and Mathilda evaluates it and prints the
result as `Out[n]=`. The counter `n` increments with each evaluation, and you
can refer back to earlier inputs and outputs later in the session.

Let's start with some arithmetic. Notice that Mathilda is *exact* by default —
it never rounds integers or fractions to floating point unless you ask it to:

```mathematica
In[1]:= 2 + 3
Out[1]= 5

In[2]:= 1/3 + 1/6
Out[2]= 1/2
```

The fraction `1/3 + 1/6` stays an exact rational and simplifies to `1/2`, rather
than turning into `0.5`. Big numbers are just as effortless, thanks to GMP.
Here's `2^100` and the factorial of 100, computed exactly:

```mathematica
In[1]:= 2^100
Out[1]= 1267650600228229401496703205376

In[2]:= 100!
Out[2]= 93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000
```

No overflow, no scientific-notation approximation — every digit is correct.

## Surface syntax essentials

Mathilda's syntax mirrors the Wolfram Language. A few conventions cover most of
what you'll type day to day.

**Function calls use square brackets.** You write `f[x]`, not `f(x)`. Round
parentheses are only for grouping arithmetic. For example, `Sin` takes its
argument in square brackets, and rational multiples of `Pi` reduce to exact
values:

```mathematica
In[1]:= Sin[Pi/6]
Out[1]= 1/2
```

**Lists use curly braces.** A list of items is written `{a, b, c}`:

```mathematica
In[1]:= {1, 2, 3}
Out[1]= {1, 2, 3}
```

**The usual operators work as you'd expect** — `+` `-` `*` `/` for arithmetic
and `^` for powers. A couple of operators are worth calling out:

- `==` builds an *equation* (it does not test equality and collapse to a
  boolean unless the two sides are concretely equal). For instance, `x^2 == 4`
  stays symbolic:

```mathematica
In[1]:= x^2 == 4
Out[1]= x^2 == 4
```

- `->` builds a *rule*, used to describe a replacement. Combined with the
  replacement operator `/.` ("replace all"), a rule rewrites parts of an
  expression. Read `expr /. a -> b` as "in `expr`, replace `a` with `b`":

```mathematica
In[1]:= {a, b, c} /. b -> 99
Out[1]= {a, 99, c}
```

Many useful operations are ordinary functions. `Range[n]` builds the list
`{1, ..., n}`, and you can give it a step:

```mathematica
In[1]:= Range[5]
Out[1]= {1, 2, 3, 4, 5}

In[2]:= Range[1, 10, 2]
Out[2]= {1, 3, 5, 7, 9}
```

`Map[f, list]` applies a function to every element of a list — here `f` stays
symbolic because we never defined it, which nicely shows what `Map` does:

```mathematica
In[1]:= Map[f, {a, b, c}]
Out[1]= {f[a], f[b], f[c]}
```

And `Factor` factors a polynomial:

```mathematica
In[1]:= Factor[x^2 - 1]
Out[1]= (-1 + x) (1 + x)
```

## Getting help

Every built-in function carries a short documentation string. To read it, type
`?` followed by the name. For example, `?Range`:

```mathematica
In[1]:= ?Range
Out[1]= "Range[n]
	generates the list {1, 2, 3, ..., n}.
Range[n, m]
	generates the list {n, n + 1, ..., m - 1, m}.
Range[n, m, d]
	uses step d."
```

The docstrings are concise but complete. Here's `?Sin`:

```mathematica
In[1]:= ?Sin
Out[1]= "Sin[z]
	gives the sine of z (argument in radians).
Sin is Listable. Numeric inputs are evaluated via libm (Real) or MPFR
(arbitrary precision); rational multiples of Pi reduce to exact values."
```

Some functions have longer entries. `?Factor`, for instance, describes its
extension options and the underlying algorithm — the output begins:

```mathematica
In[1]:= ?Factor
Out[1]= "Factor[poly] factors a polynomial over the integers.
Factor[poly, Extension -> alpha] factors over Q(alpha), where alpha is
Sqrt[c], c^(1/n) (rational c), or I.  ...
```

(truncated — the full string continues with the compositum and primitive-element
details). When you're exploring, `?Name` is the fastest way to recall what a
function does and how to call it.

## Assignment basics

Use `=` (the `Set` operator) to give a name a value. Once assigned, the name
stands in for that value everywhere in the session:

```mathematica
In[1]:= x = 5
Out[1]= 5

In[2]:= x^2 + 1
Out[2]= 26
```

Names persist for the rest of the session, so it's good practice to release a
name when you're done with it. `Clear[x]` removes the value, and `x` becomes a
plain symbol again:

```mathematica
In[1]:= x = 5
Out[1]= 5

In[2]:= Clear[x]
Out[2]= Null

In[3]:= x
Out[3]= x
```

`Clear` returns `Null` (the "no interesting value" result), and afterwards `x`
evaluates to itself — it's symbolic once more.

## Multiline input and exiting

For a long expression, end a line with a backslash `\` and Mathilda will wait
for the rest on the next line before evaluating. The two physical lines below
are treated as one input:

```mathematica
In[1]:= 1 + 2 + \
3 + 4
Out[1]= 10
```

When you're finished, leave the REPL by evaluating `Quit[]`:

```mathematica
In[1]:= Quit[]
```

You can also press **Ctrl-C** or send EOF (**Ctrl-D** on an empty line) to exit.

## Where to next

You now know enough to drive the REPL: evaluate expressions, build lists, apply
functions, look up help, and assign names. The next tutorial digs into what's
really going on under the hood.

- **[2. Expressions & evaluation](02-expressions-and-evaluation.md)** —
  everything is an expression; meet `FullForm`, `Head`, the attribute system,
  and the fixed-point evaluator.
- **[Function documentation](../documentation/index.md)** — the full reference
  for every built-in, organized by category.
