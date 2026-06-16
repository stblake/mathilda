# Arithmetic

This tutorial is a hands-on tour of Mathilda's arithmetic — the bedrock that
every other part of the system is built on. You will add, multiply, and divide
exact numbers that never overflow and never round; take roots; move deliberately
between exact, machine-precision, and arbitrary-precision values; round numbers
and pull out their parts; work with complex numbers; compute factorials,
binomials, and Fibonacci numbers; do arithmetic over whole lists; and finally
take numbers apart digit by digit in any radix. The closing examples put all of
this to work on a couple of genuinely non-trivial problems.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself (without the prompt) and you will see the same
`Out[...]`. Mathilda's output is sometimes ordered differently from a textbook —
it sorts terms into a canonical order, so constants and lower-degree terms tend
to come first — but it is always mathematically correct.

## The basic operators

Arithmetic uses the operators you expect — `+`, `-`, `*`, `/`, `^` — with the
usual precedence, so multiplication binds tighter than addition. Each operator is
just sugar for a named function (`Plus`, `Subtract`, `Times`, `Divide`,
`Power`), and you may call either form:

```mathematica
In[1]:= 2 + 3 * 4
Out[1]= 14

In[2]:= Subtract[10, 3]
Out[2]= 7

In[3]:= Divide[12, 4]
Out[3]= 3

In[4]:= 12/8
Out[4]= 3/2
```

`In[1]` evaluates the multiplication first, then the addition, giving `14`.
`In[2]` and `In[3]` show the function spellings of `-` and `/`. The last line is
the first surprise for newcomers: `12/8` does **not** become `1.5`. Division of
exact integers stays exact, and Mathilda reduces `12/8` to the fraction `3/2` in
lowest terms — a theme we return to throughout.

## Exact integers never overflow

In most programming languages an integer is a fixed-size machine word — 64 bits —
and arithmetic silently wraps around once a result grows too large. Mathilda has
no such ceiling. Integers are **arbitrary precision**: they grow to whatever size
the answer needs, using the GMP big-integer library under the hood.

```mathematica
In[1]:= 2^100
Out[1]= 1267650600228229401496703205376

In[2]:= 20!
Out[2]= 2432902008176640000

In[3]:= Head[2^100]
Out[3]= Integer
```

`2^100` is a 31-digit number and `20!` already exceeds what a 64-bit integer can
hold — both are computed exactly, every digit correct. Yet `In[3]` shows the head
is still plain `Integer`: there is no special "big integer" type to opt into.
Whether a value fits in a machine word or spans hundreds of digits is an
implementation detail you never have to think about; Mathilda promotes its
internal storage automatically.

## Exact rationals stay fractions

Divide two integers and Mathilda keeps the result as an exact **rational**,
reduced to lowest terms, rather than collapsing it to a decimal:

```mathematica
In[1]:= 1/3 + 1/6
Out[1]= 1/2

In[2]:= 3/6 + 2/4
Out[2]= 1

In[3]:= Rational[6, 4]
Out[3]= 3/2

In[4]:= Precision[1/3]
Out[4]= Infinity
```

`1/3 + 1/6` is added over a common denominator and reduced to `1/2`; `3/6 + 2/4`
is `1/2 + 1/2 = 1`, an integer, so Mathilda simplifies it all the way. `Rational`
is the explicit head behind the `/` notation, and it too auto-reduces (`6/4`
becomes `3/2`). `In[4]` asks how much precision an exact number carries:
`Precision` reports the number of reliable significant digits, and for an exact
quantity the answer is `Infinity` — there is no rounding error, so *every* digit
is reliable.

## Square roots

`Sqrt` is exact too. It simplifies perfect squares, pulls square factors out from
under the radical, and leaves genuinely irrational roots in symbolic form rather
than approximating them:

```mathematica
In[1]:= Sqrt[16]
Out[1]= 4

In[2]:= Sqrt[2]
Out[2]= Sqrt[2]

In[3]:= Sqrt[8]
Out[3]= 2 Sqrt[2]

In[4]:= Sqrt[-4]
Out[4]= 2*I
```

`Sqrt[16]` is the perfect square `4`. `Sqrt[2]` has no exact value, so Mathilda
keeps it symbolic. `Sqrt[8]` factors as `Sqrt[4 * 2]` and the perfect square
comes out front, leaving `2 Sqrt[2]`. And `Sqrt[-4]` is no error — Mathilda
works happily over the complex numbers and returns `2*I`, where `I` is the
imaginary unit. (We meet complex numbers properly below.)

## Machine vs arbitrary precision

Write a number with a decimal point and you opt into a different world. `1.0`,
`3.5`, `0.2` are **machine-precision reals** — IEEE double-precision floating
point, the same fast hardware arithmetic your CPU does natively. A single decimal
point anywhere is contagious: it turns the whole calculation approximate.

```mathematica
In[1]:= 1.0 / 3.0
Out[1]= 0.333333

In[2]:= Precision[1.0]
Out[2]= MachinePrecision

In[3]:= $MachinePrecision
Out[3]= 15.9546
```

Compare `In[1]` with the exact `1/3` above: the decimal points forced a
floating-point division, and Mathilda printed a rounded decimal. `Precision`
reports the special value `MachinePrecision` rather than a number, because the
width is fixed by the hardware at about 16 decimal digits (`$MachinePrecision` is
that count). Mathilda displays only the first six significant digits by default;
the full ≈16 digits are still stored internally.

Those ~16 digits are *all* you get, and rounding error accumulates. The classic
demonstration:

```mathematica
In[1]:= 0.1 + 0.2 - 0.3
Out[1]= 2.77556e-17

In[2]:= 1/10 + 2/10 - 3/10
Out[2]= 0
```

The answer should be exactly zero. With machine reals it is not (`In[1]`),
because `0.1`, `0.2`, and `0.3` cannot be represented exactly in binary — each is
rounded as it is entered, and the tiny errors survive. Do the same sum *exactly*
with rationals (`In[2]`) and the error vanishes completely.

You rarely type decimal points yourself: you keep things exact and ask for a
numeric value only at the end with `N`. `N[expr]` gives a machine-precision real,
while `N[expr, d]` computes `d` significant digits using arbitrary-precision
arithmetic (the MPFR library under the hood):

```mathematica
In[1]:= N[Pi]
Out[1]= 3.14159

In[2]:= N[1/7]
Out[2]= 0.142857

In[3]:= N[Pi, 50]
Out[3]= 3.1415926535897932384626433832795028841971693993751

In[4]:= N[Sqrt[2], 40]
Out[4]= 1.4142135623730950488016887242096980785697

In[5]:= N[E, 60]
Out[5]= 2.718281828459045235360287471352662497757247093699959574966968
```

`Pi` on its own stays the exact symbol `Pi`; wrapping it in `N` produces its
decimal value (`In[1]`). Fifty digits of π, forty of √2, sixty of *e* — all
correct, and you may ask for a thousand just as easily. `Precision` and
`Accuracy` confirm the extra digits are genuine. `Precision` counts reliable
significant digits; `Accuracy` counts reliable digits to the right of the decimal
point:

```mathematica
In[1]:= Precision[N[Pi, 50]]
Out[1]= 50.272

In[2]:= Accuracy[N[Pi, 50]]
Out[2]= 49.7749

In[3]:= N[2^100]
Out[3]= 1.26765e+30
```

The result of `N[Pi, 50]` carries about 50 reliable significant digits (the few
extra are internal guard digits) and, since π ≈ 3, about 49.8 of those sit after
the decimal point. `In[3]` is the cautionary note about going *down* the ladder:
applying bare `N` to a large exact integer throws away all but the first ~16 of
its 31 digits.

When you need to fix the precision or accuracy of a number directly — rather than
recompute a symbolic expression — use `SetPrecision` and `SetAccuracy`:

```mathematica
In[1]:= SetPrecision[Pi, 40]
Out[1]= 3.1415926535897932384626433832795028841971

In[2]:= SetAccuracy[1.23456, 3]
Out[2]= 1.234
```

`SetPrecision[Pi, 40]` re-evaluates π to 40-digit precision, like `N[Pi, 40]`.
`SetAccuracy[1.23456, 3]` keeps three digits past the decimal point. The takeaway
that runs through this section: **keep numbers exact for as long as you can, and
approximate only at the last step — and even then, ask for as many digits as you
actually need.**

## Rounding and number parts

Four functions turn a real into a nearby integer. `Floor` rounds down,
`Ceiling` rounds up, and `Round` rounds to the nearest integer:

```mathematica
In[1]:= Floor[3.7]
Out[1]= 3

In[2]:= Ceiling[3.2]
Out[2]= 4

In[3]:= Round[3.5]
Out[3]= 4

In[4]:= Round[2.5]
Out[4]= 2

In[5]:= Floor[-3.2]
Out[5]= -4
```

`Floor` and `Ceiling` are the integer below and above. `Round` goes to the
nearest integer, and on a halfway value it rounds to the nearest **even**
integer — that is why `Round[3.5]` is `4` but `Round[2.5]` is `2`. This
"round half to even" rule (banker's rounding) avoids the upward bias of always
rounding `.5` up. `Floor[-3.2]` is `-4` because that is the integer *below*
`-3.2`. All three work on exact rationals too, so `Floor[7/2]` is `3`.

`Sign` extracts the sign as `-1`, `0`, or `1`, and `Abs` gives the absolute value
(which, for a complex number, is its modulus):

```mathematica
In[1]:= Sign[-5]
Out[1]= -1

In[2]:= Sign[0]
Out[2]= 0

In[3]:= Abs[-7]
Out[3]= 7

In[4]:= Abs[3 - 4 I]
Out[4]= 5
```

`Sign` reports the direction of a number on the line. `Abs[-7]` is its distance
from zero. The last line is the gateway to the next section: `Abs[3 - 4 I]` is the
modulus √(3² + 4²) = `5`, the length of the complex number as a vector.

## Complex numbers

A complex number `a + b I` has a real part, an imaginary part, an argument
(angle), and a conjugate. `Re` and `Im` pick out the two parts, and `ReIm`
returns both at once as a list:

```mathematica
In[1]:= Re[3 + 4 I]
Out[1]= 3

In[2]:= Im[3 + 4 I]
Out[2]= 4

In[3]:= ReIm[3 + 4 I]
Out[3]= {3, 4}

In[4]:= Conjugate[3 + 4 I]
Out[4]= 3 - 4*I
```

`Re` and `Im` read off the rectangular coordinates; `ReIm` is the convenient
pairing of the two. `Conjugate` flips the sign of the imaginary part, reflecting
the number across the real axis. `Arg` gives the polar angle (the argument),
measured in radians and returned exactly when possible:

```mathematica
In[1]:= Arg[I]
Out[1]= 1/2 Pi

In[2]:= Arg[-1]
Out[2]= Pi

In[3]:= Arg[1 + I]
Out[3]= 1/4 Pi
```

`I` points straight up, so its argument is `π/2`; `-1` points along the negative
real axis at angle `π`; and `1 + I` lies on the 45° diagonal, argument `π/4`. You
can also build a complex number from its parts with `Complex`, and ordinary
arithmetic combines them:

```mathematica
In[1]:= Complex[2, 3]
Out[1]= 2 + 3*I

In[2]:= (2 + 3 I) (1 - I)
Out[2]= 5 + I
```

`Complex[2, 3]` is the literal `2 + 3 I`. `In[2]` multiplies two complex numbers
the way you would by hand — `(2 + 3I)(1 - I) = 2 - 2I + 3I - 3I² = 5 + I`, using
`I² = -1`.

## Combinatorial functions

Mathilda has the standard counting functions built in, all exact and all
overflow-proof. `Factorial` (also written `n!`) and its variants come first:

```mathematica
In[1]:= Factorial[5]
Out[1]= 120

In[2]:= Factorial2[6]
Out[2]= 48

In[3]:= Factorial2[7]
Out[3]= 105

In[4]:= FactorialPower[10, 3]
Out[4]= 720
```

`Factorial[5]` is `5! = 120`. `Factorial2` is the *double* factorial, the product
of every second integer down to 1 or 2: `Factorial2[6] = 6·4·2 = 48` and
`Factorial2[7] = 7·5·3·1 = 105`. `FactorialPower[10, 3]` is the falling factorial
`10·9·8 = 720`, the number of ordered ways to pick 3 items from 10.

`Binomial` counts unordered selections, and `Fibonacci`/`LucasL` generate the two
classic integer sequences:

```mathematica
In[1]:= Binomial[10, 3]
Out[1]= 120

In[2]:= Fibonacci[10]
Out[2]= 55

In[3]:= Fibonacci[100]
Out[3]= 354224848179261915075

In[4]:= LucasL[10]
Out[4]= 123
```

`Binomial[10, 3]` is "10 choose 3" = `120`, the number of 3-element subsets of a
10-element set. `Fibonacci[10]` is `55`; `Fibonacci[100]` is a 21-digit number,
computed exactly and instantly. `LucasL` produces the Lucas numbers, the
Fibonacci sequence's companion that starts `2, 1, 3, 4, 7, ...` — so
`LucasL[10] = 123`.

## List arithmetic

Arithmetic does not stop at single numbers. Four functions operate on whole
lists at once. `Total` sums a list; `Accumulate` returns its running (cumulative)
sums:

```mathematica
In[1]:= Total[{1, 2, 3, 4, 5}]
Out[1]= 15

In[2]:= Accumulate[{1, 2, 3, 4}]
Out[2]= {1, 3, 6, 10}
```

`Total` collapses the list to its sum, `15`. `Accumulate` keeps every partial sum
along the way: `1`, `1+2`, `1+2+3`, `1+2+3+4`. The two are related — the last
element of `Accumulate` is always `Total`.

`Differences` is the discrete derivative — each element minus its predecessor —
and `Ratios` is its multiplicative cousin, each element divided by its
predecessor:

```mathematica
In[1]:= Differences[{1, 4, 9, 16, 25}]
Out[1]= {3, 5, 7, 9}

In[2]:= Ratios[{1, 2, 4, 8, 16}]
Out[2]= {2, 2, 2, 2}
```

The input to `In[1]` is the list of squares; their first differences are the
consecutive odd numbers `3, 5, 7, 9` — the classic fact that the gaps between
squares grow by 2 each time. `Ratios` of the powers of two is the constant list
`2, 2, 2, 2`, exposing the geometric ratio. `Differences` even takes a second
argument for higher-order differences: `Differences[list, 2]` differences twice,
and for the squares that constant `2` falls straight out.

## Digits and radix

The final family takes numbers apart digit by digit, in any base. `IntegerDigits`
returns the digits of an integer (base 10 by default, or a base you name), and
`FromDigits` reassembles them:

```mathematica
In[1]:= IntegerDigits[12345]
Out[1]= {1, 2, 3, 4, 5}

In[2]:= IntegerDigits[255, 2]
Out[2]= {1, 1, 1, 1, 1, 1, 1, 1}

In[3]:= FromDigits[{1, 1, 1, 1}, 2]
Out[3]= 15

In[4]:= IntegerString[255, 16]
Out[4]= "ff"
```

`IntegerDigits[12345]` lists the decimal digits in order. In base 2, `255` is
eight ones (`In[2]`). `FromDigits` is the inverse: the binary digits `1111` make
the number `15`. When you want the digits as a *string* rather than a list — handy
for hexadecimal — `IntegerString` does the job, rendering `255` as `"ff"`.

Three companions summarise an integer without listing every digit.
`IntegerLength` counts digits, `IntegerExponent[n, b]` counts how many times the
base `b` divides `n`, and `DigitCount` tallies how often each digit appears:

```mathematica
In[1]:= IntegerLength[2^100]
Out[1]= 31

In[2]:= IntegerExponent[1000, 10]
Out[2]= 3

In[3]:= DigitCount[1122334455]
Out[3]= {2, 2, 2, 2, 2, 0, 0, 0, 0, 0}
```

`IntegerLength[2^100]` confirms the 31 digits we saw earlier. `IntegerExponent`
finds the three factors of 10 in `1000` (its three trailing zeros). `DigitCount`
returns a length-10 tally — the count of `1`s, `2`s, …, `9`s, and finally `0`s —
so `1122334455` shows two each of the digits 1 through 5 and none of the rest.

Reals have their own decomposition. `RealDigits` gives the significant digits plus
the position of the decimal point; `RealExponent` is the (real-valued) base-10
exponent; and `MantissaExponent` splits a number into a mantissa in `[1/10, 1)`
and an integer power of 10:

```mathematica
In[1]:= RealDigits[Pi, 10, 10]
Out[1]= {{3, 1, 4, 1, 5, 9, 2, 6, 5, 3}, 1}

In[2]:= RealExponent[12345]
Out[2]= 4.09149

In[3]:= MantissaExponent[3.14]
Out[3]= {0.314, 1}
```

`RealDigits[Pi, 10, 10]` returns the first ten decimal digits of π together with
`1`, meaning the decimal point sits after the first digit (π ≈ 3.14…, so it is
`3 × 10^0`, exponent index 1). `RealExponent[12345]` is `log₁₀(12345) ≈ 4.09`.
`MantissaExponent[3.14]` writes `3.14` as `0.314 × 10^1`.

## Application: Catalan numbers from binomials

The Catalan numbers `Cₙ = C(2n, n)/(n + 1)` count an astonishing variety of
things — balanced bracketings, triangulations of a polygon, binary trees — and
they fall straight out of `Binomial` and exact division. Generate the first nine
of them with `Table`:

```mathematica
In[1]:= Table[Binomial[2 n, n]/(n + 1), {n, 0, 8}]
Out[1]= {1, 1, 2, 5, 14, 42, 132, 429, 1430}
```

Every division `Binomial[2n, n]/(n+1)` comes out a whole number — that it always
does is a small theorem, and Mathilda's exact rational arithmetic both performs
the division and verifies the divisibility for free. The sequence
`1, 1, 2, 5, 14, 42, …` is the Catalan numbers. Their running sums are themselves
meaningful, and `Accumulate` produces them in one step:

```mathematica
In[1]:= Accumulate[Table[Binomial[2 n, n]/(n + 1), {n, 0, 8}]]
Out[1]= {1, 2, 4, 9, 23, 65, 197, 626, 2056}
```

This chains three arithmetic tools — `Binomial`, exact `/`, and `Accumulate` —
into a single pipeline, with no floating point anywhere and every digit exact.

## Application: dissecting 100!

`100!` is a 158-digit integer. Without ever converting it to an approximate
float, we can interrogate its structure entirely with the digit tools from above:

```mathematica
In[1]:= IntegerLength[100!]
Out[1]= 158

In[2]:= IntegerExponent[100!, 10]
Out[2]= 24

In[3]:= Total[DigitCount[100!]]
Out[3]= 158

In[4]:= Total[IntegerDigits[100!]]
Out[4]= 648
```

`IntegerLength` confirms the 158 digits. `IntegerExponent[100!, 10]` is `24`: the
factorial ends in 24 zeros, one for each factor of 10, governed by the count of
5s among `1..100` (since 2s are plentiful) — `20 + 4 = 24`. As a sanity check,
`Total[DigitCount[100!]]` sums the per-digit tallies back to `158`, matching the
length. And `Total[IntegerDigits[100!]]` adds up all 158 digits to get `648` — a
computation that would be tedious by hand but is exact and immediate here.

## Application: a Fibonacci identity, verified

A famous identity, *Cassini's identity*, states that
`Fₙ² − Fₙ₋₁·Fₙ₊₁ = (−1)ⁿ` for every `n`. We can verify it across a range of `n`
in one line, using only multiplication, subtraction, and `Fibonacci`:

```mathematica
In[1]:= Table[Fibonacci[n]^2 - Fibonacci[n - 1] Fibonacci[n + 1], {n, 1, 8}]
Out[1]= {1, -1, 1, -1, 1, -1, 1, -1}
```

The result alternates `1, -1, 1, -1, …` — precisely `(−1)ⁿ` — confirming the
identity for `n = 1` through `8` with exact integer arithmetic, so the check is a
proof of those cases, not a numerical approximation. A second identity,
`Fₙ² + Fₙ₊₁² = F₂ₙ₊₁`, says the sum of two consecutive Fibonacci squares is again
a Fibonacci number:

```mathematica
In[1]:= Fibonacci[10]^2 + Fibonacci[11]^2 - Fibonacci[21]
Out[1]= 0
```

The difference is exactly `0`, so `F₁₀² + F₁₁² = F₂₁` holds on the nose. As a
final flourish, the ratio of consecutive Fibonacci numbers converges to the
golden ratio — and exact arithmetic lets us pin that down to many digits:

```mathematica
In[1]:= N[Fibonacci[200]/Fibonacci[199], 25]
Out[1]= 1.6180339887498948482045868

In[2]:= N[GoldenRatio, 25]
Out[2]= 1.6180339887498948482045868
```

`Fibonacci[200]/Fibonacci[199]` is an exact rational of two 40-ish-digit
integers; asking `N` for 25 digits shows it already agrees with the golden ratio
`GoldenRatio` to all 25 — the convergence is that fast. Exact arithmetic followed
by a single deliberate `N` at the end gives the best of both worlds.

## Where to next

You can now drive Mathilda's exact integers and rationals, take roots, move
fluently between machine and arbitrary precision with `N`, `Precision`, and
`SetPrecision`, round numbers and read off their parts, manipulate complex
numbers, compute factorials and Fibonacci numbers, do arithmetic over lists, and
take any number apart digit by digit. This is the numeric foundation the rest of
Mathilda rests on.

- **[5. Number theory](05-number-theory.md)** — push these integers further:
  primes and factorisation, GCD and LCM, modular arithmetic, and the
  number-theoretic functions that build directly on the exact arithmetic above.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions. The
  [arithmetic](../documentation/arithmetic/index.md) section documents every
  function used here, including the complete options for `N`, `Round`,
  `RealDigits`, and the rest.
