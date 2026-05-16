# Arithmetic in Mathilda

In the [pattern matching](./Mathilda_patterns_blog.md) and [functional programming](./Mathilda_functional_blog.md) notes I focused on the language layer. This time I would like to look one floor down — at how Mathilda actually does arithmetic. Mathilda supports three numeric regimes:

1. **Machine precision** — IEEE-754 double-precision floats, the same 64-bit numbers your C compiler gives you.
2. **Arbitrary-precision integers** — exact integers of any size, backed by [GNU MP (GMP)](https://gmplib.org/).
3. **Arbitrary-precision reals** — `N[expr, n]` evaluates `expr` to `n` decimal digits using [MPFR](https://www.mpfr.org/) for correctly-rounded floating-point.

There is no other numeric tower. In particular, **Mathilda does not implement Mathematica's significance arithmetic** — the precision of an MPFR number stays where you put it, even through cancellations and ill-conditioned operations. I will say more about this at the end. All of the examples below were run in Mathilda; the `Out[]` lines are the actual REPL output.

## Machine integers and the slide into bignums

Integers begin life as 64-bit signed values. Once an operation would overflow, Mathilda quietly promotes the result to a GMP `mpz_t` and keeps going:

```mathematica
In[1]:= 2^62
Out[1]= 4611686018427387904

In[2]:= 2^63
Out[2]= 9223372036854775808

In[3]:= 2^64
Out[3]= 18446744073709551616

In[4]:= 2^200
Out[4]= 1606938044258990275541962092341162602522202993782792835301376
```

The transition between `EXPR_INTEGER` and `EXPR_BIGINT` is invisible to the user. After every arithmetic operation a helper `expr_bigint_normalize()` checks whether the result fits back into an `int64_t` and demotes when it does, so simple arithmetic stays fast.

The classic stress test is the factorial. `100!` has 158 decimal digits and is exact:

```mathematica
In[5]:= 100!
Out[5]= 93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000

In[6]:= 1 + Floor[N[Log[10, 100!]]]
Out[6]= 158
```

## Exact rationals

When you mix integers with division, Mathilda keeps the result as an exact rational. There is no silent conversion to float:

```mathematica
In[7]:= 1/3 + 1/4
Out[7]= 7/12

In[8]:= Total[1/Range[10]]
Out[8]= 7381/2520
```

This matters in numerical analysis. Floating-point addition is not associative — the canonical example is the failure of `0.1 + 0.2 == 0.3`:

```mathematica
In[9]:= (0.1 + 0.2) - 0.3
Out[9]= 2.77556e-17

In[10]:= 1/10 + 1/5 - 3/10
Out[10]= 0
```

The first line shows the well-known $2^{-55} \approx 2.78 \times 10^{-17}$ rounding residue; the second line is exact because `1/10`, `1/5`, and `3/10` are stored as `Rational[1,10]`, `Rational[1,5]`, `Rational[3,10]`.

## Catastrophic cancellation, the textbook way

A classic pitfall in numerical computing is computing `Sqrt[b^2 + 1] - b` for large `b`. Algebraically the answer is positive and tiny: $1/(b + \sqrt{b^2 + 1}) \approx 1/(2b)$. Numerically it can lose almost all of its significant digits.

```mathematica
In[11]:= b = 100000000000.0
Out[11]= 1e+11

In[12]:= b^2 + 1
Out[12]= 1e+22

In[13]:= Sqrt[b^2 + 1] - b
Out[13]= 0.000152588
```

The exact answer is `4.9999999999999999999998750...e-12`. Machine precision lost every digit: `1` is not even representable next to `10^22` in 53 bits of mantissa, so `b^2 + 1` collapses back to `b^2`, and the subtraction returns the rounding error of $\sqrt{b^2}$.

Switching to exact integers and then asking for `n` decimal digits at the end gives the right answer:

```mathematica
In[14]:= N[Sqrt[100000000000^2 + 1] - 100000000000, 30]
Out[14]= 4.999999999999999999999875000001e-12

In[15]:= N[1/(100000000000 + Sqrt[100000000000^2 + 1]), 30]
Out[15]= 5e-12
```

The same expression rewritten to avoid the cancellation (`1/(b + Sqrt[b^2 + 1])`) is well-conditioned and falls out cleanly at any precision.

## Newton's iteration

Newton's iteration for $\sqrt{2}$ is a `NestList` over rational numbers. Because Mathilda keeps rationals exact, the numerator and denominator double in length on each step:

```mathematica
In[16]:= NestList[(# + 2/#)/2 &, 1, 5]
Out[16]= {1, 3/2, 17/12, 577/408, 665857/470832, 886731088897/627013566048}

In[17]:= N[886731088897/627013566048, 30]
Out[17]= 1.414213562373095048801689623503

In[18]:= N[Sqrt[2], 30]
Out[18]= 1.414213562373095048801688724209
```

Five Newton steps from the seed `1` already agree with the true $\sqrt{2}$ to 21 decimal places — the iteration's quadratic convergence in action. The same iteration with a `1.0` seed converges in machine precision and stops dead at $\approx$ 16 digits.

## Number theory

GMP makes number-theoretic computations almost trivial. Mersenne-prime tests work out to enormous sizes:

```mathematica
In[19]:= PrimeQ[2^61 - 1]
Out[19]= True

In[20]:= PrimeQ[2^89 - 1]
Out[20]= True

In[21]:= PrimeQ[2^127 - 1]
Out[21]= True

In[22]:= PrimeQ[2^257 - 1]
Out[22]= False
```

`2^257 - 1` is composite — Frank Cole's famous 1903 talk produced the factorisation `761838257287 × 193707721 × ...` for the 67-bit Mersenne, and the 257-bit one is also notoriously composite. Mathilda uses GMP's Miller–Rabin for `PrimeQ`, with deterministic small-prime trial division as a precheck.

Fermat numbers are another quick test. $F_5 = 2^{2^5} + 1$ was the first one Euler showed to be composite:

```mathematica
In[23]:= 2^(2^5) + 1
Out[23]= 4294967297

In[24]:= FactorInteger[2^(2^5) + 1]
Out[24]= {{641, 1}, {6700417, 1}}
```

`FactorInteger` itself is a hybrid driver: small primes by trial division, then Pollard rho, Pollard P–1, Williams P+1, SQUFOF, CFRAC, and finally [GMP-ECM](https://gitlab.inria.fr/zimmermann/ecm) (the elliptic-curve method, which is bundled in `src/external/ecm/`).

`GCD` and `Mod` extend to bignums automatically. There is a pretty identity: $\gcd(2^a - 1, 2^b - 1) = 2^{\gcd(a,b)} - 1$. We can check it numerically:

```mathematica
In[25]:= GCD[2^200 - 1, 2^120 - 1]
Out[25]= 1099511627775

In[26]:= 2^GCD[200, 120] - 1
Out[26]= 1099511627775
```

Modular exponentiation has its own primitive — `PowerMod` — that uses GMP's `mpz_powm` rather than computing a giant power and reducing at the end. A toy RSA round-trip with key $(n=3233, e=17, d=2753)$:

```mathematica
In[27]:= PowerMod[65, 17, 3233]
Out[27]= 2790

In[28]:= PowerMod[2790, 2753, 3233]
Out[28]= 65

In[29]:= PowerMod[17, -1, 3120]
Out[29]= 2753
```

The negative second argument is a modular inverse — `PowerMod[a, -1, m]` solves $a\,d \equiv 1 \pmod{m}$ via the extended Euclidean algorithm. The same primitive recovers the private key 2753 directly from the public exponent and the totient.

A more subtle test: the Carmichael number 561 is composite but passes Fermat's little theorem for every base coprime to it.

```mathematica
In[30]:= PrimeQ[561]
Out[30]= False

In[31]:= PowerMod[2, 560, 561]
Out[31]= 1
```

Fermat alone says `True`; a real primality test (Miller–Rabin) says `False`.

## `N[expr]` and `N[expr, n]`

`N` is the bridge between exact and inexact. With one argument it produces a machine-precision double. With two arguments — `N[expr, n]` — it evaluates to `n` decimal digits of working precision using MPFR. The named constants `Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, and `Degree` all have correctly-rounded MPFR fillers behind them:

```mathematica
In[32]:= N[Pi, 100]
Out[32]= 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170681

In[33]:= N[E, 100]
Out[33]= 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274

In[34]:= N[GoldenRatio, 60]
Out[34]= 1.618033988749894848204586834365638117720309179805762862135448

In[35]:= N[EulerGamma, 60]
Out[35]= 0.5772156649015328606065120900824024310421593359399235988057672

In[36]:= N[Catalan, 60]
Out[36]= 0.9159655941772190150546035149323841107741493742816721342664983
```

`N` plumbs through to elementary functions (`Sqrt`, `Sin`, `Cos`, `Log`, `Exp`, ...) so any closed-form expression composed of them can be evaluated to arbitrary precision. The exact-trig table for rational multiples of $\pi$ continues to apply before MPFR takes over:

```mathematica
In[37]:= Tan[Pi/3]
Out[37]= Sqrt[3]

In[38]:= N[Tan[Pi/3], 50]
Out[38]= 1.73205080756887729352744634150587236694280525381038

In[39]:= N[Sin[1], 80]
Out[39]= 0.841470984807896506652502321630298999622563060798371065672751709991910404391239667
```

## Machin's formula

A nice cross-check: Machin's 1706 formula

$$\frac{\pi}{4} = 4\arctan\!\left(\frac{1}{5}\right) - \arctan\!\left(\frac{1}{239}\right)$$

agrees with the MPFR constant Pi to all 200 digits requested:

```mathematica
In[40]:= N[16 ArcTan[1/5] - 4 ArcTan[1/239], 200]
Out[40]= 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038198

In[41]:= N[Pi, 200]
Out[41]= 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038195
```

The two strings agree to 199 places — the last digit differs because both are rounded independently from a longer internal computation.

## Mixing exact and inexact

The same expression evaluated three different ways gives three different answers. Stirling's approximation $n! \approx \sqrt{2\pi n}\,(n/e)^n$ is a good comparator:

```mathematica
In[42]:= N[100!, 30]
Out[42]= 9.332621544394415268169923885628e+157

In[43]:= N[Sqrt[2 Pi 100] (100/E)^100, 30]
Out[43]= 9.324847625269343247764756127666e+157

In[44]:= N[100!, 30] / N[Sqrt[2 Pi 100] (100/E)^100, 30]
Out[44]= 1.000833677872012141849827856783
```

The ratio is $1 + 1/(12\,n) + O(1/n^2) \approx 1.0008333\ldots$ — exactly Stirling's first correction term.

## A look at chaos

The logistic map $x_{n+1} = r\,x_n(1-x_n)$ is the standard demonstration of sensitivity to initial conditions. With $r = 399/100$ and $x_0 = 2/5$ the orbit, kept exact, blows up in size:

```mathematica
In[45]:= NestList[399/100 # (1 - #) &, 4/10, 3]
Out[45]= {2/5, 1197/1250, 25312959/156250000, 1322447176215313281/2441406250000000000}
```

The numerator's bit length roughly doubles per step. After a few dozen iterations the rationals become unwieldy and you have to make a choice — drop to machine precision and lose digits to chaos, or keep them exact and watch GMP grind through ever-larger integers. Mathilda happily supports either choice.

## What Mathilda does *not* do — significance arithmetic

This is the most important caveat. Mathematica's `N[expr, n]` does not just produce a number with `n` decimal digits — it produces a number that *carries* its precision through subsequent operations and *loses* precision when those operations are ill-conditioned. Subtracting two nearly-equal high-precision numbers in Mathematica gives back a low-precision number, automatically.

Mathilda does not do this. Once an MPFR number has been created at `n` digits, every arithmetic operation runs at `n` digits. The reported `Precision` of the result will read the same `n` whether the operation was harmless or catastrophic:

```mathematica
In[46]:= x = N[E, 50]
Out[46]= 2.71828182845904523536028747135266249775724709369996

In[47]:= Precision[x]
Out[47]= 50.272

In[48]:= Precision[x - 27/10]
Out[48]= 50.272

In[49]:= N[x - 27/10, 50]
Out[49]= 0.0182818284590452353602874713526624977572470936999662
```

In Mathematica, `Precision[x - 27/10]` would drop noticeably because the subtraction destroys leading digits. In Mathilda the precision is the precision *of the representation*, not the precision *of the value*. The user is responsible for tracking error.

The deliberate choice here is simplicity. Significance arithmetic is delicate to implement correctly (Mathematica gets it subtly wrong in places too), it interacts awkwardly with symbolic evaluation, and its absence is rarely a problem if you write your computations in stable form to begin with. Where you need rigorous error bounds, computing with a generous extra precision and checking the answer at half is usually enough.

## Underneath

The numeric stack rests on three vendored or linked libraries:

- **The C standard `<math.h>`** — `sqrt`, `sin`, `log`, etc. for `EXPR_REAL` (machine doubles).
- **[GMP](https://gmplib.org/)** — `mpz_t` for arbitrary-precision integers (`EXPR_BIGINT`), and `mpz_powm`, `mpz_gcd`, `mpz_probab_prime_p` for `PowerMod`, `GCD`, `PrimeQ`.
- **[MPFR](https://www.mpfr.org/)** — correctly-rounded multi-precision floating point, gated behind `USE_MPFR=1` in the makefile, used by `N[expr, n]` and the constants registry in `src/numeric.c`.
- **[GMP-ECM](https://gitlab.inria.fr/zimmermann/ecm)** — bundled in `src/external/ecm/`, providing the elliptic-curve method as the deepest tier of `FactorInteger`.

The interesting glue is in `src/arithmetic.c` (overflow detection and bigint promotion), `src/numeric.c` (MPFR fillers, the constants table, the digit-to-bit conversion factor `log2(10) ≈ 3.32193`), and the MPFR-aware paths in `plus.c`, `times.c`, `power.c`, the trig and hyperbolic modules, and `logexp.c`. A flag (`USE_MPFR`) at compile time toggles MPFR support entirely, so Mathilda still builds and runs without it — `N[expr, n]` then prints a warning and falls back to machine precision.

## In summary

Mathilda gives you three numeric regimes — fast machine doubles, exact integers and rationals via GMP, and arbitrary-precision reals via MPFR — joined by a single bridge, `N`. The boundaries are explicit: you opt in to floating-point by writing a literal with a decimal point or calling `N`, and you stay in the exact world as long as you don't. There is no significance arithmetic — what you ask for is what you get — and the underlying libraries (GMP, MPFR, GMP-ECM, the C math library) do the heavy lifting. For a system whose binary weighs in at a few hundred kilobytes, you can take it surprisingly far.


