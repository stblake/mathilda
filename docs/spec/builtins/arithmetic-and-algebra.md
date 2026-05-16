# Arithmetic and Algebra

## Plus (+)
Symbolic sum.
- `a + b + ...`

**Features**:
- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and collects like terms (e.g., `x + 2x` becomes `3x`).
- Returns `0` if no arguments are provided.
- Returns `Overflow[]` if integer addition overflows or if any argument is `Overflow[]`.

```mathematica
In[1]:= 1 + 2 + x + 2*x
Out[1]= 3 + 3*x
```

## Total
Gives the total of elements in a list.
- `Total[list]`
- `Total[list, n]`
- `Total[list, {n}]`
- `Total[list, {n1, n2}]`

**Features**:
- `Protected`.
- `Total[list]` is equivalent to `Apply[Plus, list]`.
- `Total[list, n]` totals all elements down to level `n`.
- `Total[list, {n}]` totals elements at level `n` only.
- Supports negative levels to count from the bottom (`-1` is the last dimension).
- Handles ragged arrays correctly by summing from the inside out when multiple levels are specified.
- `Total[list, Infinity]` totals all atoms in the expression.

```mathematica
In[1]:= Total[{a, b, c, d}]
Out[1]= a + b + c + d

In[2]:= Total[{{1, 2}, {3, 4}}]
Out[2]= {4, 6}

In[3]:= Total[{{1, 2}, {3, 4}}, 2]
Out[3]= 10

In[4]:= Total[{{1, 2}, {3}}, 2]
Out[4]= 6
```

## Accumulate
Gives the running cumulative totals of the elements in a list.
- `Accumulate[list]`
- `Accumulate[list, Method -> "CompensatedSummation"]`

**Features**:
- `Protected`.
- `Accumulate[list]` has the same length as `list`, and is effectively equivalent to `FoldList[Plus, list]`.
- The head of the input is preserved, so `Accumulate[f[a, b, c]]` returns `f[a, a + b, a + b + c]`.
- Threads naturally over rows via `Listable` `Plus`, so `Accumulate` of a matrix accumulates within columns.
- Works on machine integers, GMP arbitrary-precision integers, machine-precision doubles, and symbolic expressions.
- With `Method -> "CompensatedSummation"`, Kahan compensated summation is used in double precision when every element is a machine number, reducing floating-point round-off error. For symbolic or mixed input the option is silently ignored and the standard symbolic accumulation is returned.

```mathematica
In[1]:= Accumulate[{a, b, c, d}]
Out[1]= {a, a + b, a + b + c, a + b + c + d}

In[2]:= Accumulate[{{a, b}, {c, d}, {e, f}}]
Out[2]= {{a, b}, {a + c, b + d}, {a + c + e, b + d + f}}

In[3]:= Accumulate[f[a, b, c, d]]
Out[3]= f[a, a + b, a + b + c, a + b + c + d]

In[4]:= Accumulate[{1, 2, 3, 4, 5}]
Out[4]= {1, 3, 6, 10, 15}

In[5]:= Accumulate[{1.0, 2.0, 3.0}, Method -> "CompensatedSummation"]
Out[5]= {1., 3., 6.}
```

## Times (*)
Symbolic product.
- `a * b * ...`

**Features**:
- `Flat`, `Orderless`, `Listable`.
- Combines numeric constants and groups identical bases into `Power` expressions.
- Handles `I` as `Complex[0, 1]`.
- Returns `1` if no arguments are provided.
- Returns `Overflow[]` if integer multiplication overflows or if any argument is `Overflow[]`.

```mathematica
In[1]:= 2 * x * 3 * y
Out[1]= 6*x*y

In[2]:= I * I
Out[2]= -1
```

## Power (^)
Exponentiation.
- `base ^ exp`

**Features**:
- `Listable`.
- Simplifies integer powers of integers.
- Returns `Overflow[]` if the result exceeds 64-bit integer limits.
- Reduces radicals (e.g., `8^(1/2)` becomes `2*Sqrt[2]`).
- Supports complex results for negative bases (e.g., `(-1)^(1/2)` becomes `I`).
- Distributes power over product if the exponent is an integer.

```mathematica
In[1]:= Sqrt[45]
Out[1]= 3*Sqrt[5]

In[2]:= (a * b)^2
Out[2]= a^2 * b^2
```

## Sqrt
Square root.
- `Sqrt[z]`: Internally represented as `Power[z, 1/2]`.

## Numerator
Gives the numerator of an expression.
- `Numerator[expr]`

**Features**:
- `Protected`, `Listable`.
- Picks out terms which do not have superficially negative exponents.
- Can be used on rational and complex numbers.

```mathematica
In[1]:= Numerator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-2 + x) (-1 + x)

In[2]:= Numerator[3/7 + I/11]
Out[2]= 33 + 7 I
```

## Denominator
Gives the denominator of an expression.
- `Denominator[expr]`

**Features**:
- `Protected`, `Listable`.
- Picks out terms which have superficially negative exponents.
- Can be used on rational and complex numbers.

```mathematica
In[1]:= Denominator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-3 + x)^2

In[2]:= Denominator[3/7 + I/11]
Out[2]= 77
```

## Cancel
Cancels out common factors in the numerator and denominator of an expression.
- `Cancel[expr]`
- `Cancel[expr, Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Threads over equations, inequalities, logic functions, and sums dynamically.
- Evaluates greatest common divisors via polynomial GCD derivations avoiding extraneous expansions.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `Sqrt[y]`, `y^(1/3)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial cancellation in `g`, then substitutes back.
- The algebraic-generator pass runs `Together` on the substituted form (not just GCD-cancellation), so inputs whose `g`-substituted denominator is a Plus of terms with different `g`-denominators (e.g. `1/(g^2 - 1/g)` from `1/(y^(2/3) - 1/y^(1/3))`) are handled correctly.
- Extracts algebraic-constant atoms (`Sqrt[2]`, `Sqrt[3]`, `CubeRoot[5]`, `2^(2/3)`, ...) that appear in every summand of *both* numerator *and* denominator and divides them out before the polynomial GCD step. Closes a long-standing gap where `PolynomialGCD[Sqrt[2], Sqrt[2] + Sqrt[2] x^4]` returns 1 (the integer-content recursion treats `Sqrt[2]` as having content 1), so cancellations whose only shared factor was an algebraic constant survived as-is. The pass is intentionally narrow — only `Power[integer, rational/non-integer]` factors are eligible — to avoid disturbing the rational-function intermediates that the integration dispatcher pattern-matches against.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) cancels common factors over `Q(alpha)` instead of `Q`. Implementation: lifts numerator and denominator into `Q(alpha)[x]` via the QAUPoly machinery, runs `qaupoly_gcd`, divides both sides by `g`, and re-renders. Works for single-fraction inputs; `Plus` inputs (sums of fractions) currently fall back to the no-extension path because `PolynomialQuotient` does not yet accept `Extension` (Phase 0.5 follow-up).

```mathematica
In[1]:= Cancel[(x^2 - 1) / (x - 1)]
Out[1]= 1 + x

In[2]:= Cancel[(x - y)/(x^2 - y^2) + (x^3 - 27)/(x^2 - 9)]
Out[2]= (9 + 3 x + x^2)/(3 + x) + 1/(x + y)

In[3]:= Cancel[(y - 1)/(Sqrt[y] - 1)]
Out[3]= 1 + Sqrt[y]

In[4]:= Cancel[(y - 1)/(y^(1/3) - 1)]
Out[4]= 1 + y^(1/3) + y^(2/3)

In[5]:= Cancel[1/(y^(2/3) - 1/y^(1/3))]
Out[5]= y^(1/3)/(-1 + y)

In[6]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[6]= Sqrt[2] + x

In[7]:= Cancel[(x^3 - 2)/(x - 2^(1/3)), Extension -> 2^(1/3)]
Out[7]= 2^(2/3) + 2^(1/3) x + x^2

In[8]:= Cancel[Sqrt[2]/(Sqrt[2] + Sqrt[2] x^4)]
Out[8]= 1/(1 + x^4)

In[9]:= Cancel[(Sqrt[2] x + Sqrt[2] y)/(Sqrt[2] + Sqrt[2] x)]
Out[9]= (x + y)/(1 + x)
```

## Together
Puts terms in a sum over a common denominator, and cancels factors in the result.
- `Together[expr]`
- `Together[expr, Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Makes a sum of terms into a single rational function.
- Computes lowest common multiples (LCM) of denominators securely without unconditionally destroying pre-factored bases unnecessarily.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `y^(1/3)`, `y^(2/3)`, `y^(73/24)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial pipeline in `g`, then substitutes back.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) combines into a single fraction with the standard combiner, then runs `Cancel[..., Extension -> alpha]` on the result so algebraic-coefficient cancellations fire. Effective on simple inputs like `1/(x - Sqrt[2]) + 1/(x + Sqrt[2])`, which collapses to `(2 x)/(x^2 - 2)`. Inputs whose summands themselves carry algebraic-coefficient denominators are deferred to Phase 0.5 (which will plumb `Extension` through `PolynomialLCM` / `PolynomialQuotient` / `together_recursive`).

```mathematica
In[1]:= Together[a/b + c/d]
Out[1]= (b c + a d) / (b d)

In[2]:= Together[x^2/(x^2 - 1) + x/(x^2 - 1)]
Out[2]= x / (-1 + x)

In[3]:= Together[1/x + 1/(x + 1) + 1/(x + 2) + 1/(x + 3)]
Out[3]= (2 (3 + 11 x + 9 x^2 + 2 x^3)) / (x (1 + x) (2 + x) (3 + x))

In[4]:= Together[x^2/(x - y) - x y/(x - y)]
Out[4]= x

In[5]:= Together[y^(5/8)*(y^(19/8) - y^(73/24)/(y^(2/3) - 1/y^(1/3)))]
Out[5]= -y^3 / (-1 + y)

In[6]:= Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2]), Extension -> Sqrt[2]]
Out[6]= (2 x)/(-2 + x^2)
```

## Apart
Gives the partial fraction decomposition of a rational expression.
- `Apart[expr]`
- `Apart[expr, var]`
- `Apart[expr, var, Extension -> alpha]`

**Features**:
- `Protected`, `Listable`.
- Writes `expr` as a polynomial in `var` together with a sum of ratios of polynomials with minimal denominators.
- If `var` is not specified, intelligently selects the main polynomial variable natively.
- Implements exact undetermined coefficients algebraically leveraging row-reduced identity expansions over algebraic inputs avoiding recursive fractional losses natively.
- When `Together[expr]` produces a numerator or denominator that is not polynomial in the chosen variable (e.g. fractional-power inputs whose Together'd form is `y^(1/3)/(y - 1)`), the matrix-of-coefficients algorithm cannot apply; Apart returns the `Together` form unchanged rather than synthesising a spurious zero.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) factors the denominator over `Q(alpha)` before partial-fraction decomposition runs, splitting reducible-over-extension factors (e.g. `x^2 - 2` into `(x - Sqrt[2])(x + Sqrt[2])` under `Extension -> Sqrt[2]`) and producing the corresponding linear-factor partial fractions. The pre-`Together` step is also extension-aware so any algebraic-number cancellations in numerator/denominator fire before splitting.

```mathematica
In[1]:= Apart[1/((1+x)(5+x))]
Out[1]= 1/(4 (1 + x)) - 1/(4 (5 + x))

In[2]:= Apart[(x^5-2)/((1+x+x^2)(2+x)(1-x))]
Out[2]= 2 - x + (-1 - x/3)/(1 + x + x^2) + 1/(9 (-1 + x)) - 34/(9 (2 + x))

In[3]:= Apart[(x+y)/((x+1)(y+1)(x-y)), x]
Out[3]= 2 y/((1 + y)^2 (x - y)) - (-1 + y)/((1 + x) (1 + y)^2)

In[4]:= Apart[1/(y^(2/3) - 1/y^(1/3))]
Out[4]= y^(1/3)/(-1 + y)

In[5]:= Apart[1/(x^2 - 2), x, Extension -> Sqrt[2]]
Out[5]= -1/2 1/(Sqrt[2] (Sqrt[2] + x)) + 1/4 Sqrt[2]/(-Sqrt[2] + x)
```

## Mod, Quotient, QuotientRemainder
- `Mod[n, m]`: Remainder of `n/m`.
- `Quotient[n, m]`: Integer part of `n/m`.
- `QuotientRemainder[n, m]`: Returns `{Quotient[n, m], Mod[n, m]}`.

## GCD, LCM
- `GCD[n1, n2, ...]`: Greatest common divisor.
- `LCM[n1, n2, ...]`: Least common multiple.

## PowerMod
Gives modular exponentiations, inverses, and roots.
- `PowerMod[a, b, m]`: Gives $a^b \pmod m$.
- `PowerMod[a, -1, m]`: Finds the modular inverse of `a` modulo `m`.
- `PowerMod[a, 1/r, m]`: Finds a modular `r`-th root of `a`.

**Features**:
- `Protected`, `Listable`.
- Evaluates much more efficiently than `Mod[a^b, m]`.
- Returns unevaluated if the corresponding inverse or root does not exist.
- Allows threading over lists natively.

```mathematica
In[1]:= PowerMod[2, 10, 3]
Out[1]= 1

In[2]:= PowerMod[3, -2, 7]
Out[2]= 4

In[3]:= PowerMod[3, 1/2, 2]
Out[3]= 1

In[4]:= PowerMod[2, {10, 11, 12, 13, 14}, 5]
Out[4]= {4, 3, 1, 2, 4}
```

## Factorial (!)
Gives the factorial of an integer or half-integer.
- `n!` or `Factorial[n]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- Evaluates exactly for positive integers up to $20!$.
- Yields `ComplexInfinity` for negative integers.
- Supports half-integers utilizing factors of $\sqrt{\pi}$ recursively.
- Supports trailing `!` parsed natively as a postfix operator.

```mathematica
In[1]:= 5!
Out[1]= 120

In[2]:= (1/2)!
Out[2]= Sqrt[Pi]/2

In[3]:= Factorial[0]
Out[3]= 1
```

## FactorialPower
Falling factorial. `FactorialPower[n, k]` = $n (n-1) (n-2) \cdots (n - k + 1)$.
- `FactorialPower[n, k]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- For non-negative integer $n$ and non-negative integer $k$: exact GMP product.
- For symbolic $n$ with concrete $k \le 32$: expands to an explicit product
  `Times[n, n - 1, ..., n - k + 1]` so `Expand` and `D` can act on it.
- Equivalent to $n! / (n - k)!$ when both arguments are non-negative integers.
- Used as the symbolic-order derivative of a power: $D[x^n, \{x, k\}] = \mathrm{FactorialPower}[n, k]\, x^{n-k}$.

```mathematica
In[1]:= FactorialPower[5, 3]
Out[1]= 60

In[2]:= FactorialPower[n, 3]
Out[2]= n (-2 + n) (-1 + n)

In[3]:= FactorialPower[n, 0]
Out[3]= 1

In[4]:= D[x^n, {x, k}]
Out[4]= x^(-k + n) FactorialPower[n, k]
```

## Binomial
Gives the binomial coefficient $\binom{n}{m}$.
- `Binomial[n, m]`

**Features**:
- `Protected`, `Listable`, `NumericFunction`.
- Evaluates exactly for integers, half-integers, and dynamically factors symbolic terms correctly (e.g. `Binomial[n, 4]`).
- Reduces numerical boundaries logically utilizing continuous `Gamma` interpolations.

```mathematica
In[1]:= Binomial[10, 3]
Out[1]= 120

In[2]:= Binomial[8.5, -4.2]
Out[2]= 0.0000604992

In[3]:= Binomial[9/2, 7/2]
Out[3]= 9/2

In[4]:= Binomial[n, 4]
Out[4]= 1/24 (-3 + n) (-2 + n) (-1 + n) n

In[5]:= Binomial[0, 1]
Out[5]= 0
```

## PrimeQ
- `PrimeQ[n]`: Returns `True` if `n` is a prime number, `False` otherwise.
- `PrimeQ[z]`: For a Gaussian integer `z = a + b I`, returns `True` if `z` is a Gaussian prime.

**Features**:
- `Listable`, `Protected`.
- A Gaussian integer `a + b I` is a Gaussian prime if:
  - Both `a` and `b` are nonzero and `a^2 + b^2` is an ordinary prime, or
  - One of `a`, `b` is zero and the absolute value of the other is a prime congruent to 3 mod 4.
- Returns `False` for non-integer complex numbers, rationals, reals, and strings.

```mathematica
In[1]:= PrimeQ[7]
Out[1]= True

In[2]:= PrimeQ[1 + I]
Out[2]= True

In[3]:= PrimeQ[1 + 2 I]
Out[3]= True

In[4]:= PrimeQ[3 I]
Out[4]= True

In[5]:= PrimeQ[5 I]
Out[5]= False

In[6]:= PrimeQ[2 + 2 I]
Out[6]= False
```

## PrimePi
- `PrimePi[x]`: Returns the number of primes less than or equal to `x`.

**Features**:
- `Listable`, `Protected`.
- Uses Meissel-Lehmer algorithm for efficient counting.

```mathematica
In[1]:= PrimePi[10]
Out[1]= 4

In[2]:= PrimePi[100]
Out[2]= 25

In[3]:= PrimePi[{10, 100}]
Out[3]= {4, 25}
```

## NextPrime
- `NextPrime[x]`: Smallest prime above `x`.
- `NextPrime[x, k]`: $k$-th next prime above `x` (or previous if $k$ is negative).

**Features**:
- `Protected`, `ReadProtected`.
- Supports negative $k$ for finding previous primes.
- Remains unevaluated if no such prime exists (e.g., `NextPrime[2, -1]`).

## FactorInteger
- `FactorInteger[n]`: Returns a list of prime factors and their exponents.
- `FactorInteger[n, k]`: Partial factorization, at most `k` distinct factors.
- `FactorInteger[n, Automatic]`: Pulls out easy factors using trial division.
- `FactorInteger[n, Method -> method]`: Factors `n` using a specific algorithmic method.

**Options for `Method`**:
- `"Automatic"`: (Default) Attempts factorization by sequentially executing Trial Division, Pollard's Rho, and ECM.
- `"TrialDivision"`: Extracts bounds matching the first 1000 computed primes and halts cleanly. Composite residues are preserved.
- `"PollardRho"`: Executes the Brent cycle-finding variant of Pollard's Rho algorithm targeting GMP bignums.
- `"BlakeRationalBaseDescent"`: Executes the Rational Base Descent algorithm against semiprimes `n = p q` where one factor is approximately `c (a/b)^k` for some unknown `c`, `k` and a small coprime integer pair `a > b`. With no explicit `"Base"`, auto-searches coprime partitions `a + b = j` for `j = 3, 4, ...` indefinitely until a factor is found (interruptible).
- `{"BlakeRationalBaseDescent", "Base" -> a/b}`: Same algorithm but pinned to a specific rational base `a/b`. Skips the auto-search.
- `"PollardP-1"`: Executes the Pollard $P-1$ algorithm, leveraging GMP and ECM natively.
- `"WilliamsP+1"`: Executes the Williams $P+1$ algorithm via the ECM library natively.
- `"Fermat"`: Explores factors symmetrically close to the square root boundary natively on large integers.
- `"CFRAC"`: Implements the continued fraction integer factorization method natively on GMP bignums.
- `"ShanksSquareForms"`: Implements Shanks's Square Forms Factorization (SQUFOF) natively on large integers. Halts explicitly if factors are not discovered within the loop constraints.
- `"ECM"`: Explicitly triggers Elliptic-Curve Method discovery natively.

**Features**:
- `Listable`, `Protected`.
- Supports negative integers (includes `{-1, 1}`).
- Supports rational numbers (denominator factors have negative exponents).

```mathematica
In[1]:= FactorInteger[12]
Out[1]= {{2, 2}, {3, 1}}

In[2]:= FactorInteger[-12]
Out[2]= {{-1, 1}, {2, 2}, {3, 1}}

In[3]:= FactorInteger[3/4]
Out[3]= {{2, -2}, {3, 1}}

In[4]:= FactorInteger[100, 1]
Out[4]= {{2, 2}, 25}
```

## EulerPhi
Gives the Euler totient function $\phi(n)$.
- `EulerPhi[n]`

**Features**:
- `Listable`, `Protected`.
- Counts the number of positive integers less than or equal to $n$ that are relatively prime to $n$.
- Returns 0 for $n = 0$, and handles negative integers via $\phi(-n) = \phi(n)$.

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4
```

## Rational
Represents a rational number.
- `Rational[n, d]`

**Features**:
- Automatically simplifies to lowest terms (e.g. `Rational[15, 5]` evaluates to `3`, `Rational[2, 4]` evaluates to `Rational[1, 2]`).
- Returns `Indeterminate` when `n` and `d` are both `0` (e.g. `Rational[0, 0]`).
- Returns `ComplexInfinity` when `n` is non-zero and `d` is `0` (e.g. `Rational[1, 0]`).

```mathematica
In[1]:= Rational[15, 5]
Out[1]= 3
```

## Re, Im, ReIm, Abs, Conjugate, Arg
Complex number functions.
- `Re[z]`: Real part.
- `Im[z]`: Imaginary part.
- `ReIm[z]`: Returns `{Re[z], Im[z]}`.
- `Abs[z]`: Magnitude.
- `Conjugate[z]`: Complex conjugate.
- `Arg[z]`: Phase angle.

## Solve
Attempts to solve an equation or system of equations for one or more variables.
- `Solve[expr, vars]`: Solve `expr` for `vars` over the complex numbers (default).
- `Solve[expr, vars, dom]`: Solve over the domain `dom`. Supported: `Complexes` (default), `Reals`.

**Features**:
- `HoldAll`, `Protected`.
- The initial implementation handles a single equality `lhs == rhs` that is
  polynomial in one variable *after* clearing denominators.  Multi-equation
  systems and inequalities are reserved for a future cut and currently leave
  `Solve[...]` unevaluated.
- Acts as a router that dispatches to `Solve`SolvePolynomialEquality` (see below).
- Returns the solution set as a `List` of `List` of `Rule` pairs:
  - `{}` -- no solutions.
  - `{{}}` -- tautology (full-dimensional solution set).
  - `{{x -> v1}, {x -> v2}, ...}` -- one inner list per solution.  Multiplicity
    is preserved (repeated roots appear once per unit of multiplicity).
- **Rational-equality canonicalisation**: both sides are run through `Together`
  to combine into single fractions `N1/D1 == N2/D2`, then cross-multiplied to
  `N1*D2 - N2*D1 == 0` and `Collect`-ed in the solving variable before
  dispatch.  This routes equations like `a/x + b == 0` or `1/(x-1) == 2`
  through the polynomial specialist.  Any candidate root that provably zeroes
  one of the cleared denominators is dropped as extraneous (e.g.
  `Solve[x/(x-1) == 2/(x-1), x]` returns `{{x -> 2}}`, not `{{x -> 1}, {x -> 2}}`).
  Symbolic / undetermined denominator values are kept (parametric inputs like
  `Solve[a/x + b == 0, x]` return `{{x -> -a/b}}`).
- Per-degree handling for irreducible factors:
  - Degree 1 / 2: closed-form rules.
  - Quadratic in `Reals`: discriminant-aware (zero, one, or two real roots).
  - Binomial `a*x^n + b == 0`: all n complex roots, or the real
    radical(s) in `Reals` (n odd -> 1 root; n even with `−b/a > 0` -> ±r;
    n even with `−b/a < 0` -> 0).
  - n-quadratic `a*x^(2n) + b*x^n + c == 0`: substitution `u = x^n` followed by
    two binomial sub-solves; 2n radical roots regardless of `Cubics` / `Quartics`.
  - Degree 3: held `Root[Function[t, p[t]], k]` objects unless `Cubics -> True`.
  - Degree 4: held `Root[]` objects unless `Quartics -> True` (Ferrari is
    deferred -- with `Quartics -> True` the four roots are emitted as held
    `Root[]` for now).
  - Degree ≥ 5: held `Root[]` objects per irreducible factor.

**Options**:
- `Cubics -> False`: Emit cubic roots as held `Root[]` objects (default).
  `Cubics -> True` switches to closed-form Cardano radicals.
- `Quartics -> False`: Emit quartic roots as held `Root[]` objects (default).
  Reserved: Ferrari closed form is deferred.
- `GeneratedParameters -> C`: Reserved.  Reserved name for newly introduced
  parameters in a future parametric-solution path.
- `VerifySolutions -> Automatic`: Reserved.

```mathematica
In[1]:= Solve[2 x + 3 == 0, x]
Out[1]= {{x -> -3/2}}

In[2]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[2]= {{x -> 2}, {x -> 3}}

In[3]:= Solve[x^2 + 1 == 0, x]
Out[3]= {{x -> -I}, {x -> I}}

In[4]:= Solve[x^2 + 1 == 0, x, Reals]
Out[4]= {}

In[5]:= Solve[(x-1)^2 == 0, x]
Out[5]= {{x -> 1}, {x -> 1}}

In[6]:= Solve[x^4 - 5 x^2 + 4 == 0, x]
Out[6]= {{x -> 1}, {x -> -1}, {x -> 2}, {x -> -2}}

In[7]:= Solve[x^3 + x + 1 == 0, x]
Out[7]= {{x -> Root[Function[1 + #1 + #1^3], 1]}, ...}

In[8]:= Solve[Sin[x] == 0, x]
Out[8]= Solve[Sin[x] == 0, x]

In[9]:= Solve[a/x + b == 0, x]
Out[9]= {{x -> -a/b}}

In[10]:= Solve[1/(x-1) == 2, x]
Out[10]= {{x -> 3/2}}

In[11]:= Solve[x/(x-1) == 2/(x-1), x]
Out[11]= {{x -> 2}}             (* x = 1 dropped as extraneous *)
```

## Solve`SolvePolynomialEquality
The polynomial-equality specialist invoked by `Solve`.  Reachable directly via
its context-qualified name when the caller has already classified its input.
- `Solve`SolvePolynomialEquality[lhs == rhs, var]`
- `Solve`SolvePolynomialEquality[lhs == rhs, var, dom]`

**Features**:
- `Protected`.
- Same algorithm and output shape as `Solve` for single polynomial equalities
  in one variable.  Does not parse options; the caller supplies them through
  the C-level entry point.

## Cubics
Option for `Solve` that controls whether cubic equations are solved via
explicit radical formulas.
- `Cubics -> False` (default): emit held `Root[]` objects.
- `Cubics -> True`: emit closed-form Cardano radicals.

## Quartics
Option for `Solve` that controls whether quartic equations are solved via
explicit radical formulas.
- `Quartics -> False` (default): emit held `Root[]` objects.
- `Quartics -> True`: reserved -- Ferrari is deferred; held `Root[]` objects
  are still emitted in the current implementation.

