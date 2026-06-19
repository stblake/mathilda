# Number Theory

Integer and number-theoretic functions: greatest common divisor and least common multiple (`GCD`, `LCM`, `ExtendedGCD`), modular arithmetic (`PowerMod`, `MultiplicativeOrder`, `PrimitiveRoot`), primality and factorization (`PrimeQ`, `PrimePi`, `NextPrime`, `FactorInteger`, `SquareFreeQ`), arithmetic functions (`EulerPhi`), and continued fractions (`ContinuedFraction`, `FromContinuedFraction`).

## GCD, LCM

- `GCD[n1, n2, ...]`: Greatest common divisor.
- `LCM[n1, n2, ...]`: Least common multiple.

Both fold through GMP when any argument is a bigint, so results that exceed `int64` (e.g. `LCM[20!, 10^100 + 3]`) are computed exactly rather than left symbolic.

## Divisible

- `Divisible[n, m]`: yields `True` if `n` is divisible by `m`, and `False` otherwise. Attributes: `Listable`, `Protected`.

`n` is divisible by `m` when `n` is an integer multiple of `m`; this is effectively `Mod[n, m] == 0`. `Divisible[n, m]` returns `False` unless `n` and `m` are manifestly divisible.

**Features**:
- Machine integers and GMP bigints: tested directly with `mpz_divisible_p`, so large cases such as `Divisible[10^3000 + 1, 16001]` → `True` are exact. By the GMP convention, divisibility by `0` holds iff `n == 0` (`Divisible[0, 0]` → `True`, `Divisible[6, 0]` → `False`); sign is ignored (`Divisible[10, -2]` → `True`).
- Gaussian integers, rationals, and exact numeric quantities: the quotient `n/m` is formed and evaluated; the result is `True` iff it reduces to an integer or a Gaussian integer. So `Divisible[3 + I, 1 - I]` → `True`, `Divisible[3/2, 1/2]` → `True`, `Divisible[2 Pi, Pi/2]` → `True`, while `Divisible[Sqrt[6], Sqrt[2]]` → `False`.
- `Listable`: threads element-wise over lists, e.g. `Divisible[{1, 2, 3, 4, 5, 6}, 2]` → `{False, True, False, True, False, True}`.
- Symbolic, non-numeric arguments leave the call unevaluated (e.g. `Divisible[x, 2]`).
- Diagnostics: too few arguments emit `Divisible::argm`, too many emit `Divisible::argt`; both leave the call unevaluated.

```
In[1]:= Divisible[10, 2]
Out[1]= True

In[2]:= Divisible[5, 2]
Out[2]= False

In[3]:= Divisible[3 + I, 1 - I]
Out[3]= True

In[4]:= Divisible[2 Pi, Pi/2]
Out[4]= True

In[5]:= Divisible[Sqrt[6], Sqrt[2]]
Out[5]= False

In[6]:= Divisible[{1, 2, 3, 4, 5, 6}, 2]
Out[6]= {False, True, False, True, False, True}
```

## CoprimeQ

- `CoprimeQ[n1, n2, ...]`: yields `True` if the arguments are pairwise relatively prime, and `False` otherwise. Attributes: `Listable`, `Orderless`, `Protected`.

Integers are relatively prime when their greatest common divisor is `1`. `CoprimeQ` returns `False` unless the arguments are manifestly coprime integers or Gaussian integers.

**Features**:
- Machine integers and GMP bigints, handled uniformly through `mpz_gcd`, so large cases are exact: `CoprimeQ[2^100 - 1, 3^100 - 1]` → `False` (both even), `CoprimeQ[2^127 - 1, 2^61 - 1]` → `True`. Sign is ignored; `GCD(0, n) = |n|`, so `CoprimeQ[0, 1]` → `True` but `CoprimeQ[0, 5]` → `False`.
- More than two arguments are tested *pairwise*: `CoprimeQ[6, 35, 143]` → `True`, while `CoprimeQ[2, 3, 4]` → `False` (2 and 4 share a factor).
- Gaussian integers: with `GaussianIntegers -> True`, or when any argument is an exact Gaussian integer, coprimality is tested over `Z[i]` via the Gaussian Euclidean algorithm (round-to-nearest division). `CoprimeQ[5 + I, 1 - I]` → `False` (both divisible by `1 + I`); `CoprimeQ[2, 5, GaussianIntegers -> True]` → `True`, while `CoprimeQ[2, 10, GaussianIntegers -> True]` → `False`.
- `Orderless`: argument order is irrelevant, and the `GaussianIntegers` option may appear at any position.
- `Listable`: threads element-wise over lists, e.g. `CoprimeQ[{1, 2, 3, 4, 5}, 6]` → `{True, False, False, False, True}`.
- As a `*Q` predicate it always returns a Boolean: `CoprimeQ[]` → `False`, `CoprimeQ[n]` → `True` (no pairs), and anything not a manifestly coprime integer or Gaussian integer — rationals, reals, symbols, malformed options — yields `False` (e.g. `CoprimeQ[a, b]` → `False`).

```
In[1]:= CoprimeQ[8, 11]
Out[1]= True

In[2]:= CoprimeQ[2, 4]
Out[2]= False

In[3]:= CoprimeQ[2, 3, -5, 7]
Out[3]= True

In[4]:= CoprimeQ[5 + I, 1 - I]
Out[4]= False

In[5]:= CoprimeQ[{1, 2, 3, 4, 5}, 6]
Out[5]= {True, False, False, False, True}
```

## ExtendedGCD

- `ExtendedGCD[n1, n2, ...]`: returns `{g, {r1, r2, ...}}` where `g == GCD[n1, ...]` and `g == r1 n1 + r2 n2 + ...` (a multi-argument Bézout / extended GCD). Attributes: `Listable`, `Protected` (deliberately *not* `Flat`/`Orderless`/`OneIdentity`, since the cofactor list is positional).

**Features**:
- Integer-only; both machine integers and GMP bigints are accepted, and cofactors demote back to machine integers when they fit.
- Computed by folding GMP's `mpz_gcdext` pairwise: `gcd(running_g, n_i) = s·running_g + t·n_i`, scaling the accumulated cofactors by `s` and appending `t` each step. The running gcd stays non-negative, so `g` is exactly `GCD` and the sign convention matches Mathematica.
- Threads element-wise over lists, e.g. `ExtendedGCD[3, {5, 15}]` → `{{1, {2, -1}}, {3, {1, 0}}}`.
- Edge cases: `ExtendedGCD[]` → `{0, {}}`; `ExtendedGCD[n]` → `{|n|, {±1}}`; zeros and negatives handled (`g` non-negative).
- Diagnostics: an inexact (Real) argument emits `ExtendedGCD::exact`; an exact non-integer (Rational) argument emits `ExtendedGCD::egcd`; symbolic arguments leave the call unevaluated silently.

## PowerMod

Gives modular exponentiations, inverses, and roots.
- `PowerMod[a, b, m]`: Gives $a^b \pmod m$.
- `PowerMod[a, -1, m]`: Finds the modular inverse of `a` modulo `m`.
- `PowerMod[a, 1/r, m]`: Finds a modular `r`-th root of `a`.

**Features**:
- `Protected`, `Listable`.
- Evaluates much more efficiently than `Mod[a^b, m]`.
- Integer-exponent path uses GMP `mpz_powm` / `mpz_invert`; `a`, `b`, `m`
  may all be arbitrary-precision bignums.
- Rational-exponent path `PowerMod[a, p/q, m]` solves the modular `q`-th
  root of `a^p` via Tonelli–Shanks, the coprime closed form, Hensel
  lifting, and CRT over the factorisation of `m`. The
  Adleman–Manders–Miller branch (`gcd(q, p-1) > 1` for `q > 2` over a
  large prime) is not yet supported.
- Returns unevaluated if the corresponding inverse or root does not exist
  (or cannot be computed by the implemented algorithms).
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

In[5]:= PowerMod[100, 1/2, 17 * 19 * 23]
Out[5]= 10

In[6]:= PowerMod[2, 1/2, 10^18 + 9]
Out[6]= 742174169206529574
```

## MultiplicativeOrder

Gives the multiplicative order of `k` modulo `n` — the smallest positive
integer `m` such that $k^m \equiv 1 \pmod{n}$.
- `MultiplicativeOrder[k, n]`: order of `k` modulo `n`.
- `MultiplicativeOrder[k, n, {r1, r2, ...}]`: smallest positive `m` with
  $k^m \equiv r_i \pmod{n}$ for some `i` (a multi-target discrete log).

**Features**:
- `Protected`.
- All arithmetic uses GMP `mpz_t`, so `k`, `n`, and any `r_i` may be
  arbitrary-precision bignums.
- Negative `n` is treated as `|n|`; `k` (and each `r_i`) is reduced modulo
  `n` before the search, so negative or out-of-range inputs work
  transparently.
- Returns unevaluated when `gcd(k, n) != 1`, when `n` is zero, or when
  no power of `k` lands in the residue set (3-arg form).
- The 2-arg form factors `phi(n)` and successively strips prime factors
  from `phi(n)` whose corresponding exponent still maps `k` to 1 — so
  the work is dominated by factoring `phi(n)`, not by walking the orbit.
- The 3-arg form walks `k^m mod n` for `m = 1, ..., order(k, n)`. To
  guard against pathological group sizes, the call returns unevaluated
  if the order exceeds `10^8` or does not fit in an unsigned long.
- Non-integer numeric inputs (`Real`, `Complex`, `Rational`) and symbolic
  arguments flow through unevaluated with no diagnostic.
- Diagnostic: `MultiplicativeOrder::argt` when called with anything other
  than 2 or 3 arguments.

```mathematica
In[1]:= MultiplicativeOrder[5, 8]
Out[1]= 2

In[2]:= MultiplicativeOrder[5, 7]
Out[2]= 6

In[3]:= MultiplicativeOrder[-5, 7]
Out[3]= 3

In[4]:= MultiplicativeOrder[5, 7, {3, 11}]
Out[4]= 2

In[5]:= MultiplicativeOrder[10^10000, 7919]
Out[5]= 3959

In[6]:= Select[Range[43], MultiplicativeOrder[#, 43] == EulerPhi[43] &]
Out[6]= {3, 5, 12, 18, 19, 20, 26, 28, 29, 30, 33, 34}

In[7]:= MultiplicativeOrder[10, 22]
Out[7]= MultiplicativeOrder[10, 22]
```

## PrimitiveRoot

Gives a primitive root of `n`, i.e. a generator of the multiplicative group
of integers modulo `n` that are coprime to `n`.
- `PrimitiveRoot[n]`: a primitive root of `n`.
- `PrimitiveRoot[n, k]`: the smallest primitive root of `n` greater than or
  equal to `k`.

**Features**:
- `Protected`, `Listable`.
- Returns unevaluated unless `n` is 2, 4, an odd prime power $p^k$, or
  twice an odd prime power $2 p^k$ (the moduli for which $(\mathbb{Z}/n\mathbb{Z})^*$
  is cyclic). For all other `n`, the call is left unevaluated.
- The 1-argument form returns a canonical primitive root: smallest for
  $n \in \{2, 4\}$ and odd prime powers; for $n = 2 p^k$ the formula
  $g$ if $g$ is odd else $g + p^k$ is applied, where $g$ is the smallest
  primitive root of $p^k$. This matches Mathematica's convention so that,
  e.g. `PrimitiveRoot[10] == 7` while `PrimitiveRoot[10, 1] == 3`.
- The 2-argument form walks forward from `k`; if `k > n - 1` the call is
  left unevaluated.
- All arithmetic uses GMP `mpz_t`, so machine integers, bignums, and
  symbolic bignum products like `Prime[1000000]^1000000` are handled
  uniformly. The prime-power detection iteratively strips prime exponents
  via `mpz_root`, which runs in $O(\omega(k))$ root extractions.
- Diagnostics:
  - `PrimitiveRoot::argt` if not called with 1 or 2 arguments.
  - `PrimitiveRoot::intg` if `n` (or the 2nd-arg `k` when numeric) is not
    an integer greater than 1.

```mathematica
In[1]:= PrimitiveRoot[9]
Out[1]= 2

In[2]:= PrimitiveRoot[10]
Out[2]= 7

In[3]:= PrimitiveRoot[10, 1]
Out[3]= 3

In[4]:= PrimitiveRoot[10, 4]
Out[4]= 7

In[5]:= PrimitiveRoot[{9, 7, 19}]
Out[5]= {2, 3, 2}

In[6]:= PrimitiveRoot[12]
Out[6]= PrimitiveRoot[12]
```

## PrimitiveRootList

Gives the sorted list of all primitive roots of `n` in the canonical
residues $\{1, \ldots, n-1\}$.
- `PrimitiveRootList[n]`

**Features**:
- `Protected`, `Listable`.
- Returns `{}` if `n` is not 2, 4, an odd prime power, or twice an odd
  prime power.
- Enumerates the $\varphi(\varphi(n))$ primitive roots as $g^i \bmod n$
  for $i \in \{1, \ldots, \varphi(n)\}$ with $\gcd(i, \varphi(n)) = 1$,
  where $g$ is the smallest primitive root of `n`. The list is sorted
  ascending.
- Falls back to unevaluated when $\varphi(n)$ exceeds `unsigned long`,
  since the enumeration cannot be represented.
- Non-integer numeric inputs (e.g. `11.0`, `11 + I`) flow through
  unevaluated with no diagnostic, matching Mathematica.
- Diagnostic: `PrimitiveRootList::argx` if not called with exactly 1
  argument.

```mathematica
In[1]:= PrimitiveRootList[9]
Out[1]= {2, 5}

In[2]:= PrimitiveRootList[19]
Out[2]= {2, 3, 10, 13, 14, 15}

In[3]:= PrimitiveRootList[12]
Out[3]= {}

In[4]:= Union[Table[PowerMod[2, i, 9], {i, 6}]]
Out[4]= {1, 2, 4, 5, 7, 8}
```

## PrimeQ

- `PrimeQ[n]`: Returns `True` if `n` is a prime number, `False` otherwise.
- `PrimeQ[z]`: For a Gaussian integer `z = a + b I`, returns `True` if `z` is a Gaussian prime.
- `PrimeQ[n, GaussianIntegers -> True]`: Tests primality of `n` in the
  Gaussian integers `Z[i]`. A rational integer `n` is a Gaussian prime
  iff `|n|` is prime in `Z` AND `n ≡ 3 (mod 4)` (primes `≡ 1 (mod 4)`
  split as `(a+bi)(a-bi)`, and `2 = -i(1+i)^2` is reducible).

**Features**:
- `Listable`, `Protected`.
- Always returns `True` or `False`. For non-integer / non-Gaussian
  inputs (symbols, `Sqrt[2]`, `Exp[2 Pi I/3]`, strings, etc.) returns
  `False` — `*Q` predicates never remain symbolic.
- A Gaussian integer `a + b I` is a Gaussian prime if:
  - Both `a` and `b` are nonzero and `a^2 + b^2` is an ordinary prime, or
  - One of `a`, `b` is zero and the absolute value of the other is a prime congruent to 3 mod 4.

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

In[7]:= PrimeQ[5, GaussianIntegers -> True]
Out[7]= False

In[8]:= PrimeQ[3, GaussianIntegers -> True]
Out[8]= True

In[9]:= PrimeQ[Exp[2 Pi I/3]]
Out[9]= False
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
- Accepts arbitrary-precision integers (`BigInt`). Factorization runs in GMP
  through the same trial-division / Pollard-rho / ECM cascade used by
  `FactorInteger`, so inputs of cryptographic size are tractable.
- For a prime decomposition $n = \prod p_i^{k_i}$, computes
  $\phi(n) = n \prod (1 - 1/p_i)$ as $(n / \prod p_i) \prod (p_i - 1)$
  with exact integer arithmetic.

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4

In[2]:= EulerPhi[2^89 - 1]
Out[2]= 618970019642690137449562110
```

## SquareFreeQ

- `SquareFreeQ[expr]`: Returns `True` if `expr` is a square-free polynomial or
  number, and `False` otherwise.
- `SquareFreeQ[expr, vars]`: Returns `True` if `expr` is square-free with
  respect to the variables `vars` (a `Symbol` or a `List` of `Symbol`s).
- `SquareFreeQ[..., GaussianIntegers -> True | False | Automatic]`: Tests
  square-freeness in `Z[i]` when `True`, in `Z` when `False`. `Automatic`
  (the default) enables Gaussian mode iff the input is a `Complex` literal
  with integer real and imaginary parts.

**Features**:
- `Protected`. Not `Listable` -- passing a list of inputs treats the list as
  the expression (`SquareFreeQ[{1, 2, 3}]` returns `False`).
- Always returns `True` or `False` on a structurally valid call. For inputs
  that are neither a recognised number nor a manifest polynomial -- reals,
  `Sqrt[2]`, `Sin[x]`, `Pi`, strings -- the result is `False`, never
  symbolic.
- An integer `n` is square-free iff `|n|` has no rational prime factor of
  multiplicity `>= 2`. `0` is not square-free; `+/-1` and `+/-p` (for any
  prime `p`) are.
- A rational `p/q` is square-free iff both `p` and `q` are square-free
  integers.
- A polynomial in `vars` is square-free iff for every variable `x_i` in
  `vars` that the polynomial actually depends on, `PolynomialGCD(p, dp/dx_i)`
  has degree `0` in `x_i`. Implementation routes the derivative through
  the `D` builtin and the gcd through `PolynomialGCD`.
- For `GaussianIntegers -> True` (or `Automatic` on a `Complex[Integer,
  Integer]` input), the test factors `N(z) = a^2 + b^2` over `Z` and
  classifies each rational prime by residue `mod 4`:
  - `p == 2`: the Gaussian prime above `2` is `1 + I`; its multiplicity in
    `z` equals the multiplicity of `2` in `N(z)`.
  - `p ≡ 3 (mod 4)`: `p` itself is a Gaussian prime with `N(p) = p^2`;
    its multiplicity in `z` is half the multiplicity of `p` in `N(z)`.
  - `p ≡ 1 (mod 4)`: `p` splits as `pi * conj(pi)`; ambiguous cases (the
    `e_in_norm == 2` slice) are resolved by computing `pi` via a
    Cornacchia search and stripping the multiplicity directly.
- The `Modulus -> p` option is parsed but only `Modulus -> 0` (the default
  no-modulus path) is wired in; any other value -- non-zero integer,
  non-integer, or symbolic -- emits `SquareFreeQ::modnotimpl` and leaves
  the call unevaluated until a polynomial sqfree-mod-`p` test is added.

Diagnostics:
- `SquareFreeQ[]` emits `SquareFreeQ::argb` and stays unevaluated.
- A non-`Rule` past the optional vars slot (e.g. `SquareFreeQ[1, 2, 3]`)
  emits `SquareFreeQ::nonopt` and stays unevaluated.
- `Modulus -> n` with `n != 0` emits `SquareFreeQ::modnotimpl` and stays
  unevaluated (e.g. `SquareFreeQ[x^2 + 1, Modulus -> 2]`).

```mathematica
In[1]:= SquareFreeQ[10]
Out[1]= True

In[2]:= SquareFreeQ[4]
Out[2]= False

In[3]:= SquareFreeQ[20]
Out[3]= False

In[4]:= SquareFreeQ[3 + 2 I]
Out[4]= True

In[5]:= SquareFreeQ[2, GaussianIntegers -> True]
Out[5]= False

In[6]:= SquareFreeQ[2/3]
Out[6]= True

In[7]:= SquareFreeQ[6 + 6 x + x^2]
Out[7]= True

In[8]:= SquareFreeQ[x^3 - x^2 y]
Out[8]= False

In[9]:= SquareFreeQ[x y^2, x]
Out[9]= True

In[10]:= SquareFreeQ[x y^2, y]
Out[10]= False

In[11]:= SquareFreeQ[10^70 + 3]
Out[11]= True
```

## ContinuedFraction

Gives the simple continued-fraction expansion of a number.
- `ContinuedFraction[x, n]`: a list of the first `n` terms.
- `ContinuedFraction[x]`: all terms determinable from `x`.

The list `{a1, a2, a3, ...}` corresponds to `a1 + 1/(a2 + 1/(a3 + ...))`.

**Features**:
- `Protected`, `Listable`.
- **Exact rationals** (Integer / BigInt / Rational) use the Euclidean
  algorithm and return the canonical terminating form (last term `>= 2`,
  never `{..., k-1, 1}`). A finite rational may yield fewer than `n` terms.
- **Quadratic surds** `Sqrt[d]` with `d` a non-square integer use the
  periodic surd recurrence. Without a count the result is
  `{a1, ..., {b1, ...}}`, where the bracketed block repeats cyclically; with
  a count the periodic sequence is unrolled to exactly `n` terms. (General
  quadratic irrationals are not recognised symbolically — pass an explicit
  `n` to use the numeric path.) The no-count form is declined for a `d` whose
  period would be astronomically long.
- **Inexact reals** (machine `Real` or arbitrary-precision MPFR) yield terms
  only as far as the input precision determines them, tracking the value's
  uncertainty and stopping when the next term is no longer pinned down.
- **Exact symbolic reals with an explicit `n`** (e.g. `Pi`, `Sqrt[E]`,
  `Exp[Pi Sqrt[163]]`) are numericised at adaptively increasing precision
  until `n` terms are confirmed by two consecutive evaluations.
- Left unevaluated for an exact non-rational, non-quadratic value with no
  count, or a non-real numeric value.

```mathematica
In[1]:= ContinuedFraction[47/17]
Out[1]= {2, 1, 3, 4}

In[2]:= ContinuedFraction[Sqrt[13]]
Out[2]= {3, {1, 1, 1, 1, 6}}

In[3]:= ContinuedFraction[Pi, 20]
Out[3]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2}

In[4]:= ContinuedFraction[N[Pi]]
Out[4]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}

In[5]:= ContinuedFraction[Exp[Pi Sqrt[163]], 10]
Out[5]= {262537412640768743, 1, 1333462407511, 1, 8, 1, 1, 5, 1, 4}
```

## FromContinuedFraction

The inverse of `ContinuedFraction`: reconstructs a value from its continued-fraction terms.
- `FromContinuedFraction[{a1, a2, ..., an}]`: gives `a1 + 1/(a2 + 1/(a3 + ... + 1/an))`.
- `FromContinuedFraction[{a1, ..., am, {b1, ..., bk}}]`: gives the exact
  quadratic irrational whose terms begin with the `ai` then cycle through the
  `bi` forever.

**Features**:
- `Protected` (not `Listable` — the argument is the whole term list).
- The `ai` of the finite form may be **symbolic**; the result is the convergent
  `h_n / k_n` built from the fundamental recurrence
  `h_i = a_i h_{i-1} + h_{i-2}`, `k_i = a_i k_{i-1} + k_{i-2}`, kept in nested
  (un-expanded) form — `Together` collapses it to a flat rational.
- The **periodic** form requires all `ai` and `bi` to be integers. The purely
  periodic tail solves the quadratic
  `k_{k-1} x^2 + (k_{k-2} - h_{k-1}) x - h_{k-2} = 0` (h, k the period's
  convergents); its positive root is then pushed through the leading terms by a
  Möbius transform and rationalised to a single `(P + Q Sqrt[R]) / S` in lowest
  terms.
- `FromContinuedFraction[{}]` is `0`; `FromContinuedFraction[{x}]` is `x`.
- Left unevaluated for a non-list argument, a sub-list anywhere but last, an
  empty period block, or non-integer terms in a periodic form.

```mathematica
In[1]:= FromContinuedFraction[{2, 1, 3, 4}]
Out[1]= 47/17

In[2]:= FromContinuedFraction[{a, b, c, d}]
Out[2]= (1 + a b + (a + (1 + a b) c) d)/(b + (1 + b c) d)

In[3]:= FromContinuedFraction[{8, {2, 2, 1, 7, 1, 2, 2, 16}}]
Out[3]= Sqrt[71]

In[4]:= FromContinuedFraction[{{1, 2, 3, 4}}]
Out[4]= 1/15 (9 + 2 Sqrt[39])

In[5]:= FromContinuedFraction[ContinuedFraction[Pi, 3]]
Out[5]= 333/106
```

