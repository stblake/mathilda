# Number Theory

Integer and number-theoretic functions: greatest common divisor and least common multiple (`GCD`, `LCM`, `ExtendedGCD`), modular arithmetic (`PowerMod`, `MultiplicativeOrder`, `PrimitiveRoot`, `ChineseRemainder`), primality and factorization (`PrimeQ`, `PrimePowerQ`, `Prime`, `PrimePi`, `NextPrime`, `FactorInteger`, `SquareFreeQ`), arithmetic functions (`EulerPhi`), partitions (`IntegerPartitions`, `PartitionsP`, `PartitionsQ`), and continued fractions (`ContinuedFraction`, `FromContinuedFraction`).

## GCD, LCM

- `GCD[n1, n2, ...]`: Greatest common divisor.
- `LCM[n1, n2, ...]`: Least common multiple.

Both fold through GMP when any argument is a bigint, so results that exceed `int64` (e.g. `LCM[20!, 10^100 + 3]`) are computed exactly rather than left symbolic.

- **Rationals**: `gcd(p1/q1, p2/q2) = gcd(p1, p2)/lcm(q1, q2)` and `lcm(p1/q1, p2/q2) = lcm(p1, p2)/gcd(q1, q2)`.
- **Gaussian integers and rationals**: any non-real exact complex argument switches to `Z[i]`. Each argument is written `w/d` (`w` in `Z[i]`, `d` a positive integer). The gcd is `gcd(w_i)/lcm(d_i)` and the lcm is `lcm(w_i)/gcd(d_i)`, normalised to the first-quadrant associate (`Re > 0`, `Im >= 0`). This matches the reference on 400 of 400 random pairs: `GCD[3 + I, 2]` → `1 + I`, `GCD[1/2 + I/3, 2]` → `1/6`, `LCM[3 + I, 2]` → `2 + 4 I`, `GCD[-1 - I]` → `1 + I`.
- **Attributes**: `Flat`, `Listable`, `OneIdentity`, `Orderless`, `Protected`, as in the reference. `NumericFunction` was dropped in v0.471.
- **Messages**: an inexact argument draws `GCD::exact` / `LCM::exact`. `LCM[]` draws `LCM::argm` and stays unevaluated; it used to return `1`. `GCD[]` is `0`, as in the reference.

> **Packed arrays.** Over an `int64` buffer both run on the machine words —
> `GCD[list, n]`, `LCM[list, n]` and the two-array forms — and answer with a
> packed result. A product that would overflow `int64` abandons the buffer so
> the ordinary path answers exactly in GMP; never a wrapped value. See
> [Packed arrays](packed-arrays.md).

## Divisible

- `Divisible[n, m]`: yields `True` if `n` is divisible by `m`, and `False` otherwise. Attributes: `Listable`, `Protected`.

`n` is divisible by `m` when `n` is an integer multiple of `m`; this is effectively `Mod[n, m] == 0`. `Divisible[n, m]` returns `False` unless `n` and `m` are manifestly divisible.

**Features**:
- Machine integers and GMP bigints: tested directly with `mpz_divisible_p`, so large cases such as `Divisible[10^3000 + 1, 16001]` → `True` are exact. Sign is ignored (`Divisible[10, -2]` → `True`).
- A zero divisor draws `Divisible::divz` and leaves the call unevaluated, as in the reference (`Divisible[6, 0]` and `Divisible[0, 0]`). This happens only when every argument is a number, so `Divisible[x, 0]` stays unevaluated silently. An inexact argument draws `Divisible::exact`.
- `Divisible[n, m1, m2, ...]` is a **chain**, as in the reference: `n` divisible by `m1`, `m1` by `m2`, and so on. So `Divisible[12, 2, 3]` → `False` and `Divisible[8, 2, 2]` → `True`.
- Packed `int64` buffers (against an exact Integer, or an equal-shape buffer) answer straight off the buffer with a packed `True`/`False` List (`ndint_binary_predicate`). `Divisible[Range[10^6], 7]` takes 9 ms (the reference 0.28 s). A zero divisor degrades to the element path, so `::divz` still fires.
- Gaussian integers, rationals, and exact numeric quantities: the quotient `n/m` is formed and evaluated; the result is `True` iff it reduces to an integer or a Gaussian integer. So `Divisible[3 + I, 1 - I]` → `True`, `Divisible[3/2, 1/2]` → `True`, `Divisible[2 Pi, Pi/2]` → `True`, while `Divisible[Sqrt[6], Sqrt[2]]` → `False`.
- `Listable`: threads element-wise over lists, e.g. `Divisible[{1, 2, 3, 4, 5, 6}, 2]` → `{False, True, False, True, False, True}`.
- Symbolic, non-numeric arguments leave the call unevaluated (e.g. `Divisible[x, 2]`).
- Diagnostics: one argument emits `Divisible::argmu`, none `Divisible::argm`; both leave the call unevaluated.

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
- Gaussian integers: when any argument is an exact Gaussian integer, coprimality is tested over `Z[i]` (the gcd from the shared `gi_gcd` must be a unit). `CoprimeQ[5 + I, 1 - I]` → `False` (both divisible by `1 + I`).
- **No options** (the reference). A `GaussianIntegers -> …` rule is just a non-integer argument, so `CoprimeQ[2, 5, GaussianIntegers -> True]` → `False`. The option was removed in v0.471.
- `Orderless`: argument order is irrelevant.
- Packed `int64` buffers answer straight off the buffer with a packed `True`/`False` List. `CoprimeQ[Range[10^6], 30]` takes 37 ms (the reference 0.37 s). `CoprimeQ[a, b]` compiles as `GCD == 1`.
- `Listable`: threads element-wise over lists, e.g. `CoprimeQ[{1, 2, 3, 4, 5}, 6]` → `{True, False, False, False, True}`.
- As a `*Q` predicate it always returns a Boolean. `CoprimeQ[]` → `False`. `CoprimeQ[n]` asks whether `n` is a **unit**, as in the reference: `CoprimeQ[1]`, `CoprimeQ[-1]` and `CoprimeQ[I]` → `True`, while `CoprimeQ[3]` and `CoprimeQ[0]` → `False`. It used to answer `True` for every single argument. Anything that is not an integer or Gaussian integer (rationals, reals, symbols, rules) yields `False` (e.g. `CoprimeQ[a, b]` → `False`).

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
- Integers: both machine integers and GMP bigints are accepted, and cofactors demote back to machine integers when they fit. Two machine integers skip GMP entirely, and the int64 path reproduces `mpz_gcdext`'s unique cofactors exactly (checked against the GMP fold on 500 random 10¹⁵-sized pairs).
- **Gaussian integers** (any non-real Gaussian argument): a Euclidean algorithm in `Z[i]` that starts from the operand of larger norm and rounds each quotient component to the nearest integer, ties to even (the reference's `Round`). `g` is first-quadrant normalised. `ExtendedGCD[3 + I, 2]` → `{1 + I, {-I, 2 I}}`.
  - **Deliberate difference:** the cofactors match the reference's on 99% of random pairs (397 of 400 in one sample, 399 of 400 in another). The rest are different, equally valid Bézout pairs.
  - A non-integer Gaussian rational draws `ExtendedGCD::egcd`.
- Computed by folding GMP's `mpz_gcdext` pairwise: `gcd(running_g, n_i) = s·running_g + t·n_i`, scaling the accumulated cofactors by `s` and appending `t` each step. The running gcd stays non-negative, so `g` is exactly `GCD` and the sign convention is the conventional one.
- Threads element-wise over lists, e.g. `ExtendedGCD[3, {5, 15}]` → `{{1, {2, -1}}, {3, {1, 0}}}`.
- Edge cases: `ExtendedGCD[]` → `{0, {}}`; `ExtendedGCD[n]` → `{|n|, {±1}}`; zeros and negatives handled (`g` non-negative).
- Diagnostics: an inexact (Real) argument emits `ExtendedGCD::exact`; an exact non-integer (Rational) argument emits `ExtendedGCD::egcd`; symbolic arguments leave the call unevaluated silently.

## ChineseRemainder

- `ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}]`: the smallest `x ≥ 0` satisfying every integer congruence `x mod mi == ri mod mi`, lying in `0 ≤ x < LCM[m1, m2, ...]`. Attributes: `Protected` (deliberately *not* `Listable` — the residue and modulus lists are atomic arguments).
- `ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}, d]`: the smallest such `x ≥ d`, lying in `d ≤ x < d + LCM[m1, m2, ...]`.

**Features**:
- Integer-only; machine integers and GMP bigints share one path, so a residue-number-system recovery over large coprime moduli returns the full bignum (`ChineseRemainder[data, keys]` with prime keys is the CRT-encryption idiom).
- The moduli need **not** be pairwise coprime: a solution exists iff every pair of congruences agrees modulo the gcd of their moduli. When the system is inconsistent, `ChineseRemainder` is left unevaluated, e.g. `ChineseRemainder[{1, 2}, {6, 10}]` (since `1 ≢ 2 (mod gcd(6, 10) = 2)`).
- Computed by a streaming pairwise CRT fold: an accumulator congruence `x ≡ x_acc (mod m_acc)` is merged with each `(ri, |mi|)` via `mpz_gcdext` (`x_acc += ((ri − x_acc)/g)·u·m_acc`, `m_acc = lcm`), reducing modulo the running lcm each step. The final `x = d + ((x_acc − d) mod L)` gives the smallest `x ≥ d` and handles a negative `d` and empty lists (`L = 1`).
- If all `0 ≤ ri < mi` the result satisfies `x mod mi == ri`; unreduced or negative residues are handled (the result gives `x mod mi == ri mod mi`).
- Edge cases and messages (as in the reference):
  - Two lists of numbers that are not equal-length integers with **positive** moduli draw `ChineseRemainder::pilist`. That covers a zero or negative modulus, a length mismatch, and a rational, real or complex entry. Symbolic entries leave the call unevaluated silently.
  - A non-integer offset `d` draws `ChineseRemainder::int`.
  - Empty lists stay unevaluated. They used to give `0`; the reference itself gives `ChineseRemainder[{}, {}, 5]` an `argrx` message, which Mathilda does not copy.
  - A wrong argument count draws `ChineseRemainder::argr` (one argument) or `::argrx`, with the reference's "2 arguments are expected" wording.

```
In[1]:= ChineseRemainder[{3, 4}, {4, 5}]
Out[1]= 19

In[2]:= ChineseRemainder[{2, 3, 5}, {3, 5, 7}]
Out[2]= 68

In[3]:= ChineseRemainder[{2, 3, 5}, {3, 5, 7}, 100]
Out[3]= 173

In[4]:= ChineseRemainder[{1, 2}, {6, 10}]
Out[4]= ChineseRemainder[{1, 2}, {6, 10}]

In[5]:= Mod[ChineseRemainder[{123, 456, 789}, {1000000007, 1000000009, 1000000021}], {1000000007, 1000000009, 1000000021}]
Out[5]= {123, 456, 789}
```

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
- Messages, as in the reference (each leaves the call unevaluated):
  - `PowerMod::ninv` when the base is not invertible;
  - `PowerMod::root` when no modular root exists;
  - `PowerMod::divz` for a zero modulus;
  - `PowerMod::pmod` for an inexact argument, a non-integer base or modulus, or a non-rational exponent;
  - `PowerMod::argr` / `::argrx` for a wrong argument count.
- **A negative integer modulus** stays unevaluated, as in the reference. It used to be read as `|m|`. The packed path declines it, and `Compile` already bailed.
- **Gaussian integers** (integer exponent): a Gaussian base and/or modulus is reduced with `Mod`'s nearest quotient after every product (`gi_powmod`): `PowerMod[1 + I, 3, 5]` → `-2 + 2 I`, `PowerMod[2, 10, 2 + I]` → `-1`.
  - This agrees with the reference on 399 of 400 random cases. The remaining one is a rounding tie in the final reduction, with an equally valid representative.
  - A machine route covers modulus components below 2¹⁶.
  - **Deliberate difference:** the reference leaves a few invertible Gaussian bases unevaluated under a negative exponent; Mathilda answers them.
- Machine fast path: a machine base, a non-negative exponent and a modulus below 2³² square-and-multiply in `uint64`. The `Do`-loop of 10^5 calls went from 0.050 s to 0.022 s (the reference 0.034 s).
- `Compile`: the scalar form lowers to `OP_POWMOD_I`, and the array form `PowerMod[v, e, m]` delegates to `ndint_powermod` through `A_NDFN`. That takes it out of the coverage BASELINE.
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

> **Packed arrays.** `PowerMod[list, e, m]` over an `int64` buffer is
> square-and-multiply on the machine words, with no GMP allocation per
> element. Taken for a non-negative `e` and a modulus small enough that the
> squaring step cannot overflow; a negative exponent is a modular inverse,
> which exists only for some elements, and takes the ordinary path.

## MultiplicativeOrder

Gives the multiplicative order of `k` modulo `n` — the smallest positive
integer `m` such that $k^m \equiv 1 \pmod{n}$.
- `MultiplicativeOrder[k, n]`: order of `k` modulo `n`.
- `MultiplicativeOrder[k, n, {r1, r2, ...}]`: smallest positive `m` with
  $k^m \equiv r_i \pmod{n}$ for some `i` (a multi-target discrete log).

The order is the size of the cyclic subgroup generated by `k` in the units
modulo `n`. Reversing the map — given `k` and a target `b`, find `x` with
$k^x \equiv b \pmod n$ — is the **discrete logarithm problem**, and the
three-argument form solves exactly that. Computing the order is easy (it divides
`EulerPhi[n]`); the discrete logarithm it inverts is believed hard, and that
presumed hardness underlies Diffie–Hellman key exchange and the ElGamal and DSA
schemes. See the *Handbook of Applied Cryptography* (Menezes, van Oorschot,
Vanstone, 1996), §3.6, for the standard algorithms (baby-step giant-step,
Pohlig–Hellman, index calculus).

**Features**:
- `Protected`.
- All arithmetic uses GMP `mpz_t`, so `k`, `n`, and any `r_i` may be
  arbitrary-precision bignums.
- `k` (and each `r_i`) is reduced modulo `n` before the search, so negative
  or out-of-range `k` works.
- Returns unevaluated when `gcd(k, n) != 1`, when `n <= 0` (as in
  the reference; a negative `n` used to be read as `|n|`), or when no power
  of `k` lands in the residue set (3-arg form).
- Machine path: for a machine `k` and `1 < n < 2^32` the order is found
  entirely in `uint64`, with `phi(n)` and its factorisation from the int64
  helpers the NDArray kernels use. The `Do`-loop of 10^4 calls went from
  0.15 s to 0.012 s (the reference 0.056 s), and agrees with the reference on a
  3000-entry table.
- The 2-arg form factors `phi(n)` and successively strips prime factors
  from `phi(n)` whose corresponding exponent still maps `k` to 1 — so
  the work is dominated by factoring `phi(n)`, not by walking the orbit.
- The 3-arg form is a discrete logarithm in the cyclic group generated by
  `k`, of order `d`, solved as PARI's `znlog` does
  (`src/numbertheory/nt_dlog.c`):
  - **Pohlig–Hellman** over each prime power `q^e` of `d`, one base-`q` digit
    at a time; the residues are recombined by CRT.
  - **Baby-step giant-step** for each prime `q` up to 2^36 (2^32 on the GMP
    back-end). One table per `q` is shared by every digit and every target.
  - **Pollard ρ** above that: a 32-multiplier r-adding walk with Brent cycle
    detection. It polls `TimeConstrained` and returns cleanly.
  - Arithmetic is `uint64` Montgomery for `n < 2^63` (a 128-bit remainder for
    even `n`) and GMP above.
  - Every answer is checked with one `PowerMod` before it is returned. `n` is
    factored by FactorInteger's engine, with a primality test first for the
    common prime modulus.
  - There is no size cap: the cost is about `sqrt` of the largest prime of
    `d`. A target outside `<k>` stays unevaluated, including in a non-cyclic
    `(Z/nZ)^*` such as `n = 15`. There, ρ gives up after four bounded walks;
    the chance of a wrong "no" is far below 2^-100.
  - Timing (`benchmarks/104-discrete-log`, ms, Mathilda / PARI): 24-bit
    subgroup ×40 3.8 / 27 (was 15 316); 32-bit ×10 13 / 79 (was
    unevaluated); 131-bit smooth `p` ×40 32 / 54 (was unevaluated).
- Non-integer numeric inputs (`Real`, `Complex`, `Rational`) and symbolic
  arguments flow through unevaluated with no diagnostic.
- Diagnostic: `MultiplicativeOrder::argtu` for one argument and `::argt`
  for other counts, as in the reference.

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

## JacobiSymbol

Gives the Jacobi symbol $(n/m)$.
- `JacobiSymbol[n, m]`: the Jacobi symbol of `n` with respect to `m`,
  returned as `-1`, `0`, or `1`.

**Features**:
- `Protected`, `Listable` — threads element-wise over lists and arrays.
- For prime `m` the Jacobi symbol reduces to the Legendre symbol, equal to
  `±1` according to whether `n` is a quadratic residue modulo `m`, and `0`
  when `m` divides `n`.
- Following the reference, this is the full Kronecker-symbol
  generalisation: the second argument `m` may be even or non-positive, and
  the first argument `n` may be negative.
- Computed with GMP `mpz_kronecker` in $O((\log m)^2)$ time, so `n` and `m`
  may be arbitrary-precision bignums.
- Non-integer numeric inputs and symbolic arguments flow through
  unevaluated with no diagnostic.
- Diagnostic: `JacobiSymbol::argr` for one argument and `::argrx` otherwise, as in the reference.
- Machine surfaces: two machine integers use an int64 Kronecker symbol (`ndint_kronecker_i64`, Cohen Alg. 1.4.10). That same function serves the packed / `NDArray` kernel (on `AWARE` and `INT64_OK`) and `Compile` (`OP_JACOBI_I`, at scalar and array shapes), so all of them agree with `mpz_kronecker`. `JacobiSymbol[Range[10^6], 1009]` takes 63 ms (the reference 0.29 s).

```mathematica
In[1]:= JacobiSymbol[10, 5]
Out[1]= 0

In[2]:= Table[JacobiSymbol[n, m], {n, 0, 10}, {m, 1, n, 2}]
Out[2]= {{}, {1}, {1}, {1, 0}, {1, 1}, {1, -1, 0}, {1, 0, 1},
         {1, 1, -1, 0}, {1, -1, -1, 1}, {1, 0, 1, 1, 0}, {1, 1, 0, -1, 1}}

In[3]:= JacobiSymbol[10^10 + 1, Prime[1000]]
Out[3]= 1

In[4]:= JacobiSymbol[7, 6]
Out[4]= 1

In[5]:= JacobiSymbol[{2, 3, 5, 7, 11}, 3]
Out[5]= {-1, 0, -1, 1, -1}

In[6]:= JacobiSymbol[-3, {1, 3, 5, 7}]
Out[6]= {1, 0, -1, 1}
```

## DirichletCharacter

- `DirichletCharacter[k, j, n]` — the Dirichlet character χ_{k,j}(n), modulus
  k, index j = 1..EulerPhi[k] (j = 1 principal).

Implemented in `src/special_functions/dirichlet_character.c` on the shared
engine `src/special_functions/dirichlet.{c,h}`, registered via
`dirichletcharacter_init()`. Attributes `Protected` (as in the reference); it
threads over lists in every argument itself.

- **Labelling (the reference's, verified on 148,138 values for k ≤ 100 and large
  prime powers / 2-adic / composite moduli).** Factor k with primes ascending;
  (Z/kZ)* splits into cyclic components — one per odd p^a, generated by the
  smallest primitive root mod p^a; one for 4 (generator 3); two for 2^a, a ≥ 3:
  ⟨−1⟩ then ⟨5⟩. j − 1 is the mixed-radix number over the components, the first
  component most significant, and χ(g_c) = exp(2πi d_c/ord_c). So
  `DirichletCharacter[k, j, PrimitiveRoot[k]] = Exp[2 Pi I (j-1)/EulerPhi[k]]`
  for prime k, and `DirichletCharacter[45, 8, n] =
  DirichletCharacter[9, 2, n] DirichletCharacter[5, 4, n]`.
- **Values** are exact: 0 when GCD[n, k] > 1, otherwise a product of the
  component values in the reference's normal form — the components that are units
  (±1, ±I) multiply into a unit, the rest combine into one E^(I π r),
  r ∈ (−1, 1), folded back into the unit when it lands on one:
  `DirichletCharacter[45, 8, 7] = -I E^(-2 I Pi/3)`,
  `DirichletCharacter[32, 10, 3] = -E^(3 I Pi/4)`,
  `DirichletCharacter[20!, 300, 23] = E^(I Pi/9)`.
- **Engine.** Per (k, j) a table (cached, immutable once published, safe from
  threaded Compile loops): components, digits, order, parity, conductor and the
  primitive χ*, the Euler-factor primes, and — on demand, for f ≤ 2^22 — dense
  χ* tables (exact exponents, machine values, root number) built by walking the
  component generators. Single values use Pohlig–Hellman discrete logs with
  baby-step giant-step in each prime-order subgroup (k < 2^63; a component
  order whose largest prime factor exceeds 2^44 leaves the call unevaluated).
- **Messages.** `argrx`, `intp` (position 1 or 2), `invl` (j > EulerPhi[k]),
  `int` (a non-integer n).
- **Surfaces.** An int64 packed / `NDArray` n: a real character returns an
  int64 buffer of 0, ±1; a complex one a List of exact roots of unity (memoised
  per residue). Exact-valued, so no `Compile[]` lowering (compile-coverage
  `EXEMPT`).

```mathematica
In[1]:= Table[DirichletCharacter[5, 2, n], {n, 1, 5}]
Out[1]= {1, I, -I, -1, 0}

In[2]:= Table[DirichletCharacter[8, j, n], {j, 1, 4}, {n, 1, 8}]
Out[2]= {{1, 0, 1, 0, 1, 0, 1, 0}, {1, 0, -1, 0, -1, 0, 1, 0}, {1, 0, -1, 0, 1, 0, -1, 0}, {1, 0, 1, 0, -1, 0, -1, 0}}

In[3]:= DirichletCharacter[45, 8, 7]
Out[3]= -I E^((-2/3*I) Pi)
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
  primitive root of $p^k$. This matches the reference's convention so that,
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
- Returns `{}` if `|n|` is not 2, 4, an odd prime power, or twice an odd
  prime power. The sign of `n` is ignored, as in the reference:
  `PrimitiveRootList[-7]` → `{3, 5}` (it used to give `{}`).
- Enumerates the $\varphi(\varphi(n))$ primitive roots as $g^i \bmod n$
  for $i \in \{1, \ldots, \varphi(n)\}$ with $\gcd(i, \varphi(n)) = 1$,
  where $g$ is the smallest primitive root of `n`. The list is sorted
  ascending.
- Falls back to unevaluated when $\varphi(n)$ exceeds `unsigned long`,
  since the enumeration cannot be represented.
- Non-integer numeric inputs (e.g. `11.0`, `11 + I`) flow through
  unevaluated with no diagnostic, matching the reference.
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
- Every well-formed call returns `True` or `False`. Non-integer and non-Gaussian
  inputs (symbols, `Sqrt[2]`, `Exp[2 Pi I/3]`, strings, etc.) give `False`.
- A **malformed** call stays unevaluated with the reference's message: `PrimeQ[]` →
  `PrimeQ::argx`, `PrimeQ[1, 2]` → `::nonopt`, an unknown option → `::optx`, and
  `GaussianIntegers -> Automatic` → `::opttf`. (Before v0.472 these gave `False`.)
- Packed `int64` buffers answer in one pass with a packed `True`/`False` List
  (`ndint_primeq`). Dense small values (max ≤ 2²⁶ with at least max/64 elements)
  use a byte sieve; anything else uses the scalar's deterministic Miller–Rabin.
  `PrimeQ[Range[10^6]]` went from 0.53 s to 0.004 s (the reference 0.32 s).
  `PrimeQ[n]` compiles to `OP_PRIMEQ_I`; the direct array form is in the
  coverage BASELINE with the other boolean predicates.
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


## PrimePowerQ

- `PrimePowerQ[expr]`: `True` if `expr` is a power of a prime number, `False` otherwise.
- `PrimePowerQ[m + n I]`: works over the Gaussian integers automatically. It is `True` when the argument is a unit times a power of a single Gaussian prime.
- `PrimePowerQ[n, GaussianIntegers -> True]`: tests whether `n` is a Gaussian prime power. The default is `GaussianIntegers -> Automatic`.

**Semantics** (each checked against the reference):
- A prime power is a prime or a positive integer power of a prime, up to sign. `PrimePowerQ[-8]` is `True`; `0` and `±1` are not prime powers.
- A rational or Gaussian rational is a prime power when its reciprocal is one. `1/25`, `-1/8`, `1/(2+I)` and `I/2` (that is, `1/(-2 I)`) give `True`. `2/3` and `3/(2+I)` give `False`.
- **Gaussian integers.**
  - On an axis, `z = unit·m`. Then `z` is a Gaussian prime power when `m = 2^k`, or when `m = q^j` with `q ≡ 3 (mod 4)`.
  - Off the axes, the norm must be `p^k`. With `p = 2` the answer is always `True`. With `p ≡ 1 (mod 4)` it is `True` only when `p` does not divide both parts.
  - So `3 + 4 I = -I (2 + I)^2` is `True`, but `5 (2 + I)` is `False`.
  - With `GaussianIntegers -> True`, `5` and `25` give `False`, while `2`, `3`, `8` and `9` give `True`.
- `GaussianIntegers -> False` with a complex argument gives `PrimePowerQ::par` and stays unevaluated.
- A well-formed `GaussianIntegers :> v` answers `False`, as the reference does.
- Any other exact expression (`x`, `Pi`, `"s"`) gives `False`.
- An inexact number gives `PrimePowerQ::exact` and stays unevaluated.
- Argument and option errors (`argt`, `nonopt`, `optx`, `opttfa`) also stay unevaluated.
- Every message routes through the funnel, so `Quiet` and `Check` see them.
- Attributes: `Listable`, `Protected`. Options: `{GaussianIntegers -> Automatic}`.

**Algorithm** (`mth_prime_power_base_u64` / `_mpz` in `src/facint.c`):
- **Machine integers** need no GMP:
  1. Even numbers: a power of 2 if exactly one bit is set.
  2. Trial division by the odd primes below 64. A hit strips that prime and checks for 1.
  3. One deterministic Miller–Rabin test.
  4. If composite, try exact k-th roots for k ∈ {2, 3, 5, 7}. Once no prime below 64 divides `n`, the base is at least 67, so `k ≤ 10`.
- **Bignums:**
  1. Trial division by the primes below 2^10. A hit costs one `mpz_remove`.
  2. Reduction to the primitive root: FLINT `fmpz_is_perfect_power`, or GMP `mpz_perfect_power_p` + `mpz_root` without FLINT.
  3. One `PrimeQ` test (BPSW + 25 Miller–Rabin rounds) on that root, the smallest candidate.
- **Buffers.** A packed or visible `int64` buffer is read directly (`ndint_prime_power_predicate`, `src/ndinteger.c`) and answers with a packed boolean List.
- **Compile.** `Compile[{{n, _Integer}}, PrimePowerQ[n]]` lowers to `OP_PPQ_I`, which calls the same helper.

```mathematica
In[1]:= {PrimePowerQ[8], PrimePowerQ[6], PrimePowerQ[2 + I], PrimePowerQ[1/25]}
Out[1]= {True, False, True, True}

In[2]:= PrimePowerQ[5, GaussianIntegers -> True]
Out[2]= False

In[3]:= PrimePowerQ[10^300 + 3]
Out[3]= False

In[4]:= Select[Range[100], PrimePowerQ[#] && !PrimeQ[#] &]
Out[4]= {4, 8, 9, 16, 25, 27, 32, 49, 64, 81}

In[5]:= PrimePowerQ[3.2912]
PrimePowerQ::exact: Argument 3.2912 in PrimePowerQ is not an exact number.
Out[5]= PrimePowerQ[3.2912]
```

## Prime

- `Prime[n]`: the `n`-th prime $p_n$ (`Prime[1] = 2`, `Prime[2] = 3`, ...).

**Features**:
- `Listable`, `Protected`.
- Small `n` is read directly from a sieve table of the primes below $10^6$.
- Large `n` (v0.472):
  - $p_n$ is estimated by Newton on Riemann's $R(x) \approx \mathrm{li}(x) - \mathrm{li}(\sqrt x)/2 - \mathrm{li}(x^{1/3})/3$, which is within $O(\sqrt{p_n})$.
  - **One** exact `PrimePi` is taken there; a second only if the estimate is far off.
  - The short gap to $p_n$ is closed with a segmented sieve (`primecount_nth_from`) instead of stepping `NextPrime` over up to $10^6$ primes.
  - Timings: `Prime[10^9]` went from 0.24 s to 0.018 s, `Prime[10^11]` takes 0.41 s and `Prime[10^12]` 2.4 s. The reference takes 1.4 ms for `Prime[10^9]` (from tables), but 3.0 s for `Prime[10^11]` and 15.7 s for `Prime[10^12]`.
  - Exact for $p_n \le 5 \times 10^{13}$ ($n \approx 1.6 \times 10^{12}$); beyond that the call is left unevaluated.
- A non-positive-integer argument emits `Prime::intpp`; a wrong argument count
  emits `Prime::argx`; both leave the call unevaluated.

```mathematica
In[1]:= Prime[100]
Out[1]= 541

In[2]:= Prime[{1, 3, 4, 10}]
Out[2]= {2, 5, 7, 29}

In[3]:= Prime[10^10]
Out[3]= 252097800623
```

> **Packed arrays.** `Prime[list]` over an `int64` buffer is one **sieve** to
> the largest requested index followed by a gather, rather than a separate
> prime count per element — so the whole array costs about what its largest
> single index used to.

## PrimePi

- `PrimePi[x]`: the number of primes less than or equal to `x`.
- `PrimePi[x, Method -> m]`: count with a specific algorithm.

**Features**:
- `Listable`, `Protected`.  Default `Method -> Automatic`.
- Several recognised prime-counting algorithms are available via `Method`: `"Sieve"` (segmented sieve of Eratosthenes), `"Legendre"`, `"Meissel"`, `"Lehmer"`, `"LMO"` (Lagarias–Miller–Odlyzko, segmented special leaves), `"DelegliseRivat"` and `"LucyHedgehog"`. They are independent implementations, and the tests cross-check them.
- **`"DelegliseRivat"` is the `Automatic` method** (v0.472), above the $10^6$ table. It is tuned:
  - $y = \alpha x^{1/3}$ with $\alpha = 0.0008 \ln^3 x$, which shrinks the sieve range $z = x/y$ by $\alpha$.
  - Trivial leaves ($v < p_b$) are counted, not visited.
  - Easy leaves ($v < p_{b-1}^2$) are answered from the $\pi$ table, and above $\sqrt z$ only prime $m$ are walked.
  - Hard leaves (levels with $p_{b-1} \le \sqrt z$ only) use a bitset sieve with monotone popcount sweeps in place of a Fenwick tree.
  - Leaf quotients come from a precomputed $1/m$ with an exact correction.
  - Fixing the "easy" formula for general $y$ ($\varphi = 1$ when $v < p_b$) is what made $\alpha > 1$ correct.
  - Timings, against the untuned $y = x^{1/3}$ version (one core): $10^{11}$ 0.54 → 0.04 s, $10^{12}$ 4.4 → 0.20 s, $10^{13}$ 40 → 1.06 s, $5\cdot10^{13}$ 3.4 s.
  - The reference answers `PrimePi[10^12]` in 0.088 s in a fresh kernel, and near $10^{14}$ to $10^{15}$ in under half a second, which points to precomputed checkpoints. Mathilda has no tables and stops at $5 \times 10^{13}$.
- `x` may be any real: a non-integer is floored (an exact numeric such as `Pi` or `E^10` exactly). `PrimePi[Infinity]` is `Infinity` and `PrimePi[-Infinity]` is 0. Complex and symbolic `x` stay unevaluated.
- Messages (v0.472, as in the reference): `PrimePi::argx` for no argument, `::nonopt` for a non-rule extra argument, `::optx` for an unknown option. An unrecognised `Method` value leaves the call unevaluated silently (it used to draw `PrimePi::method`).
- Options: `{Method -> Automatic, ProgressReporting :> $ProgressReporting}`. `ProgressReporting` is accepted and has no effect.

```mathematica
In[1]:= PrimePi[10^9]
Out[1]= 50847534

In[2]:= PrimePi[10^9, Method -> "LMO"]
Out[2]= 50847534

In[3]:= PrimePi[{10, 100}]
Out[3]= {4, 25}
```

## NextPrime

- `NextPrime[x]`: Smallest prime above `x`.
- `NextPrime[x, k]`: $k$-th next prime above `x` (or previous if $k$ is negative).

**Features**:
- `Protected`, **not** `Listable` (as in the reference), but a list in either argument is threaded explicitly: `NextPrime[{4, 10}, 2]` → `{7, 13}`, `NextPrime[4, {1, 2}]` → `{5, 7}`.
- **Negative primes count**, as in the reference (v0.472). The primes are ordered …, −7, −5, −3, −2, 2, 3, 5, 7, …, so `NextPrime[-10]` → `-7`, `NextPrime[2, -1]` → `-2` and `NextPrime[3, -2]` → `-2`. Before v0.472, `NextPrime[-10]` was `2` and `NextPrime[2, -1]` stayed unevaluated.
- `NextPrime[x, 0]` is `x` when `x` is a (signed) prime, and unevaluated otherwise.
- A real or rational `x` starts from `Floor[x]` going up and `Ceiling[x]` going down: `NextPrime[2.5, -1]` → `2`. Rationals of any size are exact (`NextPrime[(355/113)^720 + 123948304759283749084]` needs no `Floor`), and any other real numeric quantity is rounded by `Floor`/`Ceiling` itself: `NextPrime[Pi]` → `5`, `NextPrime[Sqrt[2]]` → `2` (v0.525). A symbolic, complex or infinite `x` stays unevaluated.
- A non-integer `k` leaves the call unevaluated. A wrong argument count draws `NextPrime::argt`.
- Checked against the reference on the full grid $x \in [-15, 15]$, $k \in [-3, 3]$, plus non-integer starts.

## FactorInteger

- `FactorInteger[n]`: the prime factorisation `{{p1, e1}, {p2, e2}, ...}`, primes ascending.
- `FactorInteger[n, k]`: a partial factorisation with **at most `k` entries**, grouped the way
  the reference groups them (all expectations below were checked against it):
  1. Primes below 2^16 are removed by trial division; what is left is the cofactor C.
  2. If C = 1, the k - 1 largest prime powers are kept and the rest merge into one entry:
     `FactorInteger[360, 2]` is `{{5, 1}, {72, 1}}`.
  3. If there are more trial primes than k - 1 and C > 1, the smallest trial primes merge, the k - 2
     largest are kept, and C stays unsplit: `FactorInteger[3 5 7 1000003 1000033, 3]` is
     `{{7, 1}, {15, 1}, {1000036000099, 1}}`.
  4. Otherwise C is split into the room left, smallest factors first, and any overshoot merges.
  5. A negative sign goes onto the merged entry, else onto an unsplit C, else onto the smallest entry,
     else onto the largest prime power (written with exponent 1): `FactorInteger[-36, 2]` is
     `{{-9, 1}, {2, 2}}`. For a rational, the sign prefers the denominator part.
  6. Entries are sorted by base.

  `FactorInteger[n, 1]` is `{{n, 1}}`; `k < 1` leaves the call unevaluated. The reference splits the
  cofactor in the order its own search happens to find factors, so when C has three or more prime
  factors and only some fit, the split can differ from the reference's (Mathilda's is smallest-first).
  Gaussian partial factorisation with `k >= 2` keeps the full answer: the reference's output there is
  inconsistent and is not reproduced.
- `FactorInteger[n, Automatic]`, `FactorInteger[n, All]`, `FactorInteger[n, Infinity]`: the same as
  `FactorInteger[n]`.
- `FactorInteger[n, GaussianIntegers -> True]`: factorisation over Z[i]. A Gaussian integer or Gaussian
  rational argument (`3 + 4 I`, `(3 + 4 I)/(1 + I)`) is factored over Z[i] automatically. The unit comes
  first, when it is not 1. Primes are given in the first quadrant (Re > 0, Im >= 0), sorted by real part
  and then imaginary part. Every output matches the reference.
- `FactorInteger[n, Method -> m]`: the algorithm is chosen by `m`, which is a method name, `Automatic`,
  or `{name, "Sub" -> v, ...}`. Every method is also a builtin of its own,
  `FactorInteger`name[n, "Sub" -> v, ...]`, taking the same `k` and `GaussianIntegers` arguments.
- Rationals: the denominator's primes get negative exponents. Negative `n` starts with `{-1, 1}`;
  `0` gives `{{0, 1}}`, `1` gives `{{1, 1}}`.
- Messages, all routed through the `Quiet`/`Check` funnel:
  - `FactorInteger::nofac`: a composite was left unfactored, once per such composite;
  - `FactorInteger::method`: unknown method;
  - `FactorInteger::mopt`: bad suboption;
  - `FactorInteger::optx`: unknown option;
  - `FactorInteger::exact`: inexact or symbolic argument, as in the reference;
  - `FactorInteger::argt`: no arguments; `FactorInteger::nonopt`: a third non-rule argument.

The implementation is in `src/numbertheory/factorinteger/` (one file per algorithm; see `fi_internal.h`).
Every randomised choice (rho constants, ECM sigmas, p+1 seeds, QS polynomials) is drawn from a generator
seeded by `n`, so the answer, the time taken and any `nofac` depend on `n` alone.

| Method (`FactorInteger`name`) | Algorithm | Suboptions |
|---|---|---|
| `Automatic` | Below 2^64 the word engine handles everything: Montgomery arithmetic, a deterministic Miller–Rabin test, micro-ECM with Suyama curves and a baby-step/giant-step stage 2, and rho/SQUFOF as backstops. Above 2^64 the tiers are: trial division to 2^16 → Pollard rho (two-limb Montgomery below 2^126) → two-limb micro-ECM (below 2^126) → a Fermat probe → one p-1 → one p+1 (three seeds from 70 digits) → ECM levels 15, 20, ... digits, switching to the quadratic sieve once the ECM effort reaches 2/7 of the composite's digit count, for composites up to 100 digits. After a split, both parts restart at the tier that succeeded. | — |
| `"TrialDivision"` | Primes are tested in word-sized products, one bignum remainder per 2–4 primes; a segmented sieve takes over above 2^20. The default bound is sqrt(n) when that is below 10^7, else 10^7. | `"Bound"` (≤ 2^40) |
| `"PollardRho"` | Brent's variant with batched gcds and backtracking; Montgomery arithmetic on one and two words, GMP above that. | `"Iterations"`, `"Starts"` |
| `"Fermat"` | Fermat's method with quadratic-residue filters modulo nine small moduli, applied 64 values of a at a time through precomputed bit masks (~0.5 ns per value; budget 2^26 values by default). | `"Iterations"` |
| `"Lehman"` | Lehman's O(n^(1/3)) method; complete below 2^42, budgeted above (2^28 values of a by default, set-up included). On bignums the multipliers are interleaved in rounds of a level c growing 2× per round: multiplier k ≥ 2 is scanned to the fraction min(1, c/k) of its window and k = 1 (plain Fermat) to min(1, 8c), so near-square factors cost ~1.1× what plain Fermat pays while factors with p/q near a small rational u/v (multiplier k = uv) are still found at any size. The sweep order does not depend on the budget, so a larger `"Iterations"` only scans further (it never finds less). | `"Iterations"` |
| `"ShanksSquareForms"` (alias `"SQUFOF"`) | SQUFOF over the 16 Gower–Wagstaff multipliers, in word arithmetic for kN < 2^124. | `"Iterations"` |
| `"CFRAC"` | Morrison–Brillhart: Knuth–Schroeppel multiplier, early-abort trial division, single large primes, GF(2) elimination. Relation collection is pipelined across threads: one thread extends the continued fraction into the next batch while every thread trial-divides the current one, and the hits are replayed in term order. The relations, restarts and divisor are therefore those of the serial method for any thread count. | `"Iterations"` |
| `"Dixon"` | Random squares drawn near sqrt(kN), so the residues are square-root sized, with the same relation machinery as CFRAC. The draws are taken in batches; the candidates are built and trial-divided in parallel, and the hits are replayed in draw order. Results and the random stream are those of the serial method. | `"Iterations"` |
| `"QuadraticSieve"` (aliases `"QS"`, `"SIQS"`, `"MPQS"`; the heads `FactorInteger`QS` and `FactorInteger`SIQS` exist, `"MPQS"` is a `Method` string only) | Native self-initialising QS: Gray-code polynomial switching, 32 KB block sieve with bucket sieving for primes ≥ 2^15, batched candidate trial division. Single large primes, plus double large primes from 72 digits, combined by cycle-finding in the large-prime graph. Polynomial families are sieved in parallel from 45 digits (one A per thread). Linear algebra is block Lanczos (Montgomery 1995) above 1000 columns, dense elimination below. | — |
| `"PollardPMinus1"` (alias `"PollardP-1"`) | GMP-ECM p-1, with a B1 back-off when every prime is caught at once; a word walk with prime-by-prime backtracking below 2^64. | `"B1"` (10^6), `"B2"` |
| `"WilliamsPPlus1"` (alias `"WilliamsP+1"`) | GMP-ECM p+1 with seeds 2/7, 6/5, 3, and the same back-off and word walk. The seeds run in parallel, one per thread, and the first success in seed order is taken. | `"B1"` (10^6), `"B2"` |
| `"ECM"` | GMP-ECM, climbing its optimal levels (B1 = 2000, 11000, 5e4, ...) with the expected curve count at each. Curves run in parallel, one per thread, for n ≥ 2^160. They use the batch parametrisation with sigmas drawn deterministically, and successes are taken in curve order, so the result does not depend on scheduling. A user `"Sigma"` is a Suyama parameter. There is a B1 back-off on gcd = n. | `"B1"`, `"B2"`, `"Curves"`, `"Sigma"` |
| `"BlakeRationalBaseDescent"` | Blake's rational base descent, for a factor close to c (a/b)^k: six variants (exact descent, the reference descent, all-hits, lattice, two-pass, Coppersmith), threaded base search. See [below](#blakerationalbasedescent). | `"Base"`, `"Variant"`, `"Order"`, `"Height"`, `"Bases"`, `"MaxCoefficient"`, `"Multiplier"`, `"Tolerance"`, `"Beta"`, `"Multiplicity"`, `"Dimension"`, `"BaseRange"`, `"Coefficient"`, `"Exponent"` |

Fallbacks: builds without GMP-ECM use FLINT's ECM and p+1 (and a native p-1); builds without either use
native mpz ECM, p-1 and p+1. The quadratic sieve and every other method are native, with no external
dependency.

Threads: the quadratic sieve, ECM, p+1, CFRAC, Dixon and BlakeRationalBaseDescent use every online CPU (up to 16) in builds with
`MATHILDA_THREADS`. `TimeConstrained` budgets are elapsed time (v0.463; they were process CPU time, so a
16-thread region used a budget up 16 times too fast). An abort is taken on the main thread as soon as the
current region ends. Long units of work stop early when the deadline passes: ECM curves and p±1 seeds through
GMP-ECM's `stop_asap` hook, QS families every 16 polynomials, and BlakeRationalBaseDescent per candidate. A
CFRAC/Dixon region lasts only a few milliseconds.

### BlakeRationalBaseDescent

`FactorInteger`BlakeRationalBaseDescent[n]` (the same as `FactorInteger[n, Method -> "BlakeRationalBaseDescent"]`)
is S. Blake's *rational base descent* (RBD), from two papers: "Rational Base Descent: A Deterministic
Algorithm for Factoring Structured Semiprimes" and "Lattice-Augmented Rational Base Descent". It is a
special-purpose method. It succeeds when the **larger** factor P of `n = P M` is close to a power of a
rational number, `P ≈ c (a/b)^k`, and is no threat to unstructured (e.g. RSA) moduli.

**Mechanism.** Starting from `Q = n`, repeat `Q <- floor(b Q / a)`. Each step divides by the base. After
`k` steps `Q` has divided `P` down to about `c`, so `Q` lies within a small constant of the integer `c M`,
and `gcd(Q - j, n)` for a small `j` exposes `M` (or `M gcd(c, P)`). The paper's theorem: if
`(a/b)^k > M` and `|P - floor(c (a/b)^k)| < (a/b)^k / M`, a factor is found at step `k`. One descent tests
every exponent and every coefficient at once, in time polynomial in `log n` for a fixed base.

The base is supplied with `"Base"` or found by enumeration. The bases that succeed form a comb of isolated
teeth of relative width ~1/(k M) around the representations of P (plus harmonic teeth at `a/d`, `d | a`).
No signal points towards a tooth from a distance, so bases are enumerated, every rational once and
skipping perfect powers `(x/y)^m`. By default they come in ascending `a / log2(a/b)`. A descent with
base `r` costs about `log n / log r` steps, and a simple base (small `a`) is a priori likelier, so this
is weighted shortest-job-first. On random targets it reaches the right base 2–2.5× sooner on average
than breadth-first Stern–Brocot, under both a uniform and a log-height prior. For an integer
base the descent truncates the base-`a` digits of `n`, and truncations inherit algebraic divisors. This
*truncation divisibility* is what recovers the cyclotomic factors of Cunningham numbers and repunits.

| Variant (`"Variant"` -> name; also the head `FactorInteger`BlakeRationalBaseDescent<name>`) | Algorithm |
|---|---|
| `"Exact"` (default) | The scaled exact descent. The recurrence runs on `S n`, with `S = 2^s ≥ 2a/(a-b)`, so the target always lies in a base-independent window of four consecutive values. Nothing is tested above `sqrt(cmax n)`, since no target `c M` can lie there; that depth is reached by one exact division. The window products are folded into an accumulator mod `n`, with one gcd per block of steps; a block is replayed candidate by candidate only when its gcd is non-trivial. |
| `"Classic"` | The paper's reference algorithm R1–R4. It tests `gcd(Q - j, n)` for `|j| ≤ ceil(2 + a/(a-b))` at every depth, searching before the termination test. |
| `"AllHits"` | The reference descent run to the end. Every distinct non-trivial gcd is recorded, and the record is refined into pairwise-coprime divisors of `n`; later calls on the cofactors split from that cache. Composite pieces are re-run with AllHits, which recovers nested algebraic divisors. |
| `"Lattice"` | Lattice-augmented descent. Howgrave-Graham's approximate-divisor detector (LLL on `N^(t-i) (Q+x)^i`, `x^j (Q+x)^t`) replaces the gcd window. It finds `c M` from `Q` within `X ≤ n^(β²)` when `M ≥ n^β`, widening the per-step tolerance from O(a/(a-b)) to `X`. The detector runs only in the depth band `X ≤ Q ≤ sqrt(cmax n) + X`, in parallel over the depths. Integer roots are found p-adically on the square-free part of each short polynomial. |
| `"TwoPass"` | The paper's practical algorithm. Pass 1 is Exact over the bases. Pass 2 is Lattice over a coarse grid of `"BaseRange" -> {lo, hi}`, with relative spacing `X/(2 c k n^β)` and each grid point realised by the simplest rational near it; without a range, pass 2 uses the Stern–Brocot bases of height 16 (or `"Height"`). |
| `"Coppersmith"` | A known base, coefficient and exponent: one detector call at `floor(c (a/b)^k)`. This is factoring with high bits known. It needs `"Base"`, `"Coefficient"` and `"Exponent"` (else `FactorInteger::rbdcop`). |

Suboptions (defaults in parentheses):
- `"Base"` -> `a/b` or a list of bases. A base is reduced, and inverted if below 1; `1` is rejected with
  `FactorInteger::mopt`. (Default: enumerate.)
- `"Order"` -> `"Weighted"` (default; ascending `a / log2(a/b)`), `"HeightDoubling"` (Stern–Brocot
  breadth-first in height bands `a ≤ 16`, `16 < a ≤ 32`, …), `"SternBrocot"` (breadth-first over the
  whole tree, coarse to fine; the order before v0.526) or `"Sum"` (the papers' `a + b` order).
  The Weighted order is streamed in doubling key bands, each a pruned Stern–Brocot walk, so a large
  `"Height"` does not materialise the tree. The TwoPass lattice pass keeps breadth-first order.
- `"Height"` -> `H`: largest numerator enumerated (1000). `"Bases"` -> `B`: enumeration budget (200000).
- `"MaxCoefficient"` -> `cmax` (2^64): sets the deferred-testing ceiling `sqrt(cmax n)`; `Infinity`
  tests every depth.
- `"Multiplier"` -> `m`: descend from `m n`, with gcds still taken with `n` (the paper's Open Question 2).
- `"Tolerance"` -> `X`, `"Beta"` -> `β` (1/2), `"Dimension"` -> `d` (12), `"Multiplicity"` -> `t`: the
  lattice detector. Given a tolerance, the least dimension whose LLL guarantee covers it is chosen.
  Otherwise `X` is the largest tolerance the guarantee covers at dimension `d`, less two bits.
- `"BaseRange"` -> `{lo, hi}` (TwoPass grid), `"Coefficient"` -> `c`, `"Exponent"` -> `k` (Coppersmith).

Behaviour:
- The suboptions describe the input. A composite cofactor that the configured search cannot split gets a
  light default search (Exact, height 64 in the Weighted order, at most 2000 bases); AllHits re-runs AllHits
  instead. Anything still composite is returned with `FactorInteger::nofac`.
- Threads: a many-base search runs blocks of 2048 bases on every thread. A lattice base runs its band two
  depths per thread per region. The lowest-indexed success wins (a hit caps the work claim), so the answer
  does not depend on the thread count, and a `TimeConstrained` abort is taken within milliseconds.
- Lattice reduction is FLINT's `fmpz_lll` when available, else a native fraction-free integral LLL
  (de Weger / Cohen 2.6.7, δ = 0.99): exact, but slower.
- The Exact window accumulator is division-free (`fi_rbd_modacc.c`). It is a running product mod `n` in
  Montgomery-reduced form, on fixed-size limbs, with the modulus context built once per `n` and shared
  by every thread. Each reduction multiplies by a unit mod `n` (odd `n`), so `gcd(acc, n)` needs no
  conversion. Folding an `s`-limb factor costs one `n × s` product and `s` REDC passes, about the
  information bound of the method.
- Performance (v0.526, 16 threads): a 1630-bit `p q` with `p = NextPrime[⌊(355/113)^720⌋ + …]` splits
  in about 2.8 s, against 63 s with the breadth-first order and the division-based accumulator. A single
  known base takes about 5 ms.
- The core (`fi_rbd.h`; `fi_rbd_descent.c`, `fi_rbd_modacc.c`, `fi_rbd_bases.c`, `fi_rbd_lattice.c`,
  `fi_rbd_lll.c`) uses only GMP (and optionally FLINT), with no global state; the Mathilda glue is
  `fi_rbd.c`.

Performance, measured on 500 Stern–Brocot bases against a random semiprime, single thread:

| n (bits) | Exact (µs per base) | the v0.454 implementation | speedup |
|---|---|---|---|
| 256 | 61 | 2240 | 37× |
| 512 | 151 | 11400 | 75× |
| 1024 | 546 | 54800 | 100× |

The paper's Example 1 (1140 bits, base 22/7) takes 2.8 ms with Exact and 27 ms with Classic. The paper's
2048-bit challenge modulus is factored from scratch, with no base given, in 106 s on 16 threads.

```mathematica
In[1]:= FactorInteger`BlakeRationalBaseDescent[1001, "Base" -> 22/7]
Out[1]= {{7, 1}, {11, 1}, {13, 1}}

(* P = NextPrime[Floor[3 (22/7)^150]] times a 40-bit prime: the base is found by search *)
In[2]:= FactorInteger`BlakeRationalBaseDescent[742245361114170698080188160556986382238609541182118959142360194921632978579427957105431]
Out[2]= {{623347347967, 1}, {1190741187132579502456731851112421085210615890827643793135496647841845583593, 1}}

(* truncation divisibility: the algebraic factors of 2^1041 + 1 *)
In[3]:= Select[Quiet[FactorInteger`BlakeRationalBaseDescentAllHits[2^1041 + 1, "Base" -> 2]][[All, 1]], # < 10^120 &]
Out[3]= {3, 2083, 95562442332919646317117537304253622533190207882011713489066201641121786503686867002917439712921903606443}

(* a planted miss of 2^28 at depth 59: invisible to the gcd window, found by the lattice detector *)
In[4]:= First /@ FactorInteger`BlakeRationalBaseDescentLattice[195539142822724777401163665279552940381006352574199957151137, "Base" -> 22/7]
Out[4]= {296454907933147663952202071141, 659591518271709676788151702157}

In[5]:= First /@ FactorInteger`BlakeRationalBaseDescentCoppersmith[195539142822724777401163665279552940381006352574199957151137, "Base" -> 22/7, "Coefficient" -> 3, "Exponent" -> 59]
Out[5]= {296454907933147663952202071141, 659591518271709676788151702157}
```

**Features**:
- `Listable`, `Protected`; `Options[FactorInteger]` is `{GaussianIntegers -> False, Method -> Automatic}`.
- Perfect powers are reduced to a primitive base before any search, so `(2^61 - 1)^2` factors correctly.
- `EulerPhi` and `facint_factor_complete` decline (stay unevaluated or report incompleteness) rather than
  build an answer on a composite that survived the search.

```mathematica
In[1]:= FactorInteger[-12]
Out[1]= {{-1, 1}, {2, 2}, {3, 1}}

In[2]:= FactorInteger[3/4]
Out[2]= {{2, -2}, {3, 1}}

In[3]:= FactorInteger[360, 2]
Out[3]= {{5, 1}, {72, 1}}

In[4]:= FactorInteger[5, GaussianIntegers -> True]
Out[4]= {{-I, 1}, {1 + 2 I, 1}, {2 + I, 1}}

In[5]:= FactorInteger`ECM[2^128 + 1]
Out[5]= {{59649589127497217, 1}, {5704689200685129054721, 1}}

In[6]:= FactorInteger`QuadraticSieve[NextPrime[10^29] NextPrime[3 10^29]]
Out[6]= {{100000000000000000000000000319, 1}, {300000000000000000000000000007, 1}}
```

## Divisors

- `Divisors[n]`: the ascending list of positive integers that divide `n`.
- `Divisors[n, GaussianIntegers -> True]`: includes Gaussian-integer divisors.

**Options**:
- `GaussianIntegers -> False` (default): ordinary integer divisors. Set to
  `True` to compute divisors in the Gaussian integers `Z[i]`. A non-real
  Gaussian-integer input (e.g. `6 + 4 I`) auto-enables Gaussian mode.

**Features**:
- `Listable`, `Protected`.
- Machine integers and GMP bigints are handled uniformly; the result promotes to
  a big-integer list when needed.
- The sign of `n` is ignored (`Divisors[-12] == Divisors[12]`).
- Divisors are computed from the prime factorization (the divisor lattice), so
  cost scales with the number of divisors rather than `Sqrt[n]`.
- In Gaussian mode each divisor is the canonical first-quadrant representative of
  its associate class (`Re > 0`, `Im >= 0`), and the list is sorted by
  `(Re, Im)`. Rational primes are lifted to `Z[i]`: `2` ramifies as `1 + I`,
  primes `p ≡ 1 (mod 4)` split via a sum-of-two-squares (Cornacchia)
  decomposition, and primes `p ≡ 3 (mod 4)` stay inert.
- `Divisors[0]`, non-integer arguments, and calls whose divisor count would
  overflow (e.g. `Divisors[100!]`) are left unevaluated.
- Messages, matching the reference and routed through the `Quiet`/`Check`
  funnel: `Divisors::argx` (no arguments), `Divisors::nonopt` (a second
  non-rule argument), `Divisors::optx` (unknown option), `Divisors::opttf`
  (`GaussianIntegers` set to anything but `True`/`False`). When an option is
  given twice, the first setting wins.

```mathematica
In[1]:= Divisors[1729]
Out[1]= {1, 7, 13, 19, 91, 133, 247, 1729}

In[2]:= Divisors[6]
Out[2]= {1, 2, 3, 6}

In[3]:= Divisors[{605, 871, 824}]
Out[3]= {{1, 5, 11, 55, 121, 605}, {1, 13, 67, 871}, {1, 2, 4, 8, 103, 206, 412, 824}}

In[4]:= Divisors[6 + 4 I]
Out[4]= {1, 1 + I, 1 + 5 I, 2, 3 + 2 I, 6 + 4 I}

In[5]:= Divisors[2, GaussianIntegers -> True]
Out[5]= {1, 1 + I, 2}

In[6]:= Divisors[3, GaussianIntegers -> True]
Out[6]= {1, 3}
```

## DivisorSigma

- `DivisorSigma[k, n]`: the divisor function `sigma_k(n)`, the sum of the
  `k`-th powers of the (positive) divisors of `n`.
- `DivisorSigma[k, n, GaussianIntegers -> True]`: the sum over Gaussian-integer
  divisors.

**Options** (`Options[DivisorSigma]` is `{GaussianIntegers -> False}`):
- `GaussianIntegers -> False` (default): ordinary integer divisors. Set to
  `True` to sum over the divisors in `Z[i]`. A non-real Gaussian-integer input
  (e.g. `3 + I`) auto-enables Gaussian mode. `Rule` and `RuleDelayed` are both
  accepted; of repeated settings the first wins.

**Features**:
- `Listable`, `NHoldAll`, `Protected` (as in the reference, so
  `N[DivisorSigma[2, n]]` keeps the exact `2`).
- Three evaluation paths, chosen by `k` (`n` an ordinary integer):
  - **Machine `Integer` `k`** — exact, straight from the prime factorisation:
    an `Integer` for `k >= 0` and, through `sigma_-k(n) = sigma_k(n)/n^k`, a
    `Rational` in lowest terms for `k < 0`. A machine answer for
    `|n| <= 2^20` comes from the same int64 helper as the packed kernel and
    `Compile` (`ndint_divisor_sigma_i64`); otherwise GMP, by a Horner sum per
    prime.
  - **Machine `Real` `k`** — in `double`, summing `1 + q + ... + q^e`
    (`q = p^k`) by Horner per prime. That is more accurate than the closed
    form: at most 5.8 ulp from a 40-digit reference over `n <= 3000`, against
    31 ulp for the closed form at `k = 3.7`. It also differs from the reference in
    one place: `DivisorSigma[0., n]` is the divisor count (`4.` for `n = 6`),
    where the reference's closed form evaluates `0/0` and gives `Indeterminate`.
    `n = 1` is the exact `1` (the empty product), as in the reference.
  - **Anything else** (rational, symbolic, complex or arbitrary-precision `k`)
    — the multiplicative formula
    `sigma_k(n) = Product_i (p_i^((e_i+1) k) - 1) / (p_i^k - 1)` for
    `n = Product_i p_i^e_i`, built as an expression and evaluated, giving
    radical and symbolic forms.
- The sign of `n` is ignored; machine integers and GMP bigints are handled
  uniformly.
- In Gaussian mode the product runs over the first-quadrant associates
  (`Re > 0`, `Im >= 0`) of the Gaussian prime factors of `n`. This is the
  multiplicative definition — note it differs from naively summing
  `d^k` over `Divisors[n, GaussianIntegers -> True]`.
- Non-integer or zero `n` is left unevaluated.
- Messages (all through the `Quiet`/`Check` funnel, call left unevaluated):
  - `argr` for one argument, `argrx` for none;
  - `nonopt` for a non-rule after the second argument;
  - `optx` for an unknown option name;
  - `opttf` for a `GaussianIntegers` value other than `True` or `False`
    (including `Automatic`).

```mathematica
In[1]:= DivisorSigma[1, 20]
Out[1]= 42

In[2]:= DivisorSigma[2, 20]
Out[2]= 546

In[3]:= DivisorSigma[0, 12]
Out[3]= 6

In[4]:= DivisorSigma[-2, 10]
Out[4]= 13/10

In[5]:= DivisorSigma[1/2, 12]
Out[5]= (2 (-1 + 2 Sqrt[2]))/((-1 + Sqrt[2]) (-1 + Sqrt[3]))

In[6]:= DivisorSigma[k, {2, 3, 6}]
Out[6]= {(-1 + 2^(2 k))/(-1 + 2^k), (-1 + 3^(2 k))/(-1 + 3^k), ((-1 + 2^(2 k)) (-1 + 3^(2 k)))/((-1 + 2^k) (-1 + 3^k))}

In[7]:= DivisorSigma[2, {1, 2, 3, 4, 5}]
Out[7]= {1, 5, 10, 21, 26}

In[8]:= DivisorSigma[1, 3 + I]
Out[8]= 2 + 6 I

In[9]:= DivisorSigma[2, 6, GaussianIntegers -> True]
Out[9]= -30 + 20 I

In[10]:= DivisorSigma[2., 6]
Out[10]= 50.

In[11]:= DivisorSigma[1, 6, GaussianIntegers -> Automatic]
DivisorSigma::opttf: Value of option GaussianIntegers -> Automatic should be True or False.
Out[11]= DivisorSigma[1, 6, GaussianIntegers -> Automatic]
```

> **Packed arrays, NDArray and Compile.** `DivisorSigma[k, list]` over an
> `int64` buffer never boxes its elements.
>
> - **Element kernel.** Each element is factored by trial division in `int64`,
>   with no GMP allocation.
> - **Sieve.** A large buffer instead builds one table of smallest prime factors
>   up to `max |n|`, at one byte per integer. That is always done below 4·10⁶,
>   is done up to 5·10⁷ once the buffer has at least `4 sqrt(max)` elements,
>   and is never done above 5·10⁷. Every element then factors in O(Ω(n)), so
>   `DivisorSigma[1, Range[10^7]]` scales linearly: 0.56 s, where it used to
>   take 18 s and the reference takes 30 s.
> - **No int64 answer.** If any element has none (`k < 0`, overflow, `n = 0`),
>   the buffer goes back to the List path, which is itself fast now.
> - **Compile.** `Compile` lowers the scalar form to `OP_DSIG_I` and the array
>   form to the kernel. Both bail to the interpreter where there is no machine
>   answer.
>
> All of these surfaces, plus the scalar builtin, call one helper
> (`ndint_divisor_sigma_i64`), so they agree bit for bit.

## MoebiusMu

- `MoebiusMu[n]`: the Möbius function `mu(n)`. For
  `n = u p_1^k_1 ... p_m^k_m` with `u` a unit and `p_i` primes, `mu(n)` is `0`
  unless every `k_i` equals `1`, in which case it is `(-1)^m`. Equivalently
  `mu(1) = 1`, `mu(n) = (-1)^m` for a product of `m` distinct primes, and
  `mu(n) = 0` when `n` has a squared prime factor.

**Features**:
- `Listable`, `Protected`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `0`, `1`, or `-1`.
- The sign of `n` is ignored (`mu(-n) = mu(n)`).
- A non-real Gaussian-integer argument `Complex[a, b]` is handled over `Z[i]`:
  the input is factored into Gaussian primes (the unit factor does not count),
  giving `0` for a repeated Gaussian prime factor and `(-1)^m` otherwise.
- `MoebiusMu[0]` is `0` (every square divides 0), as in the reference. It was
  left unevaluated before v0.471. A non-integer `n` is left unevaluated; a
  wrong argument count issues a `MoebiusMu::argx` message.
- Machine surfaces share one int64 helper (`ndint_moebius_mu_i64`): the scalar
  path for `|n| <= 2^20`, the packed / `NDArray` kernel, and `Compile`
  (`OP_MU_I`).
  - A large buffer factors off one smallest-prime-factor sieve
    (`NDUnaryKernel.whole_i`), so `MoebiusMu[Range[10^7]]` takes 0.44 s and
    scales linearly.
  - The `Do`-loop of 10^5 calls takes 38 ms (the reference 36 ms).

```mathematica
In[1]:= MoebiusMu[11]
Out[1]= -1

In[2]:= MoebiusMu[10]
Out[2]= 1

In[3]:= MoebiusMu[1440]
Out[3]= 0

In[4]:= MoebiusMu[{4, 10, 17, 20}]
Out[4]= {0, 1, -1, 0}

In[5]:= MoebiusMu[10^50 + 1]
Out[5]= -1

In[6]:= MoebiusMu[5 + 6 I]
Out[6]= -1
```

> **Packed arrays.** Runs on an `int64` buffer. `MoebiusMu[0]` is undefined,
> so an array containing `0` takes the ordinary path and leaves that element
> unevaluated exactly as the scalar does.

## MangoldtLambda

- `MangoldtLambda[n]`: the von Mangoldt function -- `Log[p]` when `n = p^k` is a
  power of a prime `p` (`k >= 1`), and `0` otherwise.

**Features**:
- `Listable`, `Protected`.
- The result is exact: `MangoldtLambda[8] = Log[2]`, `N[MangoldtLambda[9]] =
  1.09861`. `MangoldtLambda[1] = 0`, and every `n <= 0` gives `0`, as in the
  reference.
- The prime-power test is the shared `mth_prime_power_base_{u64,mpz}` (the one
  `PrimePowerQ` uses): no full factorisation, so `MangoldtLambda[(10^20+39)^3]`
  is immediate.
- An inexact argument emits `MangoldtLambda::exact` and stays unevaluated; the
  message is gated on the arithmetic-probe mute, so `NSum`'s continuity probe at
  non-integer points stays silent (and invisible to `Check`). Rationals,
  symbols and Gaussian integers stay unevaluated silently.
- No packed / `NDArray` / `Compile` surface, deliberately: the answer is an exact
  `Log[p]`, which no machine dtype holds, and an inexact argument is an error,
  not a value.
- Its Dirichlet series closes: `Sum[MangoldtLambda[k]/k^s, {k, 1, Infinity}]` is
  `-Zeta'[s]/Zeta[s]` (for `Re s > 1` at a numeric `s`; `s = 1` diverges and
  stays), and `NSum` takes the same closed form at the working precision, as it
  does for the other arithmetic functions: the terms `log p/p^s` live on prime
  powers only and the tail falls like `N^(1-s)`, so no direct or extrapolated
  sum got past ~1 digit at `s = 1.68`.

```mathematica
In[1]:= Table[MangoldtLambda[n], {n, 1, 10}]
Out[1]= {0, Log[2], Log[3], Log[2], Log[5], 0, Log[7], Log[2], Log[3], 0}
```

## LiouvilleLambda

- `LiouvilleLambda[n]`: the Liouville function `lambda(n) = (-1)^Omega(n)`, where
  `Omega(n)` is the number of prime factors of `n` counted **with multiplicity**.
  For `n = u p_1^k_1 ... p_m^k_m` with `u` a unit and `p_i` primes,
  `lambda(n) = (-1)^(k_1 + ... + k_m)`. Equivalently `lambda(1) = 1`, and
  `lambda` is `1` when `n` is a product of an even number of primes and `-1`
  otherwise.

**Features**:
- `Listable`, `Protected`.
- Completely multiplicative: `lambda(m n) = lambda(m) lambda(n)`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly); the result is always `1` or `-1`.
- The sign of `n` is ignored (`lambda(-n) = lambda(n)`).
- Gaussian integers: `LiouvilleLambda[n, GaussianIntegers -> True]`, or a
  non-real Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]`
  and counts the Gaussian prime factors with multiplicity. Because `2` factors
  as `-i (1 + i)^2` in `Z[i]`, e.g. `LiouvilleLambda[2, GaussianIntegers -> True]`
  is `1` while `LiouvilleLambda[2]` is `-1`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `LiouvilleLambda::argt` message.

```mathematica
In[1]:= LiouvilleLambda[8]
Out[1]= -1

In[2]:= LiouvilleLambda[9]
Out[2]= 1

In[3]:= LiouvilleLambda[{1, 2, 3, 4, 5, 6}]
Out[3]= {1, -1, -1, 1, -1, 1}

In[4]:= LiouvilleLambda[10^30 + 1]
Out[4]= -1

In[5]:= LiouvilleLambda[2 + I]
Out[5]= -1

In[6]:= LiouvilleLambda[8, GaussianIntegers -> True]
Out[6]= 1
```

## PrimeOmega

> v0.472: `Options[PrimeOmega]` is `{GaussianIntegers -> False}`, as in
> the reference. An inexact argument draws `PrimeOmega::exact`, and option
> errors draw `::nonopt` / `::optx` / `::opttf` through the shared
> `nt_parse_gaussian_options`. Non-integers such as `PrimeOmega[1/12]` stay
> unevaluated, as in the reference.

- `PrimeOmega[n]`: the number of prime factors of `n` counted **with
  multiplicity**, `Omega(n)`. For `n = u p_1^k_1 ... p_m^k_m` with `u` a unit
  and `p_i` primes, `PrimeOmega[n] = k_1 + ... + k_m`. This is the exponent sum
  that `LiouvilleLambda` reduces via `lambda(n) = (-1)^Omega(n)`.

**Features**:
- `Listable`, `Protected`.
- Completely additive: `Omega(m n) = Omega(m) + Omega(n)`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly).
- `PrimeOmega[1]` (and `PrimeOmega[-1]`) is `0`; the sign of `n` is ignored
  (`Omega(-n) = Omega(n)`).
- Gaussian integers: `PrimeOmega[n, GaussianIntegers -> True]`, or a non-real
  Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]` and counts
  the Gaussian prime factors with multiplicity. Because `2` factors as
  `-i (1 + i)^2` in `Z[i]`, `PrimeOmega[12, GaussianIntegers -> True]` is `5`
  (from `(1 + i)^4 3`) while `PrimeOmega[12]` is `3`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `PrimeOmega::argt` message.

```mathematica
In[1]:= PrimeOmega[30]
Out[1]= 3

In[2]:= PrimeOmega[12]
Out[2]= 3

In[3]:= PrimeOmega[{4, 12, 24}]
Out[3]= {2, 3, 4}

In[4]:= PrimeOmega[30!]
Out[4]= 59

In[5]:= PrimeOmega[5 + 9 I]
Out[5]= 2

In[6]:= PrimeOmega[12, GaussianIntegers -> True]
Out[6]= 5
```

## PrimeNu

> v0.472: `Options[PrimeNu]` is `{GaussianIntegers -> False}`. An inexact
> argument draws `PrimeNu::exact`, and option errors draw `::nonopt` /
> `::optx` / `::opttf`.
>
> **Deliberate difference:** for `PrimeNu[5, Foo -> 1]` the reference returns
> `1` with `OptionValue::nodef`, which looks like a reference bug; Mathilda
> draws `PrimeNu::optx` and stays unevaluated.

- `PrimeNu[n]`: the number of **distinct** prime factors of `n`, `nu(n)`. For
  `n = u p_1^k_1 ... p_m^k_m` with `u` a unit and `p_i` distinct primes,
  `PrimeNu[n] = m` (independent of the exponents). `nu` and `Omega`
  (`PrimeOmega`) coincide exactly when `n` is square-free; for a prime power
  `PrimeNu[n]` is `1`.

**Features**:
- `Listable`, `Protected`.
- Additive on coprime arguments: `nu(m n) = nu(m) + nu(n)` when
  `GCD[m, n] == 1`.
- Computed directly from the prime factorisation (machine integers and GMP
  bigints handled uniformly).
- `PrimeNu[1]` (and `PrimeNu[-1]`) is `0`; the sign of `n` is ignored
  (`nu(-n) = nu(n)`).
- Gaussian integers: `PrimeNu[n, GaussianIntegers -> True]`, or a non-real
  Gaussian-integer argument `Complex[a, b]`, factors `n` over `Z[i]` and counts
  the distinct Gaussian prime factors. Because a rational prime `p ≡ 1 (mod 4)`
  splits into two conjugate Gaussian primes, e.g.
  `PrimeNu[105, GaussianIntegers -> True]` is `4` (from `3`, the split pair over
  `5`, and `7`) while `PrimeNu[105]` is `3`.
- Non-integer or zero `n` is left unevaluated; a wrong argument count issues a
  `PrimeNu::argt` message.
- Relations: for a square-free `n`, `MoebiusMu[n] == (-1)^PrimeNu[n]` and
  `LiouvilleLambda[n] == (-1)^PrimeNu[n]`.

```mathematica
In[1]:= PrimeNu[24]
Out[1]= 2

In[2]:= PrimeNu[105]
Out[2]= 3

In[3]:= PrimeNu[{4, 28, 180}]
Out[3]= {1, 2, 3}

In[4]:= PrimeNu[50!]
Out[4]= 15

In[5]:= PrimeNu[3 + I]
Out[5]= 2

In[6]:= PrimeNu[105, GaussianIntegers -> True]
Out[6]= 4
```

## EulerPhi

Gives the Euler totient function $\phi(n)$.
- `EulerPhi[n]`

**Features**:
- `Listable`, `Protected`.
- Counts the number of positive integers less than or equal to $n$ that are relatively prime to $n$.
- Returns 0 for $n = 0$, and handles negative integers via $\phi(-n) = \phi(n)$.
- Messages, as in the reference: a non-integer number (`EulerPhi[1/2]`,
  `EulerPhi[2.]`, `EulerPhi[3 + I]`) draws `EulerPhi::int`, and a wrong
  argument count draws `EulerPhi::argx`. Symbolic arguments stay unevaluated
  silently.
- Machine surfaces share one int64 helper (`ndint_euler_phi_i64`): the scalar
  path for `|n| <= 2^20`, the packed / `NDArray` kernel, and `Compile`
  (`OP_PHI_I`, at scalar and array shapes). A large buffer factors off one
  smallest-prime-factor sieve, so `EulerPhi[Range[10^7]]` takes 0.57 s.
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

> **Packed arrays.** Runs on an `int64` buffer, factoring by trial division
> in `int64`. `EulerPhi[3.]` is not `EulerPhi[3]`, so a real buffer takes the
> ordinary path.

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
- `Protected`. Not `Listable`: `SquareFreeQ[{1, 2, 3}]` stays unevaluated, as
  in the reference.
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
- `Modulus -> p` (p prime) tests a univariate polynomial over GF(p): the
  polynomial is reduced mod p, and it is square-free iff
  `PolynomialGCD[f, D[f, x], Modulus -> p]` is a constant
  (`SquareFreeQ[x^5 + 1, Modulus -> 5]` is `False`, since x^5 + 1 = (x + 1)^5
  there). `Modulus -> 0` is the ordinary test. A multivariate polynomial or a
  composite modulus stays unevaluated, without a message.
- `Extension -> ...` stays unevaluated, without a message.

Diagnostics:
- `SquareFreeQ[]` emits `SquareFreeQ::argb` and stays unevaluated.
- A non-`Rule` past the optional vars slot (e.g. `SquareFreeQ[1, 2, 3]`)
  emits `SquareFreeQ::nonopt` and stays unevaluated.
- `SquareFreeQ::par` ("Inappropriate parameter: ..."), as in the reference, for
  an unknown option, a bad option value, or `Modulus -> p` on a number
  (`SquareFreeQ[12, Modulus -> 5]`). The call stays unevaluated.

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

## IntegerPartitions

- `IntegerPartitions[n]`: all partitions of `n` in reverse-lexicographic order, parts descending. `Length[IntegerPartitions[n]] == PartitionsP[n]`.
- `IntegerPartitions[n, k]`: partitions into at most `k` parts. `All` or `Infinity`: any number.
- `IntegerPartitions[n, {k}]`: exactly `k` parts.
- `IntegerPartitions[n, {kmin, kmax}]`: between `kmin` and `kmax` parts.
- `IntegerPartitions[n, {kmin, kmax, dk}]`: `kmin, kmin + dk, …` parts. A negative `dk` counts down from `kmin` to `kmax`: `{6, 2, -2}` is lengths 6, 4 and 2, and `{2, 4, -1}` is empty. The rows stay in one reverse-lexicographic order; they are not grouped by length.
- `IntegerPartitions[n, kspec, sspec]`: parts drawn only from `sspec`.
  - A list: each entry is its own part, in the reversed order of the list, so `{3, 1, 2}` and `{1, 1, 2}` give the corresponding orders and duplicates.
  - A number `s`: `Range[s]`, so `IntegerPartitions[6, All, 5/2]` uses parts 1 and 2.
  - `All`: `Range[n]`.
- `IntegerPartitions[n, kspec, sspec, m]`:
  - `m > 0`: the first `m` rows; `m < 0`: the last `|m|`.
  - `None` or `0`: `{}`.
  - A sequence `{i}`, `{i, j}` or `{i, j, step}`, as in `Take`: `{1, -1, 3}` is every third row, and `{-1, 1, -1}` is all rows reversed.
  - A sequence that runs past either end emits `IntegerPartitions::take` and answers what is there: `{2, 30}` of 11 rows is rows 2–11, and `{20}` is `{}`. An index `0` gives `{}` silently.

`n` and the parts may be rational or negative. `n` may be a bignum when a length bound or a positive `m` keeps the answer finite: `IntegerPartitions[10^20, {3}, All, 2]` is `{{99999999999999999998, 1, 1}, {99999999999999999997, 2, 1}}`. Attributes: `Protected`.

**Messages.** All go through the `Quiet`/`Check` funnel, and the first bad position is reported:
- `rat`: position 1 not rational, symbols included.
- `nninfseq`: a bad length spec (negative, non-integer, `{Infinity}`, or a zero step).
- `rseq`: a part spec that is not `All`, a rational, or a list of rationals.
- `seq`: a bad fourth argument.
- `undef`: an infinitely large answer, from a zero part, or parts of both signs with no length bound.
- `take`: a sequence past the end.
- `argb`: fewer than 1 or more than 4 arguments.

**Implementation** (`src/partitions.c`). Four engines, picked by the shape of the part spec; all produce the same order:
- **Range mode** (`sspec` is `All` or a number): a word-sized depth-first walk with no part list. At depth `d` with remainder `r`, the candidate parts are the integers `max(1, ⌈r/(hi−d)⌉) … min(previous, r, smax, r−(lo−d−1))`, largest first. The length bounds prune exactly, and the stack grows only with the depth reached. A GMP twin handles bignum `n`.
- **Ascending part list** of positive machine integers: the same walk over the list's entries. It uses binary search for the first candidate, and a reachability table (`r` a sum of the remaining parts) when the smallest part exceeds 1.
- **Any other list of positive machine integers**: an int64 count-vector walk (the count of each priority part, largest first), pruned by reachability and by the suffix extremes of the parts.
- **Rational or negative parts**: the exact `mpq_t` count-vector enumerator.

Each row is one `List` node over shared Integer nodes. The finished result is stamped as an evaluated fixed point (`eval_stamp_list_literal`), so the evaluator does not re-walk millions of literal nodes.

| workload | before | after | reference |
|---|---|---|---|
| `Length[IntegerPartitions[50]]` | 2.0 s | 0.072 s | 0.053 s |
| `Length[IntegerPartitions[60]]` | 10.7 s | 0.36 s | 0.27 s |
| `IntegerPartitions[100, {4}]` | 1.33 s | 0.0011 s | 0.0011 s |
| `IntegerPartitions[60, All, {11, 7, 5, 3, 2, 1}]` | 0.128 s | 0.0031 s | 0.0022 s |
| `IntegerPartitions[300, {5}, Prime[Range[30]]]` | 0.99 s | 0.00024 s | 0.0015 s |

```
In[1]:= IntegerPartitions[5]
Out[1]= {{5}, {4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}}

In[2]:= IntegerPartitions[50, All, {6, 9, 20}]
Out[2]= {{20, 9, 9, 6, 6}, {20, 6, 6, 6, 6, 6}}

In[3]:= IntegerPartitions[8, {6, 2, -2}, All, 3]
Out[3]= {{7, 1}, {6, 2}, {5, 3}}

In[4]:= IntegerPartitions[6, All, All, {2, 4}]
Out[4]= {{5, 1}, {4, 2}, {4, 1, 1}}
```

## PartitionsP

Gives the number `p(n)` of unrestricted partitions of the integer `n` — the
number of ways to write `n` as a sum of positive integers, where order is
ignored. `PartitionsP[n] == Length[IntegerPartitions[n]]`.

- `PartitionsP[n]`: the partition count `p(n)` as an exact (arbitrary-precision)
  integer.

**Features**:
- `Protected`, `Listable` — `PartitionsP[{2, 4, 6}]` → `{2, 5, 11}`.
- **`n < 1000`**: a table of `p(0..999)` filled at startup by Euler's
  pentagonal-number-theorem recurrence
  `p(m) = Σ_k (-1)^(k-1) [p(m - g_k) + p(m - g_k')]`, with `g_k = k(3k∓1)/2`
  (about 40k GMP additions). A loop over small `n` is therefore a copy:
  `Do[PartitionsP[i], {i, 300}]` takes 0.05 ms, against 7 ms before.
- **`n ≥ 1000`**: the **Hardy–Ramanujan–Rademacher** exact formula evaluated in
  MPFR (Johansson, arXiv:1205.5991):
  `p(n) = Σ_{k=1}^N √(3/k)·(4/(24n−1))·A_k(n)·U(C/k)`, with
  `U(x) = cosh x − sinh x / x` and `C = (π/6)√(24n−1)`.
  - The `A_k(n)` use Selberg's cosine-sum identity.
  - `N` comes from Rademacher's rigorous remainder bound, so the rounded sum is
    provably exact. Without MPFR the recurrence answers.
- **Machine surfaces**: the int64 prefix of the table (`p(n)` fits for
  `n ≤ 405`) is the one machine answer, shared by:
  - the scalar builtin;
  - the packed/NDArray kernel (`NDKU_PartitionsP`, on `AWARE` and `INT64_OK`),
    which takes 12 ms for 10⁶ elements;
  - `Compile[]` (`OP_PARTP_I`) at scalar and rank-1 shapes.

  Past `n = 405` the packed path declines and the List path answers in
  bignums, and the Compile VM bails to the interpreter.
- `p(0) = 1`; `p(n) = 0` for `n < 0`. A non-integer, symbolic or bignum argument
  stays unevaluated; the wrong number of arguments emits `PartitionsP::argx`.

```
In[1]:= Table[PartitionsP[k], {k, 0, 12}]
Out[1]= {1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77}

In[2]:= PartitionsP[100]
Out[2]= 190569292

In[3]:= PartitionsP[4096]
Out[3]= 6927233917602120527467409170319882882996950147283323368445315320451
```

## PartitionsQ

Gives the number `q(n)` of partitions of the integer `n` into **distinct**
parts — equivalently (Euler), into **odd** parts (OEIS A000009).

- `PartitionsQ[n]`: the distinct-parts partition count `q(n)` as an exact
  (arbitrary-precision) integer.

**Features**:
- `Protected`, `Listable` — `PartitionsQ[{2, 4, 6}]` → `{1, 2, 4}`.
- **`n < 1000`**: a table filled at startup by the recurrence from
  `∏(1−x^k)·∏(1+x^k) = ∏(1−x^{2k})`:
  `q(m) = r(m) + Σ_k (−1)^(k−1) [q(m − g_k) + q(m − g_k')]`. Here
  `r(m) = (−1)^k` when `m/2` is the generalized pentagonal number `g_k`, and
  `0` otherwise.
- **`n ≥ 1000`**: the **Hardy–Ramanujan–Rademacher / Hagis** exact series in
  MPFR (Hagis, *Amer. J. Math.* 85 (1963) 213–222):
  `q(n) = (π/√(24n+1)) · Σ_{k odd} (1/k)·A_k(n)·I₁(π√(48n+2)/(12k))`, with
  `A_k(n) = Σ_{gcd(h,k)=1} cos(π(s(h,k) − s(2h,k)) − 2πnh/k)`.
  - **Exact angles.** `6k·s(h,k)` is an integer, so each Dedekind sum is run by
    reciprocity in double precision and rounded. An exact GMP fallback takes
    over if the rounding is not clean. Every angle is then the exact rational
    `πR/(6k)`.
  - **Per-term precision.** `|term_k| ≤ pref·e^(B/k)`, so each term is computed
    at its own magnitude plus guard bits. Most terms are summed in double, and
    only the first few dozen `k` need MPFR.
  - **Speed.** `q(10^4)` takes 0.0027 s and `q(10^5)` 0.024 s, against 0.125 s
    and 1.72 s before (0.034 s and 0.73 s for the reference). Over the same
    range the series now beats the recurrence from `n ≈ 700`.
  - **Termination.** As before, it is convergence-driven: the sum must lie
    within `1/4` of an integer and the last term must be below `1/8`.
    Correctness is pinned by the exhaustive series-equals-recurrence check up to
    `n = 1500` in `tests/test_partitionsq.c`; that suite went from 83 s to 4 s.
- **Machine surfaces**: as for `PartitionsP`. `q(n)` fits int64 for `n ≤ 769`
  (`OP_PARTQ_I`, `NDKU_PartitionsQ`).
- `q(0) = 1`; `q(n) = 0` for `n < 0`; a non-integer, symbolic or bignum
  argument stays unevaluated; `PartitionsQ::argx` for the wrong arity.

```
In[1]:= Table[PartitionsQ[k], {k, 0, 20}]
Out[1]= {1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 32, 38, 46, 54, 64}

In[2]:= PartitionsQ[100]
Out[2]= 444793
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

