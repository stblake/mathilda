# Number Theory

This tutorial is a hands-on tour of Mathilda's number theory — the arithmetic of
the integers. You will compute greatest common divisors and least common
multiples, run the extended Euclidean algorithm to find Bézout coefficients,
work in modular arithmetic with powers and orders, test and count primes, factor
integers, evaluate classical arithmetic functions like Euler's totient, and
expand numbers as continued fractions. Everything here is exact: Mathilda carries
arbitrary-precision integers, so nothing rounds or overflows. The tutorial ends
with a complete, working RSA encryption and decryption that uses almost every
tool introduced along the way.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself (without the prompt) and you will see the same
`Out[...]`.

## Divisibility: GCD and LCM

The greatest common divisor and least common multiple are the two fundamental
measures of how integers share factors. `GCD` returns the largest integer
dividing all of its arguments; `LCM` returns the smallest positive integer they
all divide:

```mathematica
In[1]:= GCD[24, 36]
Out[1]= 12

In[2]:= LCM[4, 6]
Out[2]= 12

In[3]:= GCD[12, 18, 30]
Out[3]= 6
```

`In[1]` finds that `12` is the largest number dividing both `24` and `36`, and
`In[2]` that `12` is the smallest number both `4` and `6` divide. Both functions
take any number of arguments, so `In[3]` takes the gcd of three integers at once,
finding the `6` common to all of them.

The two are linked by a classic identity: for any pair of integers, the product
of their gcd and lcm equals the product of the numbers themselves. We can check
it directly:

```mathematica
In[1]:= GCD[24, 36] * LCM[24, 36]
Out[1]= 864

In[2]:= 24 * 36
Out[2]= 864
```

`In[1]` multiplies `GCD[24, 36] = 12` by `LCM[24, 36] = 72` to get `864`, which
`In[2]` confirms is exactly `24 × 36`. The shared factors counted once in the gcd
and the missing factors supplied by the lcm together reconstruct the full
product.

## The extended Euclidean algorithm

`ExtendedGCD` does more than report the gcd: it also returns the *Bézout
coefficients*, a pair of integers expressing the gcd as a linear combination of
the inputs. For `ExtendedGCD[a, b]` returning `{g, {s, t}}`, you always have
`s a + t b == g`:

```mathematica
In[1]:= ExtendedGCD[24, 36]
Out[1]= {12, {-1, 1}}

In[2]:= 12*2 + 36*(-1)
Out[2]= -12
```

`In[1]` reports that `GCD[24, 36] = 12` and that `(-1)·24 + 1·36 = 12`. `In[2]`
is a deliberate reminder to read the coefficients in the right order — using them
backwards gives `-12`, not the gcd. With the correct pairing,
`(-1)·24 + 1·36 = -24 + 36 = 12`, as promised.

A second example shows the coefficients can be large even for small inputs:

```mathematica
In[1]:= ExtendedGCD[240, 46]
Out[1]= {2, {-9, 47}}
```

Here `GCD[240, 46] = 2`, and the certificate is `(-9)·240 + 47·46 = -2160 + 2162
= 2`. These coefficients are precisely what is needed to invert numbers modulo
another — the workhorse behind modular division, which the RSA example at the end
relies on.

## Modular arithmetic

Working "modulo `n`" means caring only about remainders after division by `n`.
`Mod`, `Quotient`, and `QuotientRemainder` are the basic tools; the first gives
the remainder, the second the integer quotient, and the third both at once:

```mathematica
In[1]:= Mod[17, 5]
Out[1]= 2

In[2]:= Quotient[17, 5]
Out[2]= 3

In[3]:= QuotientRemainder[17, 5]
Out[3]= {3, 2}

In[4]:= Mod[-17, 5]
Out[4]= 3
```

`In[1]`–`In[3]` decompose `17 = 3·5 + 2`. `In[4]` shows Mathilda's `Mod` always
returns a result with the sign of the divisor, so `Mod[-17, 5]` is `3` (the
non-negative residue) rather than `-2` — the convention number theorists want,
where every residue class has a representative in `0, 1, …, n-1`.

Raising a number to a power modulo `n` comes up constantly, and computing the
full power first would be hopelessly large. `PowerMod[a, b, n]` reduces at every
step, so it stays fast and small:

```mathematica
In[1]:= PowerMod[2, 10, 1000]
Out[1]= 24

In[2]:= PowerMod[7, 256, 13]
Out[2]= 9

In[3]:= PowerMod[3, -1, 7]
Out[3]= 5
```

`In[1]` computes `2^10 = 1024`, reduced mod `1000` to `24`. `In[2]` evaluates a
much larger power, `7^256`, mod `13` without ever forming the giant integer. The
third form is the most useful: a *negative* exponent asks for the modular
inverse, so `In[3]` returns the number whose product with `3` is `1` mod `7` —
indeed `3·5 = 15 = 2·7 + 1`.

The *multiplicative order* of `a` mod `n` is the smallest positive power that
brings `a` back to `1`. `MultiplicativeOrder` computes it:

```mathematica
In[1]:= MultiplicativeOrder[2, 7]
Out[1]= 3

In[2]:= MultiplicativeOrder[3, 7]
Out[2]= 6
```

`In[1]` reports that `2^3 = 8 = 1` mod `7`, and no smaller power works. `In[2]`
shows that `3` has order `6` — the full size of the multiplicative group mod `7`
— which makes `3` a *primitive root*: its powers cycle through every nonzero
residue. `PrimitiveRoot` finds one such generator, and `PrimitiveRootList` lists
all of them:

```mathematica
In[1]:= PrimitiveRoot[7]
Out[1]= 3

In[2]:= PrimitiveRootList[7]
Out[2]= {3, 5}
```

`In[1]` returns `3`, the smallest primitive root mod `7`, matching the order-`6`
result above. `In[2]` shows there are exactly two primitive roots, `3` and `5`;
the count is `EulerPhi[EulerPhi[7]] = EulerPhi[6] = 2`, which the next section's
totient function explains.

## Primes

`PrimeQ` is the basic primality test, returning `True` or `False`:

```mathematica
In[1]:= PrimeQ[97]
Out[1]= True

In[2]:= PrimeQ[91]
Out[2]= False
```

`97` is prime (`In[1]`), but `91 = 7·13` is not (`In[2]`) — a number that looks
prime at a glance but isn't. To count primes, `PrimePi[x]` gives the number of
primes not exceeding `x`, and `NextPrime[x]` returns the first prime strictly
greater than `x`:

```mathematica
In[1]:= PrimePi[100]
Out[1]= 25

In[2]:= NextPrime[100]
Out[2]= 101

In[3]:= NextPrime[1000000]
Out[3]= 1000003
```

There are `25` primes up to `100` (`In[1]`), the first one past `100` is `101`
(`In[2]`), and even past a million `NextPrime` finds `1000003` instantly
(`In[3]`).

To break a composite number into its prime building blocks, use `FactorInteger`.
It returns a list of `{prime, exponent}` pairs:

```mathematica
In[1]:= FactorInteger[360]
Out[1]= {{2, 3}, {3, 2}, {5, 1}}

In[2]:= FactorInteger[1000000]
Out[2]= {{2, 6}, {5, 6}}
```

`In[1]` reads as `360 = 2^3 · 3^2 · 5`, and `In[2]` as `1000000 = 2^6 · 5^6` —
the prime factorisation of `10^6`. A related question is whether a number is
*square-free*, meaning no prime appears more than once in its factorisation.
`SquareFreeQ` answers it:

```mathematica
In[1]:= SquareFreeQ[30]
Out[1]= True

In[2]:= SquareFreeQ[12]
Out[2]= False
```

`30 = 2·3·5` has every prime to the first power, so it is square-free (`In[1]`);
`12 = 2^2·3` contains the square `4`, so it is not (`In[2]`).

## Arithmetic functions

Euler's totient `EulerPhi[n]` counts the integers in `1, …, n` that are coprime
to `n` — equivalently, the size of the multiplicative group mod `n`:

```mathematica
In[1]:= EulerPhi[10]
Out[1]= 4

In[2]:= EulerPhi[36]
Out[2]= 12
```

The four numbers coprime to `10` are `1, 3, 7, 9`, so `EulerPhi[10] = 4`
(`In[1]`); for `36 = 2^2·3^2` the totient is `36·(1-1/2)·(1-1/3) = 12` (`In[2]`).

The totient is the heart of *Euler's theorem*: if `a` is coprime to `n`, then
`a^EulerPhi[n] = 1` mod `n`. We can verify it directly with `PowerMod`:

```mathematica
In[1]:= PowerMod[3, EulerPhi[10], 10]
Out[1]= 1
```

Since `GCD[3, 10] = 1`, raising `3` to the power `EulerPhi[10] = 4` modulo `10`
must give `1` — and `In[1]` confirms it (`3^4 = 81 = 8·10 + 1`). This identity is
exactly why RSA decryption works, as we will see shortly.

Coprimality also governs binomial coefficients, which Mathilda computes exactly
with `Binomial`:

```mathematica
In[1]:= Binomial[10, 3]
Out[1]= 120
```

`Binomial[10, 3]` is the number of ways to choose `3` items from `10`, namely
`120`.

## Continued fractions

Every rational number has a finite *continued fraction* expansion — a list of
integers `{a0, a1, a2, …}` standing for `a0 + 1/(a1 + 1/(a2 + …))`.
`ContinuedFraction` produces the list, and `FromContinuedFraction` reconstructs
the number:

```mathematica
In[1]:= ContinuedFraction[415/93]
Out[1]= {4, 2, 6, 7}

In[2]:= FromContinuedFraction[{4, 2, 6, 7}]
Out[2]= 415/93
```

`In[1]` expands `415/93` as `4 + 1/(2 + 1/(6 + 1/7))`, and `In[2]` folds the same
list back into `415/93` exactly — the two operations are inverse. The expansion
mirrors the Euclidean algorithm: the partial quotients `4, 2, 6, 7` are precisely
the quotients you get dividing `415` by `93`, then `93` by the remainder, and so
on.

Irrational numbers have *infinite* expansions, so for them you ask for a fixed
number of terms with a second argument. The famous patterns appear immediately:

```mathematica
In[1]:= ContinuedFraction[Sqrt[2], 10]
Out[1]= {1, 2, 2, 2, 2, 2, 2, 2, 2, 2}

In[2]:= ContinuedFraction[Sqrt[7], 8]
Out[2]= {2, 1, 1, 1, 4, 1, 1, 1}

In[3]:= ContinuedFraction[Pi, 8]
Out[3]= {3, 7, 15, 1, 292, 1, 1, 1}
```

`Sqrt[2]` is the purely periodic `[1; 2, 2, 2, …]` (`In[1]`); `Sqrt[7]` has the
periodic block `1, 1, 1, 4` (`In[2]`) — every quadratic irrational is eventually
periodic. `Pi` (`In[3]`) shows no pattern, but its early terms are the key to
good rational approximations: truncating after the large `15` and folding back
gives the celebrated approximation `355/113`:

```mathematica
In[1]:= FromContinuedFraction[{3, 7, 15, 1}]
Out[1]= 355/113
```

This is the best rational approximation to `Pi` with a denominator under a
thousand, accurate to six decimal places — recovered here from just four
continued-fraction terms.

## Putting it together: textbook RSA

We now assemble the pieces into a complete RSA cryptosystem, the classic
public-key scheme whose security rests on the difficulty of factoring. We use the
historic textbook primes `p = 61` and `q = 53` so every step runs and is easy to
follow. The session below is continuous: each line builds on the variables
defined above it.

First, pick two primes and form the modulus `n = p q`:

```mathematica
In[1]:= p = 61
Out[1]= 61

In[2]:= q = 53
Out[2]= 53

In[3]:= n = p q
Out[3]= 3233
```

The modulus `n = 3233` is the public part everyone can see; its prime factors
`p` and `q` are the secret. Next compute the totient `EulerPhi[n]`, which for a
product of two distinct primes equals `(p-1)(q-1)`:

```mathematica
In[1]:= phi = EulerPhi[n]
Out[1]= 3120

In[2]:= phi == (p - 1)(q - 1)
Out[2]= True
```

`EulerPhi[3233] = 3120`, and `In[2]` verifies it matches `(61-1)(53-1) = 60·52 =
3120`. This `phi` is secret too — computing it requires knowing the factorisation.

Now choose a public exponent `e` coprime to `phi`. The conventional choice `17`
works:

```mathematica
In[1]:= e = 17
Out[1]= 17

In[2]:= GCD[e, phi]
Out[2]= 1

In[3]:= ExtendedGCD[e, phi]
Out[3]= {1, {-367, 2}}
```

`In[2]` confirms `GCD[17, 3120] = 1`, so `e` is a valid public exponent. The
extended gcd in `In[3]` already hints at the private key: it certifies
`(-367)·17 + 2·3120 = 1`, so `-367` is an inverse of `17` mod `3120`. The private
exponent `d` is exactly this inverse, taken as a positive residue, which
`PowerMod` gives directly:

```mathematica
In[1]:= d = PowerMod[e, -1, phi]
Out[1]= 2753

In[2]:= Mod[e d, phi]
Out[2]= 1
```

`In[1]` returns `d = 2753` (and indeed `2753 = -367 + 3120`, matching the
extended-gcd certificate). `In[2]` confirms the defining relation `e·d = 1` mod
`phi`. The public key is now `(n, e) = (3233, 17)` and the private key is
`(n, d) = (3233, 2753)`.

Finally, encrypt and decrypt a message. Take the plaintext number `msg = 42`;
encryption raises it to the public exponent mod `n`, decryption raises the
ciphertext to the private exponent:

```mathematica
In[1]:= msg = 42
Out[1]= 42

In[2]:= c = PowerMod[msg, e, n]
Out[2]= 2557

In[3]:= PowerMod[c, d, n]
Out[3]= 42

In[4]:= PowerMod[c, d, n] == msg
Out[4]= True
```

The message `42` encrypts to the ciphertext `2557` (`In[2]`), which looks nothing
like the original. Decrypting with the private exponent in `In[3]` returns `42`
exactly, and `In[4]` confirms the round trip. The reason it works is Euler's
theorem from earlier: because `e·d = 1` mod `phi`, raising the message to `e`
then to `d` is raising it to a power `≡ 1` mod `phi`, which by Euler's theorem
leaves it unchanged mod `n`. Anyone with the public key can encrypt, but only the
holder of `d` — derivable only by factoring `n` — can decrypt.

## Where to next

You can now compute gcds and lcms, run the extended Euclidean algorithm for
Bézout coefficients, do modular arithmetic with powers, inverses, and orders,
test, count, and factor primes, evaluate Euler's totient and verify Euler's
theorem, expand numbers as continued fractions, and assemble all of it into a
working public-key cryptosystem. This is the integer-arithmetic foundation that
sits beneath much of the rest of Mathilda.

- **[6. Algebra](06-algebra.md)** — move from the integers to polynomials and
  rational expressions: expand and factor, take polynomial gcds, simplify, and
  solve equations and systems.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions. The
  [number theory](../documentation/number-theory/index.md) section documents
  every function used above, including the complete options for `FactorInteger`,
  `PowerMod`, and `ContinuedFraction`.
