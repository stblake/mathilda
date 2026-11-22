# Bitwise operations

Bit-level operations on integers. Implemented in `src/bitwise/` (one builtin per
file, registered by `bitwise_init()`).

## BitLength

- `BitLength[n]`: The number of binary bits needed to represent the integer `n`.

**Features**:
- `Protected`, `Listable`. Threads element-wise over a list of integers, e.g.
  `BitLength[{0, 1, 2, 7, 8, 255, 256}]`.
- For `n > 0`, `BitLength[n]` is an efficient exact version of
  `Floor[Log[2, n]] + 1` -- it never converts through floating-point and is
  exact for arbitrarily large `n` (via GMP's `mpz_sizeinbase` in base 2).
- `BitLength[0]` is `0`.
- For `n < 0`, `BitLength[n]` is equivalent to `BitLength[BitNot[n]]`. In two's
  complement `BitNot[n] = -n - 1`, so `BitLength[-1]` is `0`, `BitLength[-2]` is
  `1`, and `BitLength[-2^k]` is `k`.
- Works for both machine integers and bignums. The packed / `NDArray` fast path
  (`src/ndinteger.c`) and the `Compile[]` lowering (`OP_BLEN_I`) both handle the
  full `int64` range, including `INT64_MIN` (`BitLength[-2^63]` is `63`), which
  the complement-based kernel computes without overflow.

```mathematica
In[1]:= BitLength[255]
Out[1]= 8

In[2]:= BitLength[256]
Out[2]= 9

In[3]:= BitLength[2^100]
Out[3]= 101

In[4]:= BitLength[100!]
Out[4]= 525

In[5]:= BitLength[{0, 1, 2, 7, 8, 255, 256}]
Out[5]= {0, 1, 2, 3, 4, 8, 9}

In[6]:= BitLength[{-1, -2, -8, -256}]
Out[6]= {0, 1, 3, 8}
```

**Diagnostics**. A non-integer argument emits `BitLength::int` and echoes the
call back unevaluated; a wrong-arity call emits `BitLength::argx`. A symbolic
argument is left unevaluated silently.

```mathematica
In[7]:= BitLength[3.5]
        BitLength::int: Integer expected at position 1 in BitLength[3.5].
Out[7]= BitLength[3.5]

In[8]:= BitLength[1, 2]
        BitLength::argx: BitLength called with 2 arguments; 1 argument is expected.
Out[8]= BitLength[1, 2]
```

## BitXor

- `BitXor[n1, n2, ...]`: The bitwise XOR of the integers `ni` -- the integer
  whose binary representation has ones at the positions where an odd number of
  the `ni` have ones.

**Features**:
- `Flat`, `Listable`, `OneIdentity`, `Orderless`, `Protected`. Takes any number
  of arguments; `BitXor[]` is `0` and `BitXor[n]` is `n`. Threads over lists
  (`BitXor[Range[120], 5]`, `BitXor[{1, 2}, 3, {4, 5}]`), including nested and
  empty ones (`BitXor[{}]` is `{}`).
- Negative integers are taken in two's complement (`BitXor[-2, -3]` is `3`).
  Integers of any size are exact: bignums go through GMP's `mpz_xor`, read in
  place and XORed straight into a presized result (no operand or result
  copies), and a result that fits in a machine integer is returned as one
  (`BitXor[2^70, 2^70 + 5]` is `5`).
- Automatic symbolic simplification, in this order:
  1. identical arguments cancel in pairs (`BitXor[x, x]` is `0`,
     `BitXor[x, y, x, y, x]` is `x`, `BitXor[x, 2.5, 2.5]` is `x`);
  2. an odd number of `-1` arguments is pulled out as `BitNot`
     (`BitXor[-1, x]` is `BitNot[x]`, `BitXor[1, -1, x]` is
     `BitNot[BitXor[1, x]]`);
  3. the remaining integers combine into one constant, dropped when it is `0`,
     pulled out as `BitNot` when it is `-1` (`BitXor[1, 2, x]` is
     `BitXor[3, x]`, `BitXor[-2, 1, x]` is `BitNot[x]`); any other negative
     constant stays inline (`BitXor[-5, x]`).

  `BitNot` appears only as an inert head.
- **Packed / NDArray**: an exact int64 kernel (`ndk_BitXor_ii`, plus
  vectorisable whole-buffer arms for array ⊕ scalar and array ⊕ array) in
  `src/ndinteger.c`, on the
  `AWARE` and `INT64_OK` lists. It handles array ⊕ scalar and array ⊕ array, and
  the n-ary form (`BitXor[v, w, 7]`, `BitXor[v, 1, 2]`) stays on the buffer, one
  pass per operand. A visible `NDArray` gives a visible `NDArray`. Operands that
  no int64 buffer can hold (a bignum, a symbol) take the exact List path.
- **Compile[]**: `OP_BXOR_I` (total, never bails), folded left for any arity.
  It lowers at scalar shape, array ⊕ scalar, array ⊕ array and the n-ary array
  form, and `CompileDiagnostics` reports `Compiled -> True` for each. Scalar,
  List, packed, `NDArray` and compiled results agree bit for bit, `INT64_MIN`
  included.

```
In[1]:= BitXor[61, 15]
Out[1]= 50

In[2]:= BitXor[3^100, 5^100]
Out[2]= 7888609052210118054116793203790437104101594041461795136781437785275168

In[3]:= BitXor[3333, 5555, 7777, 9999]
Out[3]= 8664

In[4]:= {BitXor[x, x, x], BitXor[1, 2, x], BitXor[-1, x], BitXor[y, x]}
Out[4]= {x, BitXor[3, x], BitNot[x], BitXor[x, y]}

In[5]:= Outer[BitXor, {1, 0}, {1, 0}]
Out[5]= {{0, 1}, {1, 0}}
```

**Diagnostics**. A concrete non-integer number that survives pair cancellation
(a Real, Rational or Complex) emits `BitXor::int`, naming its position, and the
call stays unevaluated. Symbolic arguments are silent. The message routes
through the message funnel, so `Quiet` and `Check` see it.

```
In[6]:= BitXor[2.5, 3]
        BitXor::int: Integer expected at position 1 in BitXor[2.5, 3].
Out[6]= BitXor[2.5, 3]
```

## BitAnd

- `BitAnd[n1, n2, ...]`: The bitwise AND of the integers `ni` -- the integer
  whose binary representation has ones at the positions where all of the `ni`
  have ones.

**Features**:
- `Flat`, `Listable`, `OneIdentity`, `Orderless`, `Protected`. Takes any number
  of arguments; `BitAnd[]` is `-1` (the identity) and `BitAnd[n]` is `n`.
  Threads over lists (`BitAnd[Range[70, 130], 17]`, `BitAnd[{1, 2}, 3, {4, 5}]`),
  including nested and empty ones (`BitAnd[{}]` and `BitAnd[{}, {}]` are `{}`).
- Negative integers are taken in two's complement (`BitAnd[-2, -3]` is `-4`).
  Integers of any size are exact: bignums go through GMP's `mpz_and`, read in
  place and ANDed straight into a presized result (no operand or result
  copies), and a result that fits in a machine integer is returned as one
  (`BitAnd[2^70 + 5, 7]` is `5`).
- Automatic symbolic simplification (AND is idempotent, `-1` is its identity
  and `0` absorbs):
  1. identical arguments collapse to one (`BitAnd[x, y, y, x]` is
     `BitAnd[x, y]`, `BitAnd[2.5, 2.5]` is `2.5`);
  2. the integers combine into one constant (`BitAnd[3, 5, x]` is
     `BitAnd[1, x]`). A constant of `0` absorbs everything, symbols and
     non-integers included (`BitAnd[1, 2, x]` and `BitAnd[0, 2.5]` are `0`); a
     constant of `-1` is dropped (`BitAnd[-1, x]` is `x`); any other constant
     stays inline (`BitAnd[-5, x]`).

  There are no complement or absorption rules: `BitAnd[x, BitNot[x]]` and
  `BitAnd[x, BitOr[x, y]]` stay as they are.
- **Packed / NDArray**: an exact int64 kernel (`ndk_BitAnd_ii`, plus
  vectorisable whole-buffer arms for array & scalar and array & array) in
  `src/ndinteger.c`, on the `AWARE` and `INT64_OK` lists. The n-ary form
  (`BitAnd[v, w, 7]`, `BitAnd[v, 7, 14]`) stays on the buffer, one pass per
  operand; a combined constant of `-1` costs no pass. A visible `NDArray` gives
  a visible `NDArray`. Operands that no int64 buffer can hold (a bignum, a
  symbol) take the exact List path. At 10^5 elements `BitAnd[v, 255]` and
  `BitAnd[v, w]` take about 0.04 and 0.05 ms, level with NumPy's `bitwise_and`.
- **Compile[]**: `OP_BAND_I` (total, never bails), folded left for any arity.
  It lowers at scalar shape, array & scalar, array & array and the n-ary array
  form, and `CompileDiagnostics` reports `Compiled -> True` for each. Scalar,
  List, packed, `NDArray` and compiled results agree bit for bit, `INT64_MIN`
  included.

```
In[1]:= BitAnd[61, 15]
Out[1]= 13

In[2]:= BitAnd[3^100, 5^100]
Out[2]= 503913279078601980753385379100266097026814443729

In[3]:= BitAnd[3333, 5555, 7777, 9999]
Out[3]= 1025

In[4]:= {BitAnd[x, y, y, x], BitAnd[3, 5, x], BitAnd[1, 2, x], BitAnd[-1, x]}
Out[4]= {BitAnd[x, y], BitAnd[1, x], 0, x}

In[5]:= Table[Sign[BitAnd[i, 2^2 + 2^3]], {i, 16}]
Out[5]= {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}
```

**Diagnostics**. None. A non-integer number (a Real, Rational or Complex) is
left in the call silently, exactly as a symbol is (`BitAnd[2.5, 3]` stays
unevaluated and `Check` takes the value branch).

## BitOr

- `BitOr[n1, n2, ...]`: The bitwise OR of the integers `ni` -- the integer
  whose binary representation has ones at the positions where any of the `ni`
  have ones.

**Features**:
- `Flat`, `Listable`, `OneIdentity`, `Orderless`, `Protected`. Takes any number
  of arguments; `BitOr[]` is `0` (the identity) and `BitOr[n]` is `n`. Threads
  over lists (`BitOr[Range[-40, 40], 17]`, `BitOr[{1, 2}, 8, {4, 5}]`),
  including nested and empty ones (`BitOr[{}]` and `BitOr[{}, {}]` are `{}`).
- Negative integers are taken in two's complement (`BitOr[-2, -3]` is `-1`).
  Integers of any size are exact: bignums go through GMP's `mpz_ior`, read in
  place and ORed straight into a presized result (no operand or result
  copies). A result that fits in a machine integer is returned as one, which
  happens whenever a negative operand sign-extends over a larger positive one
  (`BitOr[-5, 2^70]` is `-5`).
- Automatic symbolic simplification, the dual of BitAnd (OR is idempotent, `0`
  is its identity and `-1`, all bits on, absorbs):
  1. identical arguments collapse to one (`BitOr[x, y, y, x]` is
     `BitOr[x, y]`, `BitOr[2.5, 2.5]` is `2.5`);
  2. the integers combine into one constant (`BitOr[1, 2, x]` is
     `BitOr[3, x]`). A constant of `-1` absorbs everything, symbols and
     non-integers included (`BitOr[-1, 2, x]` and `BitOr[-1, 2.5]` are `-1`); a
     constant of `0` is dropped (`BitOr[0, x]` is `x`); any other constant
     stays inline (`BitOr[-5, x]`).

  There are no complement or absorption rules: `BitOr[x, BitNot[x]]` and
  `BitOr[x, BitAnd[x, y]]` stay as they are.
- **Packed / NDArray**: an exact int64 kernel (`ndk_BitOr_ii`, plus
  vectorisable whole-buffer arms for array | scalar and array | array) in
  `src/ndinteger.c`, on the `AWARE` and `INT64_OK` lists. The n-ary form
  (`BitOr[v, w, 7]`, `BitOr[v, 8, 16]`) stays on the buffer, one pass per
  operand; a combined constant of `0` costs no pass. A visible `NDArray` gives
  a visible `NDArray`. Operands that no int64 buffer can hold (a bignum, a
  symbol) take the exact List path. At 10^5 elements `BitOr[v, 255]`,
  `BitOr[v, w]` and `BitOr[u, v, w]` take about 0.035, 0.044 and 0.089 ms,
  level with NumPy's `bitwise_or` (0.034, 0.043, 0.085 ms).
- **Compile[]**: `OP_BOR_I` (total, never bails), folded left for any arity.
  It lowers at scalar shape, array | scalar, array | array and the n-ary array
  form, and `CompileDiagnostics` reports `Compiled -> True` for each, so
  auto-compilation picks it up. Scalar, List, packed, `NDArray` and compiled
  results agree bit for bit, `INT64_MIN` included.

```
In[1]:= BitOr[61, 15]
Out[1]= 63

In[2]:= BitOr[3^100, 5^100]
Out[2]= 7888609052210118054117297117069515706082347426840895402878464599718897

In[3]:= BitOr[3333, 5555, 7777, 9999]
Out[3]= 16383

In[4]:= {BitOr[x, y, y, x], BitOr[1, 2, x], BitOr[-1, 2, x], BitOr[0, x]}
Out[4]= {BitOr[x, y], BitOr[3, x], -1, x}

In[5]:= Outer[BitOr, {1, 0}, {1, 0}]
Out[5]= {{1, 1}, {1, 0}}
```

**Diagnostics**. None. A non-integer number (a Real, Rational or Complex) is
left in the call silently, exactly as a symbol is (`BitOr[2.5, 3]` stays
unevaluated and `Check` takes the value branch).

## BitNot

- `BitNot[n]`: The bitwise NOT of the integer `n` -- ones become zeros and
  vice versa.

**Features**:
- `Listable`, `Protected`. Integers are taken in two's complement with an
  unlimited number of digits, so `BitNot[n]` is exactly `-1 - n`: `BitNot[0]`
  is `-1` (all bits on) and `BitNot[7]` is `-8`. Threads over lists, including
  nested and empty ones (`BitNot[{}]` is `{}`).
- Exact at any size. Machine integers complement in int64, which is total
  (`BitNot[-2^63]` is `2^63 - 1`). Bignums go through GMP's `mpz_com` into a
  fresh node, never touching the operand; a result that fits in a machine
  integer is returned as one (`BitNot[2^63]` is `-2^63 - 1`, `BitNot[-2^63 - 1]`
  is `2^63`).
- Automatic symbolic simplification: BitNot is an involution, so
  `BitNot[BitNot[x]]` is `x`. This also completes BitXor's `-1` extraction:
  `BitXor[-1, BitNot[x]]` is `x`. There are no other rules: `BitNot[-1 - x]`
  and `BitAnd[x, BitNot[x]]` stay as they are.
- **Packed / NDArray**: an exact int64 unary kernel (`ndk_BitNot_ii`, with a
  vectorisable whole-buffer arm `ndk_BitNot_whole`) in `src/ndinteger.c`, on
  the `AWARE` and `INT64_OK` lists. A packed array stays packed; a visible
  `NDArray` gives a visible `NDArray`. A real buffer has no integer answer, so
  its elements stay unevaluated, as the scalar does. At 10^5 elements
  `BitNot[v]` takes about 0.038 ms, level with NumPy's `invert` (0.034 ms).
  At 10^7 elements it takes 22 ms against 18.5 ms; both are bound by memory
  bandwidth.
- **Compile[]**: `OP_BNOT_I` (total, never bails). It lowers at scalar and
  array shape, and `CompileDiagnostics` reports `Compiled -> True` for both.
  Auto-compilation picks it up: `Table[x + BitNot[Floor[x]], {x, 0., 99999.}]`
  takes 1.9 ms compiled and 84 ms interpreted, with the same answer. Scalar,
  List, packed, `NDArray` and compiled results agree bit for bit, `INT64_MIN`
  included.

```
In[1]:= BitNot[7]
Out[1]= -8

In[2]:= BitNot[Range[-10, 10]]
Out[2]= {9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11}

In[3]:= BitNot[2^70]
Out[3]= -1180591620717411303425

In[4]:= {BitNot[BitNot[x]], BitXor[-1, BitNot[x]]}
Out[4]= {x, x}

In[5]:= BitNot[]
BitNot::argx: BitNot called with 0 arguments; 1 argument is expected.
Out[5]= BitNot[]
```

**Diagnostics**. `BitNot::argx` for a call that does not have exactly one
argument (`BitNot[]`, `BitNot[1, 2]`). It goes through the message funnel, so
`Quiet` silences it and `Check` sees it. A non-integer number (a Real,
Rational or Complex) is left in the call silently, exactly as a symbol is
(`BitNot[2.5]` stays unevaluated and `Check` takes the value branch).

## BitSet

- `BitSet[n, k]`: The integer `n` with the bit corresponding to the
  coefficient of `2^k` set to 1 -- the bit `k` places from the right.

**Features**:
- `Listable`, `Protected`. Integers are taken in two's complement with an
  unlimited number of digits. The result grows to whatever size the bit needs
  (`BitSet[42, 9]` is `554`, `BitSet[0, 64]` is `2^64`). A negative `n`
  already has every high bit set (`BitSet[-5, 1]` is `-5`, `BitSet[-1, 100]` is
  `-1`). Threads over lists in either or both arguments (`BitSet[3, {0, 1, 5}]`
  is `{3, 3, 35}`). Lists of unequal length emit `Thread::tdlen`.
- Exact at any size. Machine integers set the bit in int64 when the answer
  fits: k ≤ 62, or any k for a negative `n`. Everything else (a bignum `n`, or
  bit 63 and up of a non-negative `n`) goes through GMP's `mpz_setbit` into a
  fresh node, never touching the operand. A result that fits a machine
  integer is returned as one (`BitSet[-2^64, 63]` is `-2^63`).
- Left unevaluated **silently**, with no message: a negative `k`
  (`BitSet[42, -1]`), a bignum `k`, any non-integer argument (`BitSet[2.5, 1]`,
  `BitSet[5, 1.]`, `BitSet[x, 1]`), and the wrong argument count (`BitSet[5]`,
  `BitSet[5, 1, 2]`). There are no symbolic rules:
  `BitSet[BitSet[x, 1], 1]` stays as it is.
- **Packed / NDArray**: an exact int64 binary kernel in `src/ndinteger.c`,
  through the machine helper `ndint_bitset_i64`, on the `AWARE` and `INT64_OK`
  lists. It has three whole-buffer arms:
  - a buffer of `n` with a scalar `k`: one OR with a constant mask;
  - a scalar `n` with a buffer of `k`;
  - two buffers.

  The last two are branch-free, with per-element shifts. On x86-64 they are
  also compiled for AVX2 (`vpsllvq`) and chosen once at start-up from the
  CPU's feature bits. Baseline SSE2 has no per-lane shift. A packed array stays
  packed and a visible `NDArray` gives a visible `NDArray`. An element with no
  int64 answer (bit 63 of a non-negative `n`, a negative `k`) abandons the
  buffer, and the List path answers exactly; nothing is ever truncated
  (`BitSet[Range[3], 63]` gives bignums). At 10^5 elements, against the NumPy
  equivalents:
  - `BitSet[v, 5]`: 0.049–0.054 ms vs `bitwise_or(v, 1 << 5)`
    0.046–0.051 ms, level.
  - `BitSet[7, kv]`: 0.061–0.069 ms vs `bitwise_or(7, left_shift(1, kv))`
    0.091–0.105 ms.
  - `BitSet[v, kv]`: 0.074–0.088 ms vs 0.117–0.129 ms.

  At 10^7 elements, `BitSet[v, kv]` takes 26 ms vs 45 ms. `BitSet[v, 5]` takes
  21 ms vs 20 ms; both are bound by memory bandwidth.
- **Compile[]**: `OP_BSET_I`, through the same `ndint_bitset_i64`. It bails
  to the interpreter where there is no int64 answer, so the bignum or the
  unevaluated form comes back exactly as uncompiled. It lowers at scalar
  shape and at array shape (array `n`, array `k`, or both), and
  `CompileDiagnostics` reports `Compiled -> True` for each. Auto-compilation
  follows. Scalar, List, packed, `NDArray` and compiled results agree bit for
  bit.

```
In[1]:= IntegerDigits[42, 2]
Out[1]= {1, 0, 1, 0, 1, 0}

In[2]:= BitSet[42, 0]
Out[2]= 43

In[3]:= BitSet[42, 4]
Out[3]= 58

In[4]:= BitSet[42, 9]
Out[4]= 554

In[5]:= IntegerDigits[%, 2]
Out[5]= {1, 0, 0, 0, 1, 0, 1, 0, 1, 0}

In[6]:= {BitSet[-5, 2], BitSet[0, 64], BitSet[42, -1]}
Out[6]= {-1, 18446744073709551616, BitSet[42, -1]}
```

## BitClear

- `BitClear[n, k]`: The integer `n` with the bit corresponding to the
  coefficient of `2^k` set to 0 -- the bit `k` places from the right.

**Features**:
- `Listable`, `Protected`. Integers are taken in two's complement with an
  unlimited number of digits. A non-negative `n` already has every high bit
  clear (`BitClear[5, 100]` is `5`). Clearing a high bit of a negative `n`
  subtracts `2^k`, so the result grows to whatever size the bit needs
  (`BitClear[-1, 0]` is `-2`, `BitClear[-1, 64]` is `-1 - 2^64`). Threads over
  lists in either or both arguments (`BitClear[7, {0, 1}]` is `{6, 5}`). Lists
  of unequal length emit `Thread::tdlen`.
- Exact at any size. Machine integers clear the bit in int64 when the answer
  fits: k ≤ 62, or any k for a non-negative `n`. Everything else (a bignum `n`,
  or bit 63 and up of a negative `n`) goes through GMP's `mpz_clrbit` into a
  fresh node, never touching the operand. A result that fits a machine integer
  is returned as one (`BitClear[2^63, 63]` is `0`).
- Left unevaluated **silently**, with no message, exactly as `BitSet` and
  `BitGet`: a negative `k` (`BitClear[5, -1]`), a bignum `k`, any non-integer
  argument (`BitClear[2.5, 1]`, `BitClear[x, 1]`), and the wrong argument count
  (`BitClear[5]`). There are no symbolic rules.
- **Packed / NDArray**: an exact int64 binary kernel in `src/ndinteger.c`,
  through the machine helper `ndint_bitclear_i64`, on the `AWARE` and
  `INT64_OK` lists. Its whole-buffer arms mirror `BitSet`'s: a buffer of `n`
  with a scalar `k` is one AND with a constant mask; a scalar `n` with a buffer
  of `k`, and two buffers, are branch-free per-element shifts, with the same
  run-time AVX2 clone on x86-64. An element with no int64 answer (bit 63 and up
  of a negative `n`, a negative `k`) abandons the buffer and the List path
  answers exactly; nothing is truncated.
- **Compile[]**: `OP_BCLR_I`, through the same `ndint_bitclear_i64`. It bails
  to the interpreter where there is no int64 answer, so the bignum or the
  unevaluated form comes back exactly as uncompiled. It lowers at scalar and at
  array shape, and `CompileDiagnostics` reports `Compiled -> True` for each.
  Scalar, List, packed, `NDArray` and compiled results agree bit for bit.

```
In[1]:= BitClear[255, 0]
Out[1]= 254

In[2]:= BitClear[255, {0, 7}]
Out[2]= {254, 127}

In[3]:= {BitClear[-1, 0], BitClear[5, 100], BitClear[2^70 + 5, 70]}
Out[3]= {-2, 5, 5}

In[4]:= BitClear[-1, 64]
Out[4]= -18446744073709551617

In[5]:= BitClear[5, -1]
Out[5]= BitClear[5, -1]
```

## BitGet

- `BitGet[n, k]`: The bit of the integer `n` corresponding to the coefficient
  of `2^k` -- the bit `k` places from the right -- as `0` or `1`.

**Features**:
- `Listable`, `Protected`. Integers are taken in two's complement with an
  unlimited number of digits. Every bit of a negative `n` above its magnitude
  is 1 (`BitGet[-1, 100]` is `1`, `BitGet[-42, 3]` is `0`), and bits above a
  non-negative `n` are 0 (`BitGet[42, 1000]` is `0`). Threads over lists in
  either or both arguments (`BitGet[5, {0, 1, 2, 3}]` is `{1, 0, 1, 0}`).
  Lists of unequal length emit `Thread::tdlen`.
- Exact at any size, and the answer is always the machine integer 0 or 1.
  Machine integers are read with a shift and a mask. The shift count is
  clamped to 63, which is exact because bit 63 is the sign bit and every
  higher bit equals it. A bignum `n` is read in place by GMP's `mpz_tstbit`,
  with no allocation and the operand untouched.
- Left unevaluated **silently**, with no message: a negative `k`
  (`BitGet[42, -1]`), a bignum `k`, any non-integer argument (`BitGet[2.5, 1]`,
  `BitGet[5, 1.]`, `BitGet[x, 1]`), and the wrong argument count (`BitGet[5]`,
  `BitGet[5, 1, 2]`). There are no symbolic rules:
  `BitGet[BitSet[x, 1], 1]` stays as it is.
- **Packed / NDArray**: an exact int64 binary kernel in `src/ndinteger.c`,
  through the machine helper `ndint_bitget_i64`, on the `AWARE` and `INT64_OK`
  lists. It has three whole-buffer arms:
  - a buffer of `n` with a scalar `k`: one constant shift and mask;
  - a scalar `n` with a buffer of `k`;
  - two buffers.

  The last two are branch-free, with per-element shifts. On x86-64 they get
  the same AVX2 clone (`vpsrlvq`) as BitSet, chosen once at start-up. Because
  any `k >= 0` has an int64 answer, the only decline is a negative `k`; those
  elements stay unevaluated and the others answer (`BitGet[Range[3], {0, -1,
  2}]` is `{1, BitGet[2, -1], 0}`). A packed array stays packed and a visible
  `NDArray` gives a visible `NDArray`. Against NumPy's `(v >> k) & 1` on
  random int64 data (x86-64, Mathilda mean time vs NumPy best of 3):

  | n    | `BitGet[v, 3]`   | `BitGet[7, kv]`  | `BitGet[v, kv]`  |
  |------|------------------|------------------|------------------|
  | 10^5 | 0.037 vs 0.056 ms | 0.044 vs 0.060 ms | 0.055 vs 0.067 ms |
  | 10^7 | 18.1 vs 23.6 ms  | 15.8 vs 23.8 ms  | 19.9 vs 27.5 ms  |
- **Compile[]**: `OP_BGET_I`, through the same `ndint_bitget_i64`. It bails
  to the interpreter only for a negative `k`, which then comes back
  unevaluated exactly as uncompiled. It lowers at scalar shape and at array
  shape (array `n`, array `k`, or both), and `CompileDiagnostics` reports
  `Compiled -> True` for each. Auto-compilation follows. Scalar, List, packed,
  `NDArray` and compiled results agree bit for bit.

```
In[1]:= IntegerDigits[42, 2]
Out[1]= {1, 0, 1, 0, 1, 0}

In[2]:= BitGet[42, 0]
Out[2]= 0

In[3]:= BitGet[42, 3]
Out[3]= 1

In[4]:= BitGet[Range[30, 60], 3]
Out[4]= {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}

In[5]:= {BitGet[-1, 100], BitGet[-2^100, 99], BitGet[42, -1]}
Out[5]= {1, 0, BitGet[42, -1]}
```

## BitShiftLeft

- `BitShiftLeft[n, k]`: shifts the binary bits in the integer `n` to the left
  by `k` places, padding with zeros on the right -- `n 2^k`.
- `BitShiftLeft[n]`: shifts one bit to the left.

## BitShiftRight

- `BitShiftRight[n, k]`: shifts the binary bits in the integer `n` to the right
  by `k` places, dropping bits shifted past the units position --
  `IntegerPart[n / 2^k]`.
- `BitShiftRight[n]`: shifts one bit to the right.

**Features** (both heads; `BitShiftRight[n, k]` is `BitShiftLeft[n, -k]` on
every argument, and one core in `src/bitwise/bitshiftleft.c` serves both):
- `Listable`, `Protected`. A negative `k` shifts the other way
  (`BitShiftLeft[5, -1]` is `2`, `BitShiftRight[5, -2]` is `20`). Threads over
  lists in either or both arguments; lists of unequal length emit
  `Thread::tdlen`.
- A right shift **truncates toward zero**, as `IntegerPart` does, so it is not
  the arithmetic shift: `BitShiftRight[-5, 1]` is `-2` (NumPy's `right_shift`
  gives `-3`), `BitShiftRight[-1, 100]` is `0`, `BitShiftRight[-2^63, 63]` is
  `-1`.
- Exact at any size. Machine integers use the shared helper
  `ndint_bitshiftleft_i64`; a left shift that leaves int64 goes to GMP's
  `mpz_mul_2exp`, and a bignum right shift to `mpz_tdiv_q_2exp`, normalised
  back to int64 when it fits (`BitShiftLeft[2^100, -100]` is the machine `1`).
  Zero shifts to `0` by any count.
- Any other number `n` is scaled by `2^k`: `BitShiftLeft[2.5, 1]` is `5.`,
  `BitShiftRight[1/2, 1]` is `1/4`. Machine reals use `ldexp` and MPFR reals
  `mpfr_mul_2si`, so the scaling is exact and keeps the precision. A complex
  `n` shifts each part, so integer parts still truncate:
  `BitShiftRight[3 + 5 I, 1]` is `1 + 2 I`, `BitShiftRight[Complex[1/2, 3], 1]`
  is `1/4 + I`. `Indeterminate` and the infinities are unchanged.
- Messages, through the Quiet/Check funnel:
  - `::intm` when `k` is a number but not a machine-sized integer
    (`BitShiftLeft[5, 1.5]`, `BitShiftLeft[5, 2.]`, `BitShiftLeft[0, 2^70]`);
  - `::argt` for 0 or 3+ arguments;
  - `General::ovfl` and `Overflow[]` for a left shift too large for GMP to
    represent (`BitShiftLeft[5, 2^63 - 1]`);
  - `General::unfl` and `Underflow[]` for a non-integer shifted right that far
    (`BitShiftRight[3/2, 2^63 - 1]`; a machine real gives `0.`).

  A symbolic argument (`BitShiftLeft[x, 2]`, `BitShiftLeft[Pi, 1]`) stays
  unevaluated with no message. Zero shifts to `0` by every count, including
  `BitShiftRight[0, -2^63]`, whose count has no int64 negation.
- **Packed / NDArray**: exact int64 kernels in `src/ndinteger.c` on the `AWARE`
  and `INT64_OK` lists -- binary for `[n, k]`, unary for the one-argument form.
  Whole-buffer arms:
  - a buffer of `n` with a scalar `k`: constant logical shifts;
  - a scalar `n` with a buffer of `k`, and two buffers: a fast single-direction
    pass over per-element counts. If any count needs the other direction, a
    general pass reruns the buffer.

  x86 has no 64-bit arithmetic right shift short of AVX-512. So the arms use
  logical shifts only:
  - a left shift fits iff `(n ^ sign) >> (63 - k)` is 0;
  - a right shift is `sign * (|n| >> k)`.

  On x86-64 the per-element-count passes use AVX2 intrinsics, picked at
  start-up. `vpsllvq`/`vpsrlvq` are defined for every count, so the lanes need
  no masking and the range test is one OR of the counts. A left shift with no
  int64 answer abandons the buffer, and the List path answers with bignums. It
  never wraps, unlike NumPy's `left_shift`. A packed array stays packed and a
  visible `NDArray` gives a visible `NDArray`. A real buffer is answered by
  the List path, element by element. The binary kernel sees a scalar `k` only
  as a double, so it cannot tell `2` from `2.`, and `2.` must give `::intm`.

  Against NumPy on random int64 data (x86-64 i9-9880H, best of 5, n = 10^5;
  all arms are memory-bound at parity by 10^7):

  | arm                      | Mathilda  | NumPy     |
  |--------------------------|-----------|-----------|
  | `BitShiftLeft[v, 3]`     | 0.033 ms  | 0.031 ms  |
  | `BitShiftRight[v, 3]`    | 0.035 ms  | 0.031 ms  |
  | `BitShiftLeft[v, kv]`    | 0.043 ms  | 0.041 ms  |
  | `BitShiftRight[v, kv]`   | 0.043 ms  | 0.045 ms  |
  | `BitShiftLeft[7, kv]`    | 0.032 ms  | 0.031 ms  |
  | `BitShiftRight[v]`       | 0.034 ms  | 0.031 ms  |

  The remaining 1-3 µs is the per-call floor every buffer head pays
  (`BitNot[v]` 0.031 ms vs `np.invert` 0.030 ms), not the kernels.
- **Compile[]**: `OP_BSHL_I` / `OP_BSHR_I`, through the same helpers; the
  one-argument form gets a constant `k = 1`. They bail to the interpreter only
  where a left shift leaves int64, and the interpreter then answers exactly.
  Both lower at scalar and array shape (array `n`, array `k`, or both), and
  `CompileDiagnostics` reports `Compiled -> True` for each. Auto-compilation
  follows. Scalar, List, packed, `NDArray` and compiled results agree bit for
  bit.

```
In[1]:= BitShiftLeft[32, 3]
Out[1]= 256

In[2]:= BitShiftLeft[2^200 - 1, 7]
Out[2]= 205688069665150755269371147819668813122841983204197482918576000

In[3]:= BitShiftLeft[{1, 2, 4, 16}, 2]
Out[3]= {4, 8, 16, 64}

In[4]:= BitShiftRight[32, 7]
Out[4]= 0

In[5]:= BitShiftRight[2^200 - 1, 7]
Out[5]= 12554203470773361527671578846415332832204710888928069025791

In[6]:= BitShiftRight[{1, 2, 4, 16}, 2]
Out[6]= {0, 0, 1, 4}

In[7]:= {BitShiftRight[-5, 1], BitShiftLeft[5, -1], BitShiftRight[3 + 5 I, 1]}
Out[7]= {-2, 2, 1 + 2*I}
```
