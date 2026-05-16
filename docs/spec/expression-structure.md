# 1. Expression Structure

Everything in Mathilda is an `Expr`, which can be one of the following types:
- **Integer**: 64-bit signed integers (e.g., `123`, `-45`).
- **Real**: Double-precision floating-point numbers (e.g., `3.14`, `1.2E-3`).
- **Rational**: Exact fractions of integers, represented as `Rational[n, d]`.
- **Complex**: Complex numbers, represented as `Complex[re, im]`. The symbol `I` represents `Complex[0, 1]`.
- **Symbol**: Named identifiers (e.g., `x`, `Plus`, `Pi`, `Infinity`, `E`).
- **String**: Character sequences enclosed in double quotes (e.g., `"hello"`).
- **Function**: Compound expressions consisting of a `head` and zero or more `arguments` (e.g., `f[x, y]`, `{1, 2, 3}`).

