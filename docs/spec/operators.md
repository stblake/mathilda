# 3. Operators and Precedence

| Operator | FullForm | Precedence | Association |
|----------|----------|------------|-------------|
| `[[...]]`| `Part` | 100 | Left |
| `f[x]`   | `f[x]` | 1000 | Left |
| `_`, `__`, `___` | `Blank` | 730 | None |
| `?`      | `PatternTest` | 680 | None |
| `@`      | `Prefix` | 620 | Right |
| `//`     | `Postfix` | 70 | Left |
| `&`      | `Function`| 90 | Left |
| `@@`     | `Apply` | 620 | Right |
| `/@`     | `Map`   | 620 | Right |
| `<>`     | `StringJoin` | 600 | Left |
| `^`      | `Power` | 590 | Right |
| `*`      | `Times` | 400 | Left |
| `/`      | `Divide`| 470 | Left |
| `+`, `-` | `Plus`  | 310 | Left |
| `->`, `→` | `Rule` | 120 | Right |
| `:>`     | `RuleDelayed` | 120 | Right |
| `==`     | `Equal` | 290 | None |
| `===`    | `SameQ` | 290 | None |
| `=`      | `Set`   | 40 | Right |
| `:=`     | `SetDelayed` | 40 | Right |
| `=.`     | `Unset` (postfix) | 40 | None |
| `+=`     | `AddTo` | 40 | Right |
| `-=`     | `SubtractFrom` | 40 | Right |
| `*=`     | `TimesBy` | 40 | Right |
| `/=`     | `DivideBy` | 40 | Right |
| `;`      | `CompoundExpression` | 10 | Left |

The Unicode rightwards arrow `→` (U+2192, Wolfram `\[Rule]`) is accepted as a
synonym for `->`, so Wolfram-Language rules and associations paste and parse
directly.



## Parenthesised comparisons break chains

A run of comparisons folds into one variadic `Inequality` at parse time, so `a < b <= c` becomes
`Inequality[a, Less, b, LessEqual, c]` rather than a nested pair. **Parentheses stop that**, and are
honoured on both operands:

| input | parse |
|---|---|
| `a >= b == c` | `Inequality[a, GreaterEqual, b, Equal, c]` |
| `(a >= b) == c` | `Equal[GreaterEqual[a, b], c]` |
| `(a > b) == (c > d)` | `Equal[Greater[a, b], Greater[c, d]]` |
| `(a < b < c) == d` | `Equal[Inequality[a, Less, b, Less, c], d]` |

This required tracking whether the left operand came from an **unparenthesised** comparison. Deciding
by inspecting the built subtree — which is what the fold did originally — cannot distinguish the two
cases, because by the time the fold looks, `(a >= b)` and a bare `a >= b` are the same node. The result
was that parentheses were ignored on the left operand only: `(2.0 >= 2.) == (1.0 > 0.)` parsed as
`Inequality[2.0, GreaterEqual, 2.0, Equal, Greater[1.0, 0.0]]` and evaluated to `2.0 == True`. The right
operand was never affected, being parsed as its own subexpression, which is why the asymmetry went
unnoticed: every symmetric-looking test passed.

`Plus`/`Times` flattening in the same part of the parser has the identical shape and is harmless, because
those heads are associative — `(a + b) + c` and `a + b + c` mean the same thing. Chaining changes
*meaning*, so comparisons need provenance tracked rather than guessed.
