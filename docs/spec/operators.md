# 3. Operators and Precedence

Precedences are numbered on a **1–10000 ladder** (`;` = 100 at the bottom,
`[[…]]` = 10000 at the top). The levels preserve the relative ordering of
Mathematica's own scheme but are spread out with generous gaps between adjacent
levels, so new operators can be inserted without renumbering. The values below
are the internal numbers used by the parser (`get_operator` in `src/parse.c`)
and, in lockstep, by the printer's parenthesiser (`get_expr_prec` in
`src/print.c`).

| Operator | FullForm | Precedence | Association |
|----------|----------|------------|-------------|
| `[[...]]`| `Part` | 10000 | Left |
| `f[x]`   | `f[x]` | 9500 | Left |
| `_`, `__`, `___` | `Blank` | 8400 | None |
| `!`, `!!` | `Factorial`, `Factorial2` | 8200 | Left |
| `?`      | `PatternTest` | 7900 | None |
| `'`      | `Derivative` | 7700 | Left |
| `++`, `--` | `Increment`, `Decrement` | 7500 | Left |
| `@*`     | `Composition` | 7200 | Left |
| `@`      | `Prefix` | 7000 | Right |
| `@@`     | `Apply` | 7000 | Right |
| `/@`     | `Map`   | 7000 | Right |
| `<>`     | `StringJoin` | 6700 | Left |
| `^`      | `Power` | 6500 | Right |
| `.`      | `Dot` | 5300 | Left |
| `/`      | `Divide`| 5000 | Left |
| `*`      | `Times` | 4500 | Left |
| `+`, `-` | `Plus`  | 3500 | Left |
| `==`, `<`, `<=`, `;;` | `Equal`, `Less`, `LessEqual`, `Span` | 3200 | None |
| `===`    | `SameQ` | 3200 | None |
| `!` (prefix) | `Not` | 3000 | Right |
| `&&`, `\|\|` | `And`, `Or` | 2800 | Left |
| `..`, `...` | `Repeated`, `RepeatedNull` | 2500 | None |
| `\|`     | `Alternatives` | 2300 | Left |
| `~~`     | `StringExpression` | 2100 | Left |
| `:`      | `Optional` | 1900 | Right |
| `/;`     | `Condition` | 1700 | Right |
| `->`, `→` | `Rule` | 1500 | Right |
| `:>`     | `RuleDelayed` | 1500 | Right |
| `/.`, `//.` | `ReplaceAll`, `ReplaceRepeated` | 1200 | Left |
| `&`      | `Function`| 1000 | Left |
| `//`     | `Postfix` | 800 | Left |
| `=`      | `Set`   | 500 | Right |
| `:=`     | `SetDelayed` | 500 | Right |
| `=.`     | `Unset` (postfix) | 500 | None |
| `+=`, `-=`, `*=`, `/=` | `AddTo`, `SubtractFrom`, `TimesBy`, `DivideBy` | 500 | Right |
| `>>`, `>>>` | `Put`, `PutAppend` | 300 | Left |
| `;`      | `CompoundExpression` | 100 | Left |

(`_`/`__`/`___` blanks are recognised structurally in `parse_primary`, not via
the operator table; the value above is their documented rank.)

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
