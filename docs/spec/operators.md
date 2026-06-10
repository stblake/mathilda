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
| `==`     | `Equal` | 290 | None |
| `===`    | `SameQ` | 290 | None |
| `=`      | `Set`   | 40 | Right |
| `:=`     | `SetDelayed` | 40 | Right |
| `=.`     | `Unset` (postfix) | 40 | None |
| `;`      | `CompoundExpression` | 10 | Left |

