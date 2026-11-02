# 3. Operators and Precedence

Precedences are numbered on a **1–10000 ladder** (`;` = 100 at the bottom,
`[[…]]` = 10000 at the top). The levels preserve the relative ordering of
the reference's own scheme but are spread out with generous gaps between adjacent
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
| `**`     | `NonCommutativeMultiply` | 5900 | Left |
| `.`      | `Dot` | 5300 | Left |
| `/`      | `Divide`| 5000 | Left |
| `*`      | `Times` | 4500 | Left |
| `+`, `-` | `Plus`  | 3500 | Left |
| `+` (prefix) | `Plus` (one argument: `+x` is `Plus[x]`, which evaluates to `x`) | 3500 | Right |
| `==`, `<`, `<=`, `;;` | `Equal`, `Less`, `LessEqual`, `Span` | 3200 | None |
| `===`    | `SameQ` | 3200 | None |
| `!` (prefix) | `Not` | 3000 | Right |
| `&&`     | `And` | 2800 | Left |
| `\|\|`   | `Or` | 2700 | Left |
| `..`, `...` | `Repeated`, `RepeatedNull` (postfix; a single token that wins over a decimal point, so `0..` is `Repeated[0]`) | 2500 | None |
| `\|`     | `Alternatives` | 2300 | Left |
| `~~`     | `StringExpression` | 2100 | Left |
| `:`      | `Pattern` (`sym:p`) / `Optional` (`p:v`) — see below | 2200 / 2150 | Left |
| `/;`     | `Condition` | 1700 | Left (`a /; b /; c` is `Condition[Condition[a, b], c]`) |
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
| `<<` (prefix) | `Get` | tightest | — |
| `;`      | `CompoundExpression` | 100 | Left |

(`_`/`__`/`___` blanks are recognised structurally in `parse_primary`, not via
the operator table; the value above is their documented rank.)

**`:` binds its left operand at atom level.** As in the reference, the
left operand of `:` is the single token right before it: a bare symbol gives
`Pattern` (`a:p`), a blank token gives `Optional` (`x_:v`, `_h:v`). It is never
the expression built so far. The right operand is read at the Pattern (2200) or
Optional (2150) rank, so it takes in everything that binds tighter, `|`, `..`,
`==` and `+` included. One further `:` after `sym:p` supplies a default:

| Input | FullForm |
|-------|----------|
| `n_.*a:_g` | `Times[Optional[Pattern[n, Blank[]]], Pattern[a, Blank[g]]]` |
| `n x_:0` | `Times[n, Optional[Pattern[x, Blank[]], 0]]` |
| `x + a:_g == b` | `Plus[x, Pattern[a, Equal[Blank[g], b]]]` |
| `x?f:_` | `PatternTest[x, Pattern[f, Blank[]]]` |
| `x|a:b` | `Alternatives[x, Pattern[a, b]]` |
| `a:_g:1` | `Optional[Pattern[a, Blank[g]], 1]` |
| `x:a:b:c:d:e` | `Optional[Pattern[x, a], Optional[Pattern[b, c], Pattern[d, e]]]` |
| `a:_g -> 1`, `a:_ /; c` | `Rule[Pattern[a, Blank[g]], 1]`, `Condition[Pattern[a, Blank[]], c]` |

Before v0.475 `:` was an ordinary low-precedence infix operator, so `n_.*a:_g`
parsed as `Optional[n_.*a, _g]`, which is a silently wrong pattern (issue #102).
A left operand that the reference rejects (`(a):_g`, `f[x]:0`) is still accepted, as
`Pattern`/`Optional` of the whole left expression.

The raw Unicode rightwards arrow `→` (U+2192) is accepted as a synonym for `->`, so rules and
associations pasted from rendered text parse directly. This is a deliberate deviation: the
reference's `\[Rule]` is the private-use U+F522, and it reads U+2192 (`\[RightArrow]`) as
`RightArrow[a, b]`. `\[Rule]` itself, raw or escaped, is `Rule`.

## Named characters in source

A non-ASCII character in source code (outside a string literal) is read as the reference
reads it. It can be written raw (`μ`), as a named character (`\[Mu]`), or as a hex escape
(`\:03bc`, `\|0003bc`); all spellings give the same token, and a symbol name stores the
UTF-8 character. The class of every named character comes from the reference parser
(`tools/classify_named_chars.wls` writes the table, which `tools/gen_named_chars.py --classes`
turns into `src/strings/named_char_class.inc`):

| class | examples | reads as |
|---|---|---|
| letter | `\[Mu]`, `\[CapitalSigma]`, `\[CurlyEpsilon]`, `\[Dagger]`, `\[Prime]` | part of an identifier: `a\[Mu]1` and `\[Mu]$` are each one symbol, `\[Alpha] \[Beta]` is a product |
| constant letter | `\[Pi]`, `\[ImaginaryI]`, `\[ExponentialE]`, `\[Infinity]`, `\[Degree]` | alone, the constant (`2\[Pi]` is `Times[2, Pi]`); inside a name, a letter (`a\[Pi]` is one symbol) |
| operator with an ASCII form | `\[Rule]`, `\[RuleDelayed]`, `\[Equal]`, `\[NotEqual]`, `\[LessEqual]`, `\[Times]`, `\[And]`, `\[Not]` | exactly that ASCII operator (`a < b \[LessEqual] c` is one `Inequality`) |
| infix operator | `\[CenterDot]`, `\[Cross]`, `\[Element]`, `\[LeftRightArrow]`, `\[Perpendicular]`, `\[Function]`, `\[Union]`, ... (about 230) | its head, n-ary when the reference groups it so (`a\[CenterDot]b\[CenterDot]c` is `CenterDot[a, b, c]`, while `(a\[CenterDot]b)\[CenterDot]c` nests), otherwise left- or right-associative |
| prefix / postfix operator | `\[Del]`, `\[Sqrt]`, `\[ForAll]` / `\[Transpose]`, `\[Conjugate]` | its head applied to the operand |
| whitespace | `\[InvisibleSpace]`, `\[VeryThinSpace]`, `\[NoBreak]`, `\[MediumSpace]` | a space (`a\[InvisibleSpace]b` is `Times[a, b]`) |
| line break | `\[IndentingNewLine]`, `\[NewLine]` | a newline, which ends a top-level statement as `\n` does |
| not source | `\[PartialD]`, brackets such as `\[LeftDoubleBracket]` | a syntax error |

Operator precedences come from the reference too and are mapped onto the 1–10000 ladder
between the operators both share (`\[CenterDot]` is 4571, just tighter than `*` at 4500).
A character missing from the table, such as an accented or CJK letter, is a letter.

Known differences: `FullForm` prints a non-ASCII symbol as the character itself, where the
reference prints `\[Mu]`. Box input (`\!\(\*SuperscriptBox[...]\)`) is not read.



## Linear box syntax: `\!\( ... \)`

`\( ... \)` groups box text, `\*` introduces a box expression, and an
optional `\!` in front asks for the boxes to be read as the expression they
display. FeynCalc's tests write powers and fractions this way:
`\!\(\*SuperscriptBox[\(mb\), \(2\)]\) - x` is `mb^2 - x`.

- Mathilda has no box language, so a group is lowered to InputForm text and
  parsed (`src/linear_boxes.c`).
  - Plain text is copied, and a nested `\( X \)` becomes `( X )`.
  - `\ ` (a backslash and a space) is a space token in box text and is lowered
    to a space (v0.6139): `\(3  p^\ 2\)` reads as `3 p^2`, as in the reference.
  - `SuperscriptBox[a, b]` becomes `a^b`, `SubscriptBox` becomes
    `Subscript`, `FractionBox[a, b]` becomes `a/b`, and `SqrtBox[a]`
    becomes `Sqrt[a]`. `OverscriptBox`, `UnderscriptBox` and
    `SubsuperscriptBox` map to their heads.
  - `RowBox[{...}]` concatenates its elements, `FormBox`, `StyleBox`,
    `TagBox` and `AdjustmentBox` show their first argument, and
    `InterpretationBox[b, expr]` stands for `expr`.
- Juxtaposed groups multiply, as in the reference.
- In a string literal, `\!`, `\(`, `\)`, `\*`, `\@`, `\^`, `\&`, `\%`, `\+`,
  `\/`, `\_` and ``\` `` are the reference's private-use characters (U+F7C1,
  U+F7C9, U+F7C0, U+F7C8, ...).
  - So `StringLength["\(\)"]` is `2`, and `ToExpression` of such a string
    reads the boxes.
  - Every form prints them back as their escapes.

## Number literals: `e` exponents

The reference writes an exponent only as `*^` (`1.5*^-3`). Mathilda also
accepts `e`/`E` (`1.5e-3`, `2e5`) as an extension, but only when digits follow
it: `2ep` is `2 ep` and `3E` is `3 E`, as the reference reads them (it used to
make `2ep` the Real `2.` times `p`, v0.591).

## Flat operator chains

`+`, `*` (and implicit multiplication), `**` and `<>` build one n-ary node per
chain: `a b c` is `Times[a, b, c]`. Only the chain itself is flattened: a
parenthesised or explicitly written operand stays one argument, as in the
reference (`(a b) c` and `Times[a, b] c` are `Times[Times[a, b], c]`,
`(a + b) - c` is `Plus[Plus[a, b], Times[-1, c]]`).

Prefix minus binds tighter than `*` and is one operand too: `-a b` is
`Times[Times[-1, a], b]` and `-(x + y) z` is `Times[Times[-1, Plus[x, y]], z]`,
so the negated sum distributes before the product forms. Binary minus puts its
`-1` into a bare product: `a - b c` is `Plus[a, Times[-1, b, c]]`, while
`a - (b c)` is `Plus[a, Times[-1, Times[b, c]]]`.

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

## Statements in a file (`Get`, `-file`)
A file is read one statement at a time, and each statement is evaluated before
the next is parsed. That way a `BeginPackage`/`Begin` prologue governs the
symbols created by later statements.
- A top-level `;` ends a statement, with one exception. If the next token on
  the same line can only continue an expression, the `;` is the
  `CompoundExpression` operator with an implicit `Null`, and that compound is
  the operator's left operand.
- The tokens that continue a statement are the binary and postfix operators:
  `/;`, `//`, `/.`, `//.`, `->`, `:>`, `:=`, `=`, `==`, `===`, `&&`, `||`, `&`,
  `@@`, `/@`, `@`, `^`, `*`, `<`, `<=`, `>`, `>=`, `<>`, `~~`, `|`, `::`, `-=`
  and `+=`.
- So `f[x_] := g[x];/; x > 0` reads as
  `Condition[CompoundExpression[f[x_] := g[x], Null], x > 0]`. Evaluating it
  defines `f` unconditionally.
- An operator at the start of the next line is a syntax error.
- The same rule holds for `;` inside one expression (`ToExpression`, a
  parenthesised body): `(a; b; -> c)` is
  `Rule[CompoundExpression[a, b, Null], c]`. `;` associates left, so the whole
  chain is built before the operator binds, and a parenthesised compound stays
  nested: `(a; b); c` is `CompoundExpression[CompoundExpression[a, b], c]`.
- On a syntax error, `Get` emits `Syntax::sntx` ("Invalid syntax in or before
  "...` (line n of "file")"), stops reading, and returns `Null`. The message
  goes through the funnel, so `Quiet` and `Check` see it. An expression left
  unfinished at end of file counts as a syntax error.
