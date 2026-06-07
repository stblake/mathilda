# Pattern matching and rules

In the previous tutorial you saw that everything in Mathilda is an expression —
a tree with a head and arguments. **Pattern matching** is how you describe the
*shape* of such a tree without writing it out in full, and **rules** are how you
rewrite one shape into another. Together they are the engine behind almost
everything symbolic Mathilda does: simplification, differentiation, solving, and
the functions you define yourself.

Every example below was run in the real REPL. Type each `In[...]` line yourself
(without the prompt) to follow along.

## Blanks: matching anything

The most basic pattern is the **blank**, written `_`. It matches any single
expression. The function `MatchQ[expr, pattern]` returns `True` or `False`
depending on whether `expr` fits `pattern`:

```mathematica
In[1]:= MatchQ[5, _]
Out[1]= True

In[2]:= MatchQ[5, _Integer]
Out[2]= True

In[3]:= MatchQ[x, _Integer]
Out[3]= False

In[4]:= MatchQ[3.5, _Real]
Out[4]= True
```

A bare `_` matches everything, so `In[1]` is `True`. Writing a symbol *after* the
underscore — `_Integer`, `_Real` — restricts the match to expressions with that
**head**. `5` has head `Integer`, so `In[2]` matches; the symbol `x` has head
`Symbol`, not `Integer`, so `In[3]` does not.

This head test works for any head at all, including your own functions and
structural heads like `List`:

```mathematica
In[1]:= MatchQ[f[a], _f]
Out[1]= True

In[2]:= MatchQ[{1, 2}, _List]
Out[2]= True

In[3]:= MatchQ[3, _Symbol]
Out[3]= False
```

Recall from tutorial 2 that `{1, 2}` is really `List[1, 2]`, so `_List` matches
it. `3` has head `Integer`, not `Symbol`, so `In[3]` is `False`. Whenever you
want to know "what would match this pattern?", reach for `MatchQ` first.

## Named patterns: capturing what matched

A blank on its own only answers yes/no. To *use* the thing that matched, give the
blank a name: `x_` means "match anything, and call it `x`". Inside a rule, the
name then stands for whatever was captured. The replacement operator `/.`
(read "slash-dot", short for `ReplaceAll`) applies a rule to an expression:

```mathematica
In[1]:= f[a, b] /. f[x_, y_] -> x + y
Out[1]= a + b

In[2]:= {2, 4, 6} /. x_ -> x^2
Out[2]= {4, 16, 36}
```

In `In[1]`, the pattern `f[x_, y_]` matches `f[a, b]`, binding `x` to `a` and
`y` to `b`; the right-hand side `x + y` then becomes `a + b`. In `In[2]`, the
pattern `x_` matches each element of the list in turn — `/.` walks the whole
expression looking for matches — so every element is squared.

Two underscores with the same name must match the **same** value. This is what
makes patterns powerful: `f[x_, x_]` matches `f[a, a]` but not `f[a, b]`.

## Sequences: `__` and `___`

A single blank matches exactly one expression. To match a *run* of arguments,
use a **double blank** `__` (`BlankSequence`, one or more) or a **triple blank**
`___` (`BlankNullSequence`, zero or more):

```mathematica
In[1]:= f[a, b, c] /. f[x__] -> {x}
Out[1]= {a, b, c}

In[2]:= h[] /. h[x__] -> matched
Out[2]= h[]

In[3]:= h[] /. h[x___] -> matched
Out[3]= matched
```

In `In[1]`, `x__` swallows all three arguments of `f` as a sequence, which the
right-hand side wraps into a list. The difference between `__` and `___` shows up
when there are *no* arguments: `h[]` has nothing for `x__` to bind to, so the
match in `In[2]` fails and the expression is returned unchanged; but `x___`
happily matches the empty sequence, so `In[3]` succeeds. Use `___` whenever an
"optional, possibly-empty" run of arguments is what you mean.

## Conditions and tests

Often you want a pattern to match only when some extra condition holds. There are
two ways to say this.

A **condition**, written `pattern /; test`, matches only when `test` evaluates to
`True`. A **pattern test**, written `pattern?predicate`, matches only when
`predicate[match]` is `True`. They overlap, but `/;` lets you write an arbitrary
boolean expression while `?` is a compact way to apply a single predicate:

```mathematica
In[1]:= {1, 2, 3, 4} /. x_ /; EvenQ[x] -> 0
Out[1]= {1, 0, 3, 0}

In[2]:= {1, 2, 3, 4} /. x_?EvenQ -> 0
Out[2]= {1, 0, 3, 0}

In[3]:= {1, 2, 3, 4, 5} /. x_ /; x > 3 -> big
Out[3]= {1, 2, 3, big, big}
```

`In[1]` and `In[2]` say the same thing two ways: replace every even element with
`0`. The odd elements fail the test, so the rule leaves them alone. `In[3]` uses
a comparison directly in the condition — `x > 3` — to replace only the elements
larger than three. The pattern variable `x` is available inside both `/;` and
`?`, which is what lets the test inspect the matched value.

## Transformation rules

A rule has two halves: a left-hand pattern and a right-hand replacement. Mathilda
has two kinds.

- `lhs -> rhs` (`Rule`, the arrow) evaluates the right-hand side **immediately**,
  when the rule is *built*.
- `lhs :> rhs` (`RuleDelayed`, colon-arrow) holds the right-hand side and
  evaluates it **each time the rule is applied**, after the pattern variables are
  bound.

For most rules the two behave identically, because the interesting work happens
through the bound pattern variables. The difference matters when the right-hand
side depends on something that can *change* between building the rule and
applying it:

```mathematica
In[1]:= val := dynamic
Out[1]= Null

In[2]:= dynamic = 100
Out[2]= 100

In[3]:= r = (x_ -> val)
Out[3]= x_ -> 100

In[4]:= dynamic = 200
Out[4]= 200

In[5]:= {a} /. r
Out[5]= 100

In[6]:= {a} /. (x_ :> val)
Out[6]= 200
```

Watch `Out[3]`: building the `->` rule evaluated `val` right then and froze it at
`100`, so even after `dynamic` becomes `200` the stored rule `r` still produces
`100` (`Out[5]`). The `:>` rule in `In[6]` leaves `val` unevaluated until it is
applied, so it sees the current value `200`. A good rule of thumb: use `:>`
whenever the right-hand side should be recomputed per match (it almost always
should when it contains a pattern variable inside a function call), and reach for
`->` only when you specifically want the value captured once.

## One pass versus a fixed point

`/.` makes a **single top-down pass**: once a subexpression is rewritten, Mathilda
does not re-scan the result for further matches. Its repeating cousin `//.`
(`ReplaceRepeated`) applies the rule **over and over until nothing changes** —
that is, to a fixed point:

```mathematica
In[1]:= {x, x^2, x^3} /. x -> 2
Out[1]= {2, 4, 8}

In[2]:= x + y /. {x -> 1, y -> 2}
Out[2]= 3
```

`In[1]` shows a plain symbol rule (`x -> 2`, no blank needed) replacing every `x`.
`In[2]` shows that the right operand of `/.` can be a *list* of rules applied
together. Now contrast the two replacement operators on a nested expression:

```mathematica
In[1]:= f[f[f[x]]] /. f[a_] -> a
Out[1]= f[f[x]]

In[2]:= f[f[f[x]]] //. f[a_] -> a
Out[2]= x
```

The rule `f[a_] -> a` strips one layer of `f`. With `/.`, the outermost `f` is
rewritten to `f[f[x]]` and the single pass stops — `Out[1]` still has two `f`s.
With `//.`, Mathilda keeps re-applying the rule to the result until no `f`
remains, peeling all the way down to `x`. Use `//.` whenever a rewrite can expose
new opportunities for the same rule.

## Defining your own functions

A definition like `square[x_] := x^2` is nothing more than a stored rule attached
to the symbol `square` (its **DownValues**, in the language of tutorial 2).
Mathilda consults these rules every time it evaluates a call to `square`. Note
the delayed assignment `:=` — it holds the body until the function is actually
called:

```mathematica
In[1]:= square[x_] := x^2
Out[1]= Null

In[2]:= square[7]
Out[2]= 49

In[3]:= add[x_, y_] := x + y
Out[3]= Null

In[4]:= add[3, 10]
Out[4]= 13
```

Defining a function returns `Null` (there is no result to show — you are storing
a rule, not computing a value), which is why `Out[1]` and `Out[3]` are blank-ish.
The calls in `In[2]` and `In[4]` then match the stored patterns and compute.

Because definitions are just rules and Mathilda evaluates to a fixed point,
**recursion works with no extra machinery**. Give a base case and a recursive
case, and the evaluator chains them automatically:

```mathematica
In[1]:= fac[0] = 1
Out[1]= 1

In[2]:= fac[n_] := n fac[n - 1]
Out[2]= Null

In[3]:= fac[5]
Out[3]= 120

In[4]:= fac[10]
Out[4]= 3628800
```

Here `fac[0] = 1` is an *immediate* assignment (`=`) — there is nothing to delay,
so it stores `1` directly and even echoes it as `Out[1]`. The recursive case uses
`:=`. When you call `fac[5]`, the recursive rule fires repeatedly until the
argument reaches `0`, at which point the base case stops the recursion. Mathilda
tries the more specific rule (`fac[0]`) before the general one (`fac[n_]`), so the
base case always wins when it applies.

## Pulling matches out of a list

Patterns are not only for rewriting — they are also a query language. `Cases`
returns every element of a list that matches a pattern, and `Count` returns how
many there are:

```mathematica
In[1]:= Cases[{1, 2, 3, 4, 5}, _?EvenQ]
Out[1]= {2, 4}

In[2]:= Count[{1, 2, 3, 4, 5}, _?EvenQ]
Out[2]= 2

In[3]:= Cases[{1, a, 2, b, 3}, _Integer]
Out[3]= {1, 2, 3}

In[4]:= Count[{1, a, 2, b, 3}, _Symbol]
Out[4]= 2
```

`In[1]` keeps just the even numbers, using the same `_?EvenQ` pattern test you met
earlier. `In[3]` filters by head — `_Integer` selects the integers and drops the
symbols `a` and `b` — while `In[4]` counts the symbols instead. Anything you can
express as a pattern, you can search a list for. This pairing of "match to
rewrite" (`/.`, `//.`) and "match to query" (`Cases`, `Count`) covers a huge
fraction of day-to-day symbolic work.

## Where to next

You can now describe the shape of any expression with blanks, capture pieces with
named patterns, guard matches with conditions and tests, rewrite with rules, and
define your own (even recursive) functions. This is the heart of how Mathilda
works — and how *you* extend it.

- **[4. Machine & arbitrary precision](04-machine-and-arbitrary-precision-arithmetic.md)**
  — meet Mathilda's three kinds of numbers: exact integers and rationals, fast
  machine-precision reals, and arbitrary-precision arithmetic with `N`.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions, including the complete details of `MatchQ`, `Cases`,
  `Count`, `Replace`, and the rule operators used above.
