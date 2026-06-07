# Expressions and the evaluation model

This tutorial is about the single idea that the whole of Mathilda is built on:
**everything is an expression**. Once you see that numbers, symbols, lists, and
function calls are all the same kind of thing — a tree — the rest of the system
(pattern matching, rules, attributes, evaluation) follows naturally.

Every example below was run in the real REPL. Type each `In[...]` line yourself
(without the prompt) to follow along.

## Everything is an expression

When you type `a + b`, Mathilda does not store it as "an addition". It stores a
tree whose top node is the symbol `Plus` and whose branches are `a` and `b`.
`FullForm` shows you that internal tree, with no pretty-printing in the way:

```mathematica
In[1]:= FullForm[a + b]
Out[1]= Plus[a, b]

In[2]:= FullForm[{1, 2, 3}]
Out[2]= List[1, 2, 3]

In[3]:= FullForm[a/b]
Out[3]= Times[a, Power[b, -1]]

In[4]:= FullForm[3 x^2]
Out[4]= Times[3, Power[x, 2]]
```

Look closely at what these reveal:

- A **sum** `a + b` is really `Plus[a, b]`.
- A **list** `{1, 2, 3}` is really `List[1, 2, 3]`. Curly braces are just
  surface syntax for `List`.
- A **division** `a/b` is not a special "divide" node at all — it is
  `Times[a, Power[b, -1]]`. Mathilda has no `Divide` head internally; division
  is multiplication by a reciprocal.
- `3 x^2` is `Times[3, Power[x, 2]]` — a product of `3` and `x` raised to `2`.

So `Plus`, `Times`, `Power`, and `List` are ordinary function heads, exactly
like a function `f` you might write yourself. There is nothing privileged about
them. This uniformity is the reason a single set of tools can operate on
*anything* in the language.

## Heads: what kind of thing is this?

Every expression has a **head** — the symbol at the top of its tree. `Head`
returns it:

```mathematica
In[1]:= Head[5]
Out[1]= Integer

In[2]:= Head[{1, 2, 3}]
Out[2]= List

In[3]:= Head[a + b]
Out[3]= Plus

In[4]:= Head[x]
Out[4]= Symbol
```

Even atoms have heads: an integer's head is `Integer`, and a bare symbol's head
is `Symbol`. For a compound expression like `a + b`, the head is the function
symbol that builds it — here, `Plus`. The head is how the evaluator decides what
to do with an expression: it looks up the head's definitions and attributes.

## Structure is uniform: Length and Part

Because a sum and a list are both just expressions with arguments, the same
structural tools work on both. `Length` counts the arguments:

```mathematica
In[1]:= Length[{10, 20, 30}]
Out[1]= 3

In[2]:= Length[a + b + c]
Out[2]= 3
```

`a + b + c` has length `3` for the same reason `{10, 20, 30}` does — both are a
head applied to three arguments (`Plus[a, b, c]` and `List[10, 20, 30]`).

`Part`, written `expr[[i]]`, extracts the *i*-th argument, and again it does not
care whether you hand it a list or a sum:

```mathematica
In[1]:= {10, 20, 30}[[2]]
Out[1]= 20

In[2]:= (a + b + c)[[2]]
Out[2]= b
```

The second argument of the sum is `b`, just as the second element of the list is
`20`. One mechanism, every expression.

## The evaluation model: rewriting to a fixed point

Mathilda evaluates by **repeated rewriting**. It applies every rule it can,
producing a new expression, then evaluates *that*, and keeps going until nothing
changes — a "fixed point". You never have to chain steps together by hand; the
loop does it for you.

Define two functions and then call the first:

```mathematica
In[1]:= f[x_] := g[x]
Out[1]= Null

In[2]:= g[x_] := x^2
Out[2]= Null

In[3]:= f[3]
Out[3]= 9
```

`f[3]` is not the answer — it rewrites to `g[3]`, which rewrites to `3^2`, which
evaluates to `9`. All three rewrites happen in a single call because evaluation
runs to a fixed point automatically. (The `Out= Null` lines are just the result
of making a definition — assignments return nothing.)

The same definitions keep working for any input:

```mathematica
In[1]:= f[5]
Out[1]= 25
```

## Attributes drive evaluation

How does the evaluator know to flatten nested sums, or to sort the arguments of
a product? It doesn't have that logic hard-coded per function. Instead each
symbol carries **attributes** — flags the evaluator consults before it processes
a call. Ask for them with `Attributes`:

```mathematica
In[1]:= Attributes[Plus]
Out[1]= {Flat, Listable, NumericFunction, OneIdentity, Orderless, Protected}
```

Three of these directly change how expressions are evaluated:

### Orderless — canonical ordering

`Orderless` means the arguments are sorted into a canonical order, so
commutative operations look the same no matter how you wrote them. That is why
`b + a` comes back reordered:

```mathematica
In[1]:= b + a
Out[1]= a + b
```

You typed the `b` first, but `Plus` is `Orderless`, so Mathilda sorts the
arguments. This is not cosmetic — it means `a + b` and `b + a` are literally the
same expression internally, which is exactly what pattern matching needs.

### Flat — associativity and flattening

`Flat` means nested calls with the same head are flattened into one. Writing a
sum with explicit parentheses produces a single flat `Plus`, not a `Plus`
containing another `Plus`:

```mathematica
In[1]:= FullForm[a + (b + c)]
Out[1]= Plus[a, b, c]
```

The inner `(b + c)` is absorbed: `a + (b + c)` becomes `Plus[a, b, c]`, a
three-argument sum. The same applies to `Times`.

### Listable — threading over lists

`Listable` means the function automatically maps over the elements of a list
argument. `Sqrt` is `Listable`, so applying it to a list applies it to each
element:

```mathematica
In[1]:= Sqrt[{1, 4, 9}]
Out[1]= {1, 2, 3}
```

You did not write a loop or a `Map`. Because `Sqrt` is `Listable`, the
evaluator threaded it over the list for you, giving the square root of each
element.

## Holding evaluation

Sometimes you want an expression *not* to evaluate — to keep it as a structure
to inspect or transform. `Hold` does exactly that: it wraps its argument and
prevents it from being evaluated.

```mathematica
In[1]:= Hold[1 + 1]
Out[1]= Hold[1 + 1]
```

The `1 + 1` inside is left untouched. To let it evaluate again, strip the wrapper
with `ReleaseHold`:

```mathematica
In[1]:= ReleaseHold[Hold[1 + 1]]
Out[1]= 2
```

`Hold` works because of an attribute called `HoldAll`, which tells the evaluator
to leave a function's arguments *unevaluated* before the function itself runs:

```mathematica
In[1]:= Attributes[Hold]
Out[1]= {HoldAll, Protected}
```

`HoldAll` is the same mechanism that lets control-flow and definition constructs
receive their arguments raw — for instance, when you write `f[x_] := g[x]`, the
right-hand side must *not* be evaluated at definition time, and a hold attribute
is what makes that possible.

## Where to next

You now have the mental model that the rest of Mathilda depends on: expressions
are trees, every tree has a head, structure is uniform, evaluation rewrites to a
fixed point, and attributes plus holding control how that happens.

- **[3. Pattern matching & rules](03-pattern-matching-and-rules.md)** — now that
  you know expressions are trees, learn how to match against their shape with
  blanks (`_`, `x_`), apply transformation rules (`->`, `:>`, `/.`, `//.`), and
  define your own functions.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions, organised by category, with details on `Head`, `Part`,
  `Attributes`, `Hold`, and everything else used above.
