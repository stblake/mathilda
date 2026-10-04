# LinearModelFit

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`LinearModelFit[data] fits a linear model with an intercept and returns a PredictorFunction carrying its coefficients. Wolfram's version returns a FittedModel with regression diagnostics (RSquared, standard errors, ANOVA); those are not implemented and are not approximated. The coefficients are the same ones Predict finds.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= m = LinearModelFit[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]
Out[1]= PredictorFunction["LinearRegression", <>]

In[2]:= {m[10.], m["Coefficients"]}
Out[2]= {21.0, {1.0, 2.0}}

In[3]:= LinearModelFit[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]
Out[3]= {1.0, 2.0, 3.0}
```

### Applications (4)

An intercept plus a slope

```mathematica
In[4]:= m = LinearModelFit[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]
Out[4]= PredictorFunction["LinearRegression", <>]
```

Apply the fitted line

```mathematica
In[5]:= m[10.]
Out[5]= 21.0
```

The same coefficients Predict finds

```mathematica
In[6]:= m["Coefficients"]
Out[6]= {1.0, 2.0}
```

Multiple regressors, last column the response

```mathematica
In[7]:= LinearModelFit[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]
Out[7]= {1.0, 2.0, 3.0}
```

## Implementation notes

**Algorithm.** `builtin_linear_model_fit` fits an ordinary least-squares line with an
intercept through the *same* `ml_ols` kernel `Predict` uses — the normal equations
`(AᵀA) c = Aᵀy` with `A = [1 | x]`, solved by Gaussian elimination with partial
pivoting — and returns a `PredictorFunction["LinearRegression", {intercept, c₁, …}, dim,
0]`. Data is read by `ml_read_training`, so a list of rules `{features -> value, …}` or
a matrix whose last column is the response both work. A singular system (collinear
features, or fewer observations than parameters) declines, since there is no unique fit.

The object is applied and queried exactly like any other `PredictorFunction`:
`m[x]` evaluates `intercept + c · x`, and `m["Coefficients"]` reads the fitted vector —
the same coefficients `Predict` finds.

**Data structures.** The shared positional, method-tagged model representation (see
`src/ml/predict.h`); the coefficient `List` is the only payload.

**Complexity / limits.** `O(n · p² + p³)` with `p = dim + 1`. Wolfram's `LinearModelFit`
returns a `FittedModel` whose properties are regression diagnostics (`RSquared`,
standard errors, ANOVA); those are a separate piece of work and are **not** approximated
here — only the coefficients, under a name a reader of regression code will look for.

**Attributes:** `Protected`.

## References

**See also:** [PredictorFunction](../../other-advanced/PredictorFunction/), [Predict](../../machine-learning/Predict/)

- T. Hastie, R. Tibshirani and J. Friedman, *The Elements of Statistical Learning*, 2nd ed. (Springer, 2009), §3.2 (least squares).
- Source: [`src/ml/predict.c`](https://github.com/stblake/mathilda/blob/main/src/ml/predict.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_predict.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_predict.c)

## Notes & additional examples

### Notes

`LinearModelFit[data]` fits a linear model with an intercept and returns a
`PredictorFunction` carrying its coefficients — the very ones `Predict` finds, through
the same least-squares kernel. It is applied and queried exactly like any other
`PredictorFunction`.

Wolfram's `LinearModelFit` returns a `FittedModel` whose properties are regression
diagnostics — `RSquared`, standard errors, ANOVA. Those are a separate piece of work and
are deliberately **not** approximated here; only the coefficients are offered, under a
name a reader of regression code will look for. A singular system returns unevaluated.
