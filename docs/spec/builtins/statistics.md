# Statistics


## Median
Gives the median estimate of the elements in data.
- `Median[data]`: gives the median estimate $\hat{q}_{1/4}$ of the elements in `data`.
- `Median[dist]`: gives the median of the distribution `dist`.

**Features**:
- `Protected`.
- Median is a robust location estimator, which means it not very sensitive to outliers.
- For `VectorQ` data $\{x_1, \dots, x_n\}$, the median can be thought of as the "middle value". Formally, when data is sorted as $\{x_{(1)}, \dots, x_{(n)}\}$, the median is given by the center element $x_{((n+1)/2)}$ if $n$ is odd and the mean of the two center elements $(x_{(n/2)} + x_{(n/2+1)})/2$ if $n$ is even.
- For `MatrixQ` data, the median is computed for each column vector. `Median` for a tensor gives columnwise medians at the first level.
- `Median` requires numeric values.

```mathematica
In[1]:= Median[{1, 2, 3, 4, 5, 6, 7}]
Out[1]= 4

In[2]:= Median[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[2]= 9/2

In[3]:= Median[{1, 2, 3, 4}]
Out[3]= 5/2

In[4]:= Median[{Pi, E, 2}]
Out[4]= E

In[5]:= Median[{1., 2., 3., 4.}]
Out[5]= 2.5

In[6]:= Median[{{1, 11, 3}, {4, 6, 7}}]
Out[6]= {5/2, 17/2, 5}

In[7]:= Median[{{{3, 7}, {2, 1}}, {{5, 19}, {12, 4}}}]
Out[7]= {{4, 13}, {7, 5/2}}

In[8]:= Median[{a, b, c}]
Median::rectn: Rectangular array of real numbers is expected at position 1 in Median[{a, b, c}].
Out[8]= Median[{a, b, c}]
```
## Mean
Gives the mean estimate of the elements in data.
- `Mean[data]`

**Features**:
- `Protected`.
- Supports numerical and symbolic data.
- For vectors, computes $(1/n) \sum x_i$.
- For matrices, computes means of elements in each column.

```mathematica
In[1]:= Mean[{1, 2, 3, 4}]
Out[1]= 5/2

In[2]:= Mean[{{a, u}, {b, v}, {c, w}}]
Out[2]= {1/3 (a + b + c), 1/3 (u + v + w)}
```

## RootMeanSquare
Gives the root mean square of values in `list`.
- `RootMeanSquare[list]`

**Features**:
- `Protected`.
- Gives the square root of the second sample moment.
- For a list `{x1, x2, ...}`, it computes `Sqrt[1/n Total[{x1^2, x2^2, ...}]]`.
- Handles both numerical and symbolic data.
- Works column-wise on matrices.

```mathematica
In[1]:= RootMeanSquare[{a, b, c, d}]
Out[1]= 1/2 Sqrt[a^2 + b^2 + c^2 + d^2]

In[2]:= RootMeanSquare[{{1, 2}, {5, 10}, {5, 2}, {4, 8}}]
Out[2]= {1/2 Sqrt[67], Sqrt[43]}

In[3]:= RootMeanSquare[{1, 2, 3, 4}]
Out[3]= Sqrt[15/2]

In[4]:= RootMeanSquare[{Pi, E, 2}]
Out[4]= Sqrt[1/3 (4 + E^2 + Pi^2)]

In[5]:= RootMeanSquare[{1., 2., 3., 4.}]
Out[5]= 2.73861
```

## Variance
Gives the unbiased variance estimate of the elements in data.
- `Variance[data]`

**Features**:
- `Protected`.
- For vectors, computes $(1/(n-1)) \sum (x_i - \hat{\mu}) \overline{(x_i - \hat{\mu})}$.
- For matrices, computes variances of elements in each column.

```mathematica
In[1]:= Variance[{1, 2, 3}]
Out[1]= 1

In[2]:= Variance[{{5.2, 7}, {5.3, 8}, {5.4, 9}}]
Out[2]= {0.01, 1}
```

## StandardDeviation
Gives the standard deviation estimate of the elements in data.
- `StandardDeviation[data]`

**Features**:
- `Protected`.
- Equivalent to `Sqrt[Variance[data]]`.
- For matrices, computes standard deviations of elements in each column.

```mathematica
In[1]:= StandardDeviation[{1, 2, 3}]
Out[1]= 1
```

## MovingAverage
Gives the moving average over a list, with either a uniform window length or a list of weights.
- `MovingAverage[list, r]`: averages runs of `r` consecutive elements.
- `MovingAverage[list, {w_1, w_2, ..., w_r}]`: weighted moving average with effective weights $w_i / \sum_j w_j$.

**Features**:
- `Protected`.
- Output length is `Length[list] - r + 1`.
- Stays unevaluated when `r < 1`, when `r > Length[list]`, when the second argument is non-integer / non-list, or when the first argument is not a `List`.
- Exact rational arithmetic for integer / rational data; bignums (arbitrary-precision integers) handled natively. Real-valued data or weights yield approximate output. Symbolic data and weights are supported.
- The unweighted form delegates to `Mean` for each window, so it inherits `Mean`'s exact / numeric / symbolic dispatch.

```mathematica
In[1]:= MovingAverage[{1, 5, 7, 3, 6, 2}, 3]
Out[1]= {13/3, 5, 16/3, 11/3}

In[2]:= MovingAverage[{1.2, 5.2, 3.4, 4.5, 2.3, 4.5}, 3]
Out[2]= {3.26667, 4.36667, 3.4, 3.76667}

In[3]:= MovingAverage[{a, b, c, d, e}, 2]
Out[3]= {1/2 (a + b), 1/2 (b + c), 1/2 (c + d), 1/2 (d + e)}

In[4]:= MovingAverage[{a, b, c, d, e}, {1, 2}]
Out[4]= {1/3 a + 2/3 b, 1/3 b + 2/3 c, 1/3 c + 2/3 d, 1/3 d + 2/3 e}

In[5]:= MovingAverage[{2^100, 2^101, 2^102, 2^103}, 2]
Out[5]= {1901475900342344102245054808064, 3802951800684688204490109616128, 7605903601369376408980219232256}

In[6]:= MovingAverage[{1, 2, 3, 4, 5}, 6]
Out[6]= MovingAverage[{1, 2, 3, 4, 5}, 6]
```

## MovingMedian
Gives the moving median of a list, taken over spans of `r` consecutive elements.
- `MovingMedian[list, r]`: medians of runs of `r` elements.

**Features**:
- `Protected`.
- Output length is `Length[list] - r + 1`.
- Operates on real-valued vectors and matrices. For matrix input, each window of `r` consecutive rows is reduced via `Median`, yielding a column-wise median vector per window.
- Exact rationals, bignums (arbitrary-precision integers), machine-precision reals, and `NumericQ`-real symbolic constants (`Pi`, `E`, ...) are all supported. Even-window medians yield exact rational midpoints when the data is exact.
- Stays unevaluated when `r < 1`, when `r > Length[list]`, when `r` is non-integer, or when the first argument is not a `List`.
- Non-numeric data triggers the `MovingMedian::arg1` message and the expression remains unevaluated.

```mathematica
In[1]:= MovingMedian[{1, 2, 5, 6, 1, 4, 3}, 3]
Out[1]= {2, 5, 5, 4, 3}

In[2]:= MovingMedian[{{1, 2}, {5, 3}, {1, 4}, {3, 2}, {5, 5}}, 2]
Out[2]= {{3, 5/2}, {3, 7/2}, {2, 3}, {4, 7/2}}

In[3]:= MovingMedian[N[{1, 5, 7, 3, 6, 2}], 3]
Out[3]= {5.0, 5.0, 6.0, 3.0}

In[4]:= MovingMedian[{1, 2, 3, 4}, 2]
Out[4]= {3/2, 5/2, 7/2}

In[5]:= MovingMedian[{2^100, 2^101, 2^102, 2^103}, 2]
Out[5]= {1901475900342344102245054808064, 3802951800684688204490109616128, 7605903601369376408980219232256}

In[6]:= MovingMedian[{a, b, c}, 2]
MovingMedian::arg1: The first argument {a, b, c} must be a vector or matrix of real values.
Out[6]= MovingMedian[{a, b, c}, 2]
```

## ExponentialMovingAverage
Gives the exponential moving average of a list with smoothing constant `alpha`.
- `ExponentialMovingAverage[list, alpha]`: produces the recurrence $y_1 = x_1$, $y_{i+1} = y_i + \alpha (x_{i+1} - y_i)$.

**Features**:
- `Protected`.
- Output has the same length as `list`.
- Two evaluation strategies: a fast O(n) double-precision path activates when at least one element of `list` or `alpha` itself is a machine-precision real and every other entry is a real-valued numeric (Integer, Real, Rational); otherwise a symbolic recurrence path is taken so exact rationals, bignums (arbitrary-precision integers), and symbolic alpha all work.
- The smoothing constant `alpha` is typically a number between 0 and 1 but may be any expression; with `alpha = 0` the output is constant at `x_1`, and with `alpha = 1` the output equals the input.
- Stays unevaluated when the first argument is not a `List`, when the list is empty, or when the call has the wrong arity.

```mathematica
In[1]:= ExponentialMovingAverage[Range[10], 1/3]
Out[1]= {1, 4/3, 17/9, 70/27, 275/81, 1036/243, 3773/729, 13378/2187, 46439/6561, 158488/19683}

In[2]:= ExponentialMovingAverage[N[{1, 5, 7, 3, 6, 2}], 1/2]
Out[2]= {1.0, 3.0, 5.0, 4.0, 5.0, 3.5}

In[3]:= ExponentialMovingAverage[{a, b, c, d}, 0]
Out[3]= {a, a, a, a}

In[4]:= ExponentialMovingAverage[{a, b, c, d}, 1]
Out[4]= {a, b, c, d}

In[5]:= ExponentialMovingAverage[{a, b}, x]
Out[5]= {a, a + x (-a + b)}

In[6]:= ExponentialMovingAverage[{2^100, 2^200}, 1/2]
Out[6]= {1267650600228229401496703205376, 803469022129495137770981046171215126561215611592144769253376}
```

