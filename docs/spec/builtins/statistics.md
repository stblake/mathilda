# Statistics

**NDArray fast paths.** `Mean`, `Median`, `Variance`, `Moment`, `CentralMoment`,
`Skewness`, `Kurtosis`, `StandardDeviation`, `RootMeanSquare`, `Quartiles`,
`MovingAverage`, `MovingMedian`, `ExponentialMovingAverage`, `Covariance` and
`Correlation` operate directly on an `NDArray`'s flat buffer at C speed,
returning a scalar or a lower-rank `NDArray` (matrix reductions are columnwise).
`Covariance` / `Correlation` reduce two vectors to a scalar with a threaded
centered inner product and the matrix forms to a `p×q` matrix with a BLAS gram
($A_c^\top B_c$ via `cblas_dsyrk` / `cblas_dgemm`); they are also lowered inside
`Compile[]`. Cases outside the fast domain (complex order statistics, weighted
moving-average specs, integer/complex covariance, …) fall back to the exact
`List` result. See `src/ndreduce.c` and `src/linalg/ndcorrcov.c`.


## Median

- `Median[dist]` of a `NormalDistribution`, `UniformDistribution`,
  `ExponentialDistribution` (`Log[2]/lambda`) or 1-D `SmoothKernelDistribution`:
  its closed form (see machine-learning.md, "Distributions").
- `Median[data]`: gives the median estimate $\hat{q}_{1/2}$ of the elements in `data`.
- `Median[dist]`: gives the median of the distribution `dist`.

**Features**:
- `Protected`.
- Median is a robust location estimator, which means it not very sensitive to outliers.
- For `VectorQ` data $\{x_1, \dots, x_n\}$, the median can be thought of as the "middle value". Formally, when data is sorted as $\{x_{(1)}, \dots, x_{(n)}\}$, the median is given by the center element $x_{((n+1)/2)}$ if $n$ is odd and the mean of the two center elements $(x_{(n/2)} + x_{(n/2+1)})/2$ if $n$ is even.
- For `MatrixQ` data, the median is computed for each column vector. `Median` for a tensor gives columnwise medians at the first level.
- `Median` requires REAL numeric values; a complex element (including an evaluated `Complex[re, im]` such as `2 + I`, whose imaginary part is nonzero) is rejected with `Median::rectn`. The same real-only gate is shared by `Quartiles`, `Quantile`, `InterquartileRange`, `MeanDeviation`, `MedianDeviation` and `MovingMedian`; it does not yet see through a numeric head, so `Sqrt[2 + I]` still passes.

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

`Median[dist]` of a distribution is `Quantile[dist, 1/2]` when Quantile evaluates
for that distribution; otherwise it stays unevaluated.

## Mean

- `Mean[x]`, `Mean[5]` and `Mean[{}]` print `Mean::rectn` and stay unevaluated (audit Mean#3).
- `Mean[dist]` of a `NormalDistribution`, `UniformDistribution` or 1-D
  `SmoothKernelDistribution`: its closed form (see machine-learning.md, "CDF, Mean, ...").
Gives the mean estimate of the elements in data.
- `Mean[data]`

**Features**:
- `Protected`.
- Supports numerical and symbolic data.
- **Machine overflow / underflow** (v0.6595). The machine-double fast path
  declines a result outside the machine range to the symbolic path, which keeps
  it at machine precision in MPFR: `Mean[{1.*^308, 1.*^308}]` is `1.*^308` (it
  was `Overflow[]`); the packed / NDArray reduction does the same.
- For vectors, computes $(1/n) \sum x_i$.
- For matrices, computes means of elements in each column.

```mathematica
In[1]:= Mean[{1, 2, 3, 4}]
Out[1]= 5/2

In[2]:= Mean[{{a, u}, {b, v}, {c, w}}]
Out[2]= {1/3 (a + b + c), 1/3 (u + v + w)}
```

A ragged or mixed-depth list (`Mean[{{1, 2}, {3}}]`) is `Mean::rectn` and stays
unevaluated, as do `Variance` and `StandardDeviation`.

## RootMeanSquare
Gives the root mean square of values in `list`.
- `RootMeanSquare[list]`

**Features**:
- `Protected`.
- Gives the square root of the second sample moment.
- **Machine overflow / underflow** (v0.6595). The machine-double fast path
  declines a result outside the machine range to the symbolic path, which keeps
  it at machine precision in MPFR: `RootMeanSquare[{1.*^308, 1.*^308}]` is `1.*^308`
  (it was `Overflow[]`), and squares that underflow to 0 from nonzero data
  (`RootMeanSquare[{1.*^-200, 1.*^-200}]`, which was `0.`) decline too; the packed / NDArray reduction does the same.
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

- `Variance[dist]` of a `NormalDistribution`, `UniformDistribution` or 1-D
  `SmoothKernelDistribution`: its closed form (see machine-learning.md, "CDF, Mean, ...").
Gives the unbiased variance estimate of the elements in data.
- `Variance[data]`

**Features**:
- `Protected`.
- For vectors, computes $(1/(n-1)) \sum (x_i - \hat{\mu}) \overline{(x_i - \hat{\mu})}$.
- Machine overflow in the squares declines to the symbolic path (v0.6595), as
  for `Mean`; `Moment`, `CentralMoment` and `Covariance` follow the same rule.
- For matrices, computes variances of elements in each column.

```mathematica
In[1]:= Variance[{1, 2, 3}]
Out[1]= 1

In[2]:= Variance[{{5.2, 7}, {5.3, 8}, {5.4, 9}}]
Out[2]= {0.01, 1}
```

Fewer than two elements (or rows) is `Variance::shlen`. Distinct machine data whose
squared deviations underflow is no zero variance: it takes the arbitrary-precision path.

## Moment
Gives the raw (power) moment of data.
- `Moment[data, r]`: the order-`r` raw moment $\mu_r = (1/n) \sum_i x_i^r$ (the sum of `r`-th powers, divided by `n`).
- `Moment[data, {r_1, ..., r_m}]`: the multivariate mixed raw moment $(1/n) \sum_i \prod_j x_{ij}^{r_j}$.

**Features**:
- `NHoldAll`, `Protected`.
- The raw moment is `CentralMoment` without the mean subtraction; `Moment[data, 1]` is `Mean[data]`, and `Moment[data, 0]` is `1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[Moment[#, r]&, x, 1]`); because there is no mean to subtract, `Mean[data^r]` threads correctly at every rank.
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

```mathematica
In[1]:= Moment[{1, 2, 3, 4}, 2]
Out[1]= 15/2

In[2]:= Moment[{1., 2., 3., 4.}, 2]
Out[2]= 7.5

In[3]:= Moment[{Pi, E, 2}, 1]
Out[3]= 1/3 (2 + E + Pi)

In[4]:= Moment[{{1, 2}, {3, 4}, {5, 6}}, 3]
Out[4]= {51, 96}

In[5]:= Simplify[Moment[{{a, b}, {c, d}}, {1, 2}]]
Out[5]= 1/2 (a b^2 + c d^2)
```

## CentralMoment
- Ragged data (`{{1, 2}, {3}}`) is one `CentralMoment::rectn`, unevaluated (audit CentralMoment#2).
Gives the central moment (moment about the mean) of data.
- `CentralMoment[data, r]`: the order-`r` central moment $\tilde{\mu}_r = (1/n) \sum_i (x_i - \hat{\mu}_1)^r$, where $\hat{\mu}_1$ is `Mean[data]`.
- `CentralMoment[data, {r_1, ..., r_m}]`: the multivariate mixed central moment $(1/n) \sum_i \prod_j (x_{ij} - \hat{\mu}_{1,j})^{r_j}$.

**Features**:
- `Protected`.
- A central moment is `Variance` without the $n/(n-1)$ bias correction: it divides by `n` (not `n-1`), raises to the power `r` (not a square), and needs only `n >= 1`.
- For a matrix or array the moment is taken columnwise over the first axis (equivalent to `ArrayReduce[CentralMoment[#, r]&, x, 1]`).
- `CentralMoment[data, 0]` is `1` (the mean of the empty products `(x - mu)^0`, with `0^0 = 1` as in any sum of powers), exact for exact data and `1.` for approximate data (audit CentralMoment#1; the exact path used to raise `Power::indet`).
- Exact input yields exact output; approximate input yields approximate output; symbolic data is handled symbolically.
- Fast path on `NDArray`/packed real buffers (`ndred_central_moment`); an integer buffer degrades to the exact `Rational` `List` result, like `Variance`.
- Lowerable inside `Compile[]` for a real vector and integer order (participates in auto-compilation).

```mathematica
In[1]:= CentralMoment[{1, 2, 3, 4}, 4]
Out[1]= 41/16

In[2]:= CentralMoment[{1., 2., 3., 4.}, 2]
Out[2]= 1.25

In[3]:= CentralMoment[{{1, 2}, {3, 4}, {5, 6}}, 2]
Out[3]= {8/3, 8/3}

In[4]:= Simplify[CentralMoment[{{a, b}, {c, d}}, {2, 2}]]
Out[4]= 1/16 (a - c)^2 (b - d)^2
```

## Skewness
Gives the coefficient of skewness (a measure of asymmetry) for the elements in data.
- `Skewness[data]`

**Features**:
- `Protected`.
- Equivalent to `CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2)`.
- For a matrix, gives the columnwise skewnesses.
- Handles numerical and symbolic data; exact input gives exact output (a radical in general).
- Fast path on `NDArray`/packed real buffers (`ndred_skewness`) and lowerable inside `Compile[]`; an integer buffer degrades to the exact `List` result.

```mathematica
In[1]:= Skewness[{1, 2, 3, 10}]
Out[1]= 18/25 Sqrt[2]

In[2]:= Skewness[{1., 2., 3., 4., 5.}]
Out[2]= 0.
```

A machine vector takes the buffer kernel, so a short list answers what its packed
form does; when the moments leave the double range (data of size `1.*^-200` or
`1.*^120`) the deviations are rescaled by their largest first, the ratio being
scale-free. `Kurtosis` likewise.

## Kurtosis
Gives the coefficient of kurtosis (peak/tail versus flank concentration) for the elements in data.
- `Kurtosis[data]`

**Features**:
- `Protected`.
- Equivalent to `CentralMoment[data, 4] / CentralMoment[data, 2]^2` (Pearson kurtosis, not the excess form).
- For a matrix, gives the columnwise kurtoses.
- Handles numerical and symbolic data; exact input gives exact output.
- Fast path on `NDArray`/packed real buffers (`ndred_kurtosis`) and lowerable inside `Compile[]`; an integer buffer degrades to the exact `List` result.

```mathematica
In[1]:= Kurtosis[{1, 2, 3, 4, 5}]
Out[1]= 17/10

In[2]:= Kurtosis[{1, 2, 4, 8}]
Out[2]= 25141/13225
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

Too short (`{1}`, `{{1, 2}}`) is `StandardDeviation::shlen`; a non-list argument or
a call whose `Variance` stays unevaluated stays unevaluated (no `Sqrt[Variance[x]]`).
Machine data uses the buffer kernel's reduction, rescaled when the squares
under- or overflow (`StandardDeviation[{1.*^-200, 2.*^-200, 3.*^-200}]` is `1.*^-200`).

## Covariance
Gives the covariance between vectors, or the cross/auto-covariance matrix of matrices.
- `Covariance[v, w]`: covariance between the vectors `v` and `w` (a scalar).
- `Covariance[a, b]`: `p×q` cross-covariance matrix between the columns of the `n×p` and `n×q` matrices `a` and `b`.
- `Covariance[a]`: `p×p` auto-covariance matrix of the columns of `a`, i.e. `Covariance[a, a]`.

**Features**:
- `Protected`.
- For vectors, the unbiased estimate $\hat{\sigma}_{vw} = \frac{1}{n-1}\sum_i (v_i - \hat{\mu}_v)\overline{(w_i - \hat{\mu}_w)}$; the conjugate is on the **second** argument, so exact / complex / symbolic inputs yield exact / complex / symbolic output.
- For matrices, element $(i,j)$ is the covariance of column $i$ of `a` with column $j$ of `b`; `Covariance[a]` is symmetric.
- NDArray / packed real data uses a threaded centered inner product (vectors) or a BLAS gram (matrices); an integer sample degrades to the exact `List` path. Lowered inside `Compile[]`.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Covariance[]` reports `Covariance::argb`. Two vectors or matrices with different numbers of observations (rows) report `Covariance::vctmat` (`Correlation::vctmat`); the matrix form used to return a zero matrix (audit Covariance#1).

```mathematica
In[1]:= Covariance[{1, 3/2}, {2, 11}]
Out[1]= 9/4

In[2]:= Covariance[{2 + I, 3 - 2 I, 5 + 4 I}, {I, 1 + 2 I, 10 - 5 I}]
Out[2]= -7/3 + (56 I)/3

In[3]:= Covariance[{{1, 2}, {3, 4}, {5, 7}}]
Out[3]= {{4, 5}, {5, 19/3}}
```

## Correlation
Gives the correlation between vectors, or the cross/auto-correlation matrix of matrices.
- `Correlation[v, w]`: correlation between the vectors `v` and `w` (a scalar).
- `Correlation[a, b]`: `p×q` cross-correlation matrix between the columns of `a` and `b`.
- `Correlation[a]`: `p×p` auto-correlation matrix of the columns of `a`.

**Features**:
- `Protected`.
- A normalized covariance, $\rho_{vw} = \sigma_{vw} / (\sigma_v\,\sigma_w)$ with $\sigma_{vw} = \mathtt{Covariance}[v,w]$ and $\sigma_v = \mathtt{StandardDeviation}[v]$; $-1 \le \rho_{vw} \le 1$ for real data.
- The auto-correlation matrix `Correlation[a]` is symmetric with a unit diagonal (exact `1` for exact/symbolic data, `1.` for real data).
- Shares `Covariance`'s NDArray / packed / `Compile[]` fast paths.
- Stays unevaluated for a single vector, mismatched shapes, or fewer than two observations. `Correlation[]` reports `Correlation::argb`.

```mathematica
In[1]:= Correlation[{5, 3/4, 1}, {2, 1/2, 1}]
Out[1]= 2 Sqrt[3/13]

In[2]:= Correlation[{1.5, 3, 5, 10}, {2, 1.25, 15, 8}]
Out[2]= 0.475976

In[3]:= Correlation[{{a, b}, {c, d}}][[1, 1]]
Out[3]= 1
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


## Quantile

- `Quantile[dist]`, `Quantile[dist, q]` of a `NormalDistribution`, `UniformDistribution` or 1-D
  `SmoothKernelDistribution`: its closed form (see machine-learning.md, "CDF, Mean, ...").
Gives the q-th quantile estimate of the elements in data.
- `Quantile[data, q]`: the q-th quantile of `data` (left-continuous: for sorted data $\{x_{(1)},\dots,x_{(n)}\}$, $x_{(\lceil nq \rceil)}$, edge-clamped).
- `Quantile[data, {q1, q2, ...}]`: a list of quantiles, one per $q_i$.
- `Quantile[data, q, {{a, b}, {c, d}}]`: the general parameterized quantile: with $h = a + (n+b)q$ and $w = c + d\,(h - \lfloor h\rfloor)$, the result is the convex combination $(1-w)\,x_{(\lfloor h\rfloor)} + w\,x_{(\lceil h\rceil)}$, edge-clamped; at integer $h$ this is $x_{(h)}$ for any $c, d$.

**Features**:
- `Protected`.
- Default parameters are `{{0, 0}, {1, 0}}` (the reference's Type-1 / left-continuous inverse CDF), so `Quantile[{1, 2, 3, 4}, 1/2]` is `2` while `Median[{1, 2, 3, 4}]` is `5/2` — the two deliberately differ on even-length data.
- `Quartiles[data]` is `Quantile[data, {1/4, 1/2, 3/4}, {{1/2, 0}, {0, 1}}]`. A second argument that is not `{{a, b}, {c, d}}` prints `Quartiles::parm` and the call stays unevaluated (audit Quartiles#4).
- Exact input gives exact output. The data is ordered by **numeric value** (certified: exact, then MPFR at increasing precision when two values are close), not by the canonical `Sort` order, which put every explicit number before `Sqrt[2]`, `E` or `Pi`: `Quantile[{2, Sqrt[3]}, 0]` is `Sqrt[3]` (audit Quantile#1). `Median`, `Quartiles`, `InterquartileRange` and `MedianDeviation` share the same order (audit Median#1/#2, InterquartileRange#1, MedianDeviation#1).
- Ragged or mixed-depth data (`{{1, 2}, {3}}`, `{{1, 2}, 3}`) gives `Quantile::rectn` and stays unevaluated; it used to recurse through an unevaluated `Transpose` to `$RecursionLimit` (audit Quantile#2). `Quartiles`, `InterquartileRange`, `Median` and `MedianDeviation` likewise give their own `::rectn` (audit Quartiles#1, InterquartileRange#2).
- For $w \in [0, 1]$ the interpolation is evaluated as the convex combination $(1-w)x_{(j)} + w\,x_{(j+1)}$, not as $x_{(j)} + w\,(x_{(j+1)} - x_{(j)})$. The two are equal in exact arithmetic, but the difference form overflows on `Real` data whose neighbours straddle zero near the machine range: `Quantile[{-1.0*10^308, 1.0*10^308}, 1/2, {{1/2, 0}, {0, 1}}]` is `0.0`, not `Infinity`. Outside $[0, 1]$ — which `{{a,b},{c,d}}` permits, though no standard quantile type uses it — the difference form is used instead, because there the convex form is the one that can produce `NaN`. At $w = 0$ or $w = 1$ the element is selected outright rather than computed.
- For `MatrixQ` data the quantile is computed per column.
- `Quantile` requires REAL numeric data and numeric $q \in [0, 1]$ (`Quantile::q100` otherwise); symbolic arguments stay unevaluated. A `Complex[re, im]` element is rejected with `Quantile::rectn` when `im` is nonzero — including an already-evaluated complex such as `2 + I` (`Complex[2, 1]`), which carries no literal `I` to search for — and accepted when `im` is zero, since `Complex[x, 0]` at MPFR precision is a real number. Not yet caught: a complex value nested under a numeric head, such as `Sqrt[2 + I]`.
- An `NDArray` argument is materialised to the exact `List` path (no buffer fast path yet).

```mathematica
In[1]:= Quantile[{1, 2, 3, 4}, 1/2]
Out[1]= 2

In[2]:= Quantile[{1, 2, 3, 4}, {1/4, 3/4}]
Out[2]= {1, 3}

In[3]:= Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 1}}]
Out[3]= 5/2

In[4]:= Quantile[{3, 1, 4, 2}, 1/2]
Out[4]= 2

In[5]:= Quantile[{{1, 2}, {3, 4}}, 1/2]
Out[5]= {1, 2}

In[6]:= Quantile[{1, 2, 3}, 2]
Quantile::q100: The quantile q is expected to be a number between 0 and 1 in Quantile[{1, 2, 3}, 2].
Out[6]= Quantile[{1, 2, 3}, 2]
```

`Quartiles[list, {{a, b}, {c, d}}]` needs numeric parameters; a symbolic one leaves
the call unevaluated.

## InterquartileRange
- `InterquartileRange[data, {{a, b}, {c, d}}]` is `q3 - q1` under those quantile parameters (audit InterquartileRange#4).
Gives the difference between the upper and lower quartiles of data.
- `InterquartileRange[data]`: $\hat{q}_{3/4} - \hat{q}_{1/4}$, using the `Quartiles` parameterization `{{1/2, 0}, {0, 1}}`.

**Features**:
- `Protected`.
- A robust scale estimator: insensitive to outliers beyond the quartiles.
- For `MatrixQ` data the IQR is computed per column.
- Requires numeric data (`InterquartileRange::rectn` otherwise).

```mathematica
In[1]:= InterquartileRange[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[1]= 4

In[2]:= InterquartileRange[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]
Out[2]= {6, 6, 6}
```

## MeanDeviation
Gives the mean absolute deviation from the mean of the elements in data.
- `MeanDeviation[data]`: `Mean[Abs[data - Mean[data]]]`.

**Features**:
- `Protected`.
- Exact input gives exact output: `MeanDeviation[{1, 2, 3, 4}]` is `1`.
- For `MatrixQ` data the deviation is computed per column.
- Requires numeric data (`MeanDeviation::rectn` otherwise).

```mathematica
In[1]:= MeanDeviation[{1, 2, 3, 4}]
Out[1]= 1

In[2]:= MeanDeviation[{1/2, 3/2}]
Out[2]= 1/2
```

## MedianDeviation
Gives the median absolute deviation from the median of the elements in data.
- `MedianDeviation[data]`: `Median[Abs[data - Median[data]]]` (the MAD, a robust scale estimator).

**Features**:
- `Protected`.
- Exact input gives exact output: `MedianDeviation[{1, 2, 3, 4}]` is `1`.
- For `MatrixQ` data the deviation is computed per column.
- Requires numeric data (`MedianDeviation::rectn` otherwise).

```mathematica
In[1]:= MedianDeviation[{1, 2, 3, 4}]
Out[1]= 1

In[2]:= MedianDeviation[{1, 2, 3, 10}]
Out[2]= 1
```
