# Random Number Generation

## RandomInteger
Gives pseudorandom integers.
- `RandomInteger[{imin, imax}]`: gives a pseudorandom integer in the range {imin, ..., imax}.
- `RandomInteger[imax]`: gives a pseudorandom integer in the range {0, ..., imax}.
- `RandomInteger[]`: pseudorandomly gives 0 or 1.
- `RandomInteger[range, n]`: gives a list of n pseudorandom integers.
- `RandomInteger[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom integers.

**Features**:
- `Protected`.
- RandomInteger[{imin, imax}] chooses integers in the range {imin, ..., imax} with equal probability.
- RandomInteger[] gives 0 or 1 with probability 1/2.
- RandomInteger gives a different sequence of pseudorandom integers whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Returns bignums when the range exceeds 64-bit integer limits.

```mathematica
In[1]:= SeedRandom[42]; RandomInteger[]
Out[1]= 1

In[2]:= SeedRandom[42]; RandomInteger[10]
Out[2]= 6

In[3]:= SeedRandom[42]; RandomInteger[{1, 6}]
Out[3]= 4

In[4]:= SeedRandom[42]; RandomInteger[{0, 9}, 5]
Out[4]= {6, 9, 4, 1, 3}

In[5]:= SeedRandom[42]; Dimensions[RandomInteger[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomInteger[{-10, -5}]
Out[6]= -6

In[7]:= IntegerQ[RandomInteger[10^20]]
Out[7]= True
```

## RandomReal
Gives pseudorandom real numbers.
- `RandomReal[]`: gives a pseudorandom real number in the range 0 to 1.
- `RandomReal[{xmin, xmax}]`: gives a pseudorandom real number in the range xmin to xmax.
- `RandomReal[xmax]`: gives a pseudorandom real number in the range 0 to xmax.
- `RandomReal[range, n]`: gives a list of n pseudorandom reals.
- `RandomReal[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom reals.

**Features**:
- `Protected`.
- RandomReal[{xmin, xmax}] chooses reals with a uniform probability distribution in the range xmin to xmax.
- RandomReal gives a different sequence of pseudorandom reals whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness for full double-precision mantissa coverage.
- Accepts integer, real, rational, and bigint range arguments.

```mathematica
In[1]:= SeedRandom[42]; RandomReal[]
Out[1]= 0.376082

In[2]:= SeedRandom[42]; RandomReal[10]
Out[2]= 3.76082

In[3]:= SeedRandom[42]; RandomReal[{-1, 1}]
Out[3]= -0.247836

In[4]:= SeedRandom[42]; Length[RandomReal[{0, 1}, 5]]
Out[4]= 5

In[5]:= SeedRandom[42]; Dimensions[RandomReal[{0, 1}, {3, 4}]]
Out[5]= {3, 4}

In[6]:= SeedRandom[42]; RandomReal[{0, 1}, 0]
Out[6]= {}

In[7]:= RandomReal[x]
Out[7]= RandomReal[x]
```

## SeedRandom
Resets the pseudorandom generator.
- `SeedRandom[n]`: seeds the generator with integer n.
- `SeedRandom[]`: reseeds from system entropy.

**Features**:
- `Protected`.
- After `SeedRandom[n]`, the sequence of pseudorandom numbers generated will be the same each time.
- Accepts bignums as seeds.

```mathematica
In[1]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[1]= {1, 1, 0}

In[2]:= SeedRandom[42]; {RandomInteger[], RandomInteger[], RandomInteger[]}
Out[2]= {1, 1, 0}
```

## RandomComplex
Gives pseudorandom complex numbers.
- `RandomComplex[]`: gives a pseudorandom complex number with real and imaginary parts in the range 0 to 1.
- `RandomComplex[{zmin, zmax}]`: gives a pseudorandom complex number in the rectangle with corners given by the complex numbers zmin and zmax.
- `RandomComplex[zmax]`: gives a pseudorandom complex number in the rectangle whose corners are the origin and zmax.
- `RandomComplex[range, n]`: gives a list of n pseudorandom complex numbers.
- `RandomComplex[range, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom complex numbers.

**Features**:
- `Protected`.
- `RandomComplex[{zmin, zmax}]` chooses complex numbers uniformly in the rectangle with corners at `zmin` and `zmax`.
- RandomComplex gives a different sequence of pseudorandom complex numbers whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Uses 53 bits of randomness per component for full double-precision mantissa coverage.
- Accepts integer, real, rational, and complex range arguments. When the range has no imaginary component, the result simplifies to a real.

```mathematica
In[1]:= SeedRandom[42]; Head[RandomComplex[]]
Out[1]= Complex

In[2]:= SeedRandom[42]; RandomComplex[2 + 3 I]
Out[2]= 0.752164 + 1.30654 I

In[3]:= SeedRandom[42]; Length[RandomComplex[{0, 1 + I}, 5]]
Out[3]= 5

In[4]:= SeedRandom[42]; Dimensions[RandomComplex[{0, 1 + I}, {3, 4}]]
Out[4]= {3, 4}

In[5]:= RandomComplex[x]
Out[5]= RandomComplex[x]
```

## RandomChoice
Gives pseudorandom choices from a list of elements.
- `RandomChoice[{e1, e2, ...}]`: gives a pseudorandom choice of one of the ei.
- `RandomChoice[list, n]`: gives a list of n pseudorandom choices.
- `RandomChoice[list, {n1, n2, ...}]`: gives an n1 x n2 x ... array of pseudorandom choices.
- `RandomChoice[{w1, w2, ...} -> {e1, e2, ...}]`: gives a pseudorandom choice weighted by the wi.
- `RandomChoice[wlist -> elist, n]`: gives a list of n weighted choices.
- `RandomChoice[wlist -> elist, {n1, n2, ...}]`: gives an n1 x n2 x ... array of weighted choices.

**Features**:
- `Protected`.
- `RandomChoice[{e1, e2, ...}]` chooses with equal probability between all of the ei.
- RandomChoice gives a different sequence of pseudorandom choices whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Weighted selection uses cumulative weight binary search for efficient O(log n) per choice.
- Weights must be non-negative real numbers with a positive total.

```mathematica
In[1]:= SeedRandom[42]; RandomChoice[{a, b, c, d, e}]
Out[1]= c

In[2]:= SeedRandom[42]; RandomChoice[{a, b, c}, 5]
Out[2]= {b, c, b, a, a}

In[3]:= SeedRandom[42]; Dimensions[RandomChoice[{1, 2, 3}, {3, 4}]]
Out[3]= {3, 4}

In[4]:= RandomChoice[{1, 0, 0} -> {a, b, c}]
Out[4]= a

In[5]:= RandomChoice[{1, 0} -> {x, y}, 5]
Out[5]= {x, x, x, x, x}

In[6]:= RandomChoice[x]
Out[6]= RandomChoice[x]
```

## RandomSample
Gives a pseudorandom sample of elements without replacement.
- `RandomSample[{e1, e2, ...}, n]`: gives a pseudorandom sample of n of the ei.
- `RandomSample[{w1, w2, ...} -> {e1, e2, ...}, n]`: gives a pseudorandom sample of n of the ei chosen using weights wi.
- `RandomSample[{e1, e2, ...}]`: gives a pseudorandom permutation of the ei.
- `RandomSample[list, UpTo[n]]`: gives a sample of n of the ei, or as many as are available.

**Features**:
- `Protected`.
- `RandomSample[{e1, e2, ...}, n]` never samples any of the ei more than once.
- `RandomSample[{e1, e2, ...}, n]` samples each of the ei with equal probability.
- `RandomSample[{e1, e2, ...}, UpTo[n]]` gives a sample of n of the ei, or as many as are available.
- RandomSample gives a different sequence of pseudorandom choices whenever you run PicoCAS. You can start with a particular seed using SeedRandom.
- Requesting n greater than the list length (without UpTo) returns unevaluated.
- Uses the Fisher-Yates shuffle for uniform sampling without replacement.
- Weighted sampling removes selected elements and renormalizes weights.

```mathematica
In[1]:= SeedRandom[42]; RandomSample[{a, b, c, d, e}, 3]
Out[1]= {d, e, b}

In[2]:= Sort[RandomSample[{1, 2, 3, 4, 5}, 5]]
Out[2]= {1, 2, 3, 4, 5}

In[3]:= Length[RandomSample[{a, b, c, d, e}]]
Out[3]= 5

In[4]:= RandomSample[{a, b, c}, 0]
Out[4]= {}

In[5]:= Length[RandomSample[{a, b, c, d, e}, UpTo[10]]]
Out[5]= 5

In[6]:= RandomSample[{1, 0, 0} -> {a, b, c}, 1]
Out[6]= {a}

In[7]:= Sort[RandomSample[{1, 1, 0} -> {a, b, c}, 2]]
Out[7]= {a, b}

In[8]:= Sort[RandomSample[{1, 2, 3} -> {a, b, c}]]
Out[8]= {a, b, c}

In[9]:= RandomSample[{a, b}, 5]
Out[9]= RandomSample[{a, b}, 5]

In[10]:= RandomSample[x]
Out[10]= RandomSample[x]
```

