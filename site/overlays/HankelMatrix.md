### Worked examples

```mathematica
In[1]:= HankelMatrix[3]
Out[1]= {{1, 2, 3}, {2, 3, 0}, {3, 0, 0}}
```

Give an explicit first column; the matrix is constant along every
antidiagonal and zero-padded below the secondary diagonal:

```mathematica
In[1]:= HankelMatrix[{a, b, c}]
Out[1]= {{a, b, c}, {b, c, 0}, {c, 0, 0}}
```

A first column plus a last row builds a rectangular catalecticant — here the
shared corner is taken from the column:

```mathematica
In[1]:= HankelMatrix[{1, 2, 3, 4}, {4, 5, 6}]
Out[1]= {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {4, 5, 6}}
```

Hankel determinants encode sequence properties. The catalecticant of the
Fibonacci numbers is a perfect power of 2:

```mathematica
In[1]:= Det[HankelMatrix[{1, 1, 2, 3, 5, 8}]]
Out[1]= -262144
```

### Notes

A Hankel matrix is constant along its antidiagonals. With a single integer `n`, the first row and column are the integers `1..n` and the lower-right triangle is zero-filled. Supplying a first column (and optionally a last row) lets you build the catalecticant of any sequence; `Det[HankelMatrix[...]]` then gives that sequence's Hankel determinant.
