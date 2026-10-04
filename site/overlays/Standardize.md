### Worked examples

```mathematica
In[1]:= Standardize[{1., 2., 3., 4., 5.}]  (* a flat list is n observations of one variable *)
```

```mathematica
In[1]:= Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}]  (* each column is standardised independently *)
```

```mathematica
In[1]:= Standardize[{{1., 5.}, {2., 5.}, {3., 5.}}]  (* a constant column becomes exactly 0, not Indeterminate *)
```

### Notes

Columns are variables and rows are observations: a flat list is `n` observations of a
*single* variable, not one observation of `n`.

The divisor is `n − 1`, the sample standard deviation, so `Standardize[x]` agrees to the
last digit with `(x − Mean[x])/StandardDeviation[x]` written out by hand. A constant
column carries no information, so it is shifted to exactly `0` rather than divided by its
zero standard deviation — which would propagate `Indeterminate` through every reduction
over the row.
