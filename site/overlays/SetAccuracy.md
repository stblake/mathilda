### Worked examples

```mathematica
In[1]:= SetAccuracy[1.5, 30]
Out[1]= 1.5

In[2]:= Accuracy[SetAccuracy[1.5, 30]]
Out[2]= 30.2279
```

Promoting an exact constant to a fixed accuracy yields its high-accuracy decimal
expansion — here 30 digits past the point for `Pi`:

```mathematica
In[1]:= SetAccuracy[Pi, 30]
Out[1]= 3.141592653589793238462643383279
```

Forcing accuracy *beyond* the genuine precision of a machine number exposes the
binary round-off in its tail (the digits past `123.456` are noise):

```mathematica
In[1]:= SetAccuracy[123.456, 20]
Out[1]= 123.45600000000000306954
```

### Notes

`SetAccuracy[x, n]` returns a value equal to `x` with `n` digits of accuracy (digits past the decimal point); use `Accuracy` to confirm, since the printed form often looks unchanged. It is the absolute-magnitude counterpart to `SetPrecision`.
