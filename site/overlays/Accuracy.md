### Worked examples

```mathematica
In[1]:= Accuracy[1.5]
Out[1]= 15.7785

In[2]:= Accuracy[5]
Out[2]= Infinity

In[3]:= Accuracy[N[Pi, 30]]
Out[3]= 29.6058
```

### Notes

`Accuracy` gives the number of digits to the right of the decimal point, i.e. the digit count relative to the absolute (not relative) magnitude. Exact numbers have `Infinity` accuracy; for a fixed precision, accuracy decreases as the magnitude grows.
