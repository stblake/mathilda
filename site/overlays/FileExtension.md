### Worked examples

```mathematica
In[1]:= FileExtension["report.txt"]
Out[1]= "txt"

In[2]:= FileExtension["/tmp/data/report.txt"]
Out[2]= "txt"
```

### Notes

`FileExtension` returns the extension (without the dot), taken after the final dot in the file name.
