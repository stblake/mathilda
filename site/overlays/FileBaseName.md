### Worked examples

```mathematica
In[1]:= FileBaseName["/tmp/data/report.txt"]
Out[1]= "report"

In[2]:= FileBaseName["archive.tar.gz"]
Out[2]= "archive.tar"
```

### Notes

`FileBaseName` strips any directory components and the final extension, so a name with multiple dots keeps all but the last.
