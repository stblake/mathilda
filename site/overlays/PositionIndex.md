### Worked examples

```mathematica
In[1]:= PositionIndex[{a, b, a, c, b, a}]  (* each value -> the positions it occurs at *)
```

### Notes

The result maps each distinct element to the sorted list of 1-based positions where
it appears, in first-appearance order of the values. It is the inverse view of a list
and runs in one O(n) hash pass. `PositionIndex[assoc]` instead maps each distinct
value to the list of keys that hold it.
