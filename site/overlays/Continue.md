### Worked examples

```mathematica
In[1]:= r = 0; Do[If[EvenQ[i], Continue[]]; r += i, {i, 10}]; r  (* skip the even i and sum the odd ones *)
```

```mathematica
In[1]:= r = 0; For[i = 1, i <= 6, i++, If[i == 3, Continue[]]; r += i]; r  (* in For, Continue still runs the increment step *)
```

```mathematica
In[1]:= Continue[]  (* inert outside any loop *)
```

### Notes

`Continue[]` abandons the rest of the current loop body and moves to the next
iteration of the innermost `Do`, `For` or `While`. What "next iteration" means
differs by loop: `Do` advances its iterator and re-tests, `For` evaluates its
increment step and re-tests, and `While` re-evaluates its test.

Like `Break`, it is a head-detected marker and takes effect only at a loop
boundary. Outside any loop it prints `Continue::nofwd` and returns the inert
`Hold[Continue[]]`.
