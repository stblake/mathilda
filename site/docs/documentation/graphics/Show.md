# Show

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Show[graphics, opts...]`**

Displays graphics (a Graphics\[...\] object) in an interactive window and returns it, merging any given options into its option list.

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Options (1)

```mathematica
In[1]:= Show[Graphics[{Point[{0,0}]}], Axes -> True]
Out[1]= -Graphics-
```

## Implementation notes

- `Protected`.
- Declines to evaluate (stays unevaluated) if its first argument isn't a
  `Graphics[...]` expression, or any trailing argument isn't a `Rule`.
- When Raylib isn't compiled in, prints a one-line message instead of
  opening a window; the option merge still happens and the merged
  `Graphics[...]` is still returned.

**Attributes:** `Protected`.

## See also

[Plot](../../graphics/Plot/), [Frame](../../other-advanced/Frame/), [FrameTicks](../../other-advanced/FrameTicks/), [FrameStyle](../../other-advanced/FrameStyle/), [AspectRatio](../../other-advanced/AspectRatio/), [ImageSize](../../other-advanced/ImageSize/), [Tan](../../elementary-functions/Tan/), [Sec](../../elementary-functions/Sec/)

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
