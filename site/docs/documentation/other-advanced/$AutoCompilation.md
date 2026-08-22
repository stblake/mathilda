# $AutoCompilation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$AutoCompilation`**

controls whether Mathilda compiles numeric bodies to bytecode behind the scenes. True by default; set it to False to force every such body through the interpreter.

<details>
<summary>Notes</summary>

Covers both automatic mechanisms: the adapter that compiles a held body once for many sample points (Plot, Plot3D, Table, NIntegrate, NSum, FindRoot, the plot samplers) and the numeric-loop compiler for Do, For, While, Map, Nest, Fold and FixedPoint bodies. Compile\[\] and any CompiledFunction the user built explicitly are NOT affected -- those were asked for. A compiled body is contracted to give the interpreter's answer, so this changes speed and nothing else; it exists so the two paths can be compared. Reads back False in a session started with the environment variable MATHILDA\_NO\_AUTOCOMPILE set. Only True or False is accepted.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
