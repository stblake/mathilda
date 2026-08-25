# $RaylibVerbose

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`$RaylibVerbose`**

controls whether the Raylib graphics backend prints its diagnostic log (window and OpenGL initialisation lines) to the terminal. False by default, so opening a Show/Plot window or exporting an image stays quiet; set it to True to see Raylib's full trace log.

<details>
<summary>Notes</summary>

Affects only the terminal chatter from the renderer, never the graphics produced. Has no effect in a build without graphics (USE\_GRAPHICS=0). Only True or False is accepted.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
