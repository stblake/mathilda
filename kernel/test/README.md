# Kernel tests

Two layers, both validated on macOS (xeus 6.0.6, AppleClang 17, gcc-16):

## 1. `ffi_cell_test.c` — the C cell-evaluation layer (no xeus needed)

Exercises `mathilda_ffi_eval_cell`, `mathilda_ffi_is_complete`, and
`mathilda_ffi_complete` directly against `libmathilda.a`, asserting the same
event semantics as the sidecar's NDJSON pipe mode (`make check-pipe-protocol`).

```bash
# from the repo root
make libmathilda.a
LDF=$(make --dry-run --always-make Mathilda 2>/dev/null \
      | grep -E '\-o Mathilda' | head -1 | grep -oE '([^ ]+\.o )(.*)$' | sed -E 's/^.*\.o //')
gcc-16 -std=c99 -g -I./src kernel/test/ffi_cell_test.c libmathilda.a $LDF -o /tmp/ffi_cell_test
MATHILDA_NO_WINDOW=1 /tmp/ffi_cell_test ./src
```

## 2. `kernel_test.py` — the full Jupyter kernel (needs the xeus env)

Launches the installed `xmathilda` kernel through `jupyter_client` and checks
`execute_request` (LaTeX result, multi-statement cell, `Print` stream,
`Power::infy` message, syntax-error status, Plotly plot), `complete_request`,
`inspect_request`, `is_complete_request`, `kernel_info_request`, and `%`
history.

```bash
# with the xeus env active and the kernelspec installed (see ../README.md)
JUPYTER_PATH="$CONDA_PREFIX/share/jupyter" MATHILDA_NO_WINDOW=1 \
    python kernel/test/kernel_test.py
```
