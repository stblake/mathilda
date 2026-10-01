# xmathilda — a xeus-based Jupyter kernel for Mathilda

A native C++ Jupyter kernel that runs Mathilda in process. It embeds the
evaluator through Mathilda's C ABI (`src/ffi/mathilda_ffi.h`) and lets
[xeus](https://github.com/jupyter-xeus/xeus) handle the Jupyter wire protocol,
so there is **no Python dependency at run time** — only the Jupyter *front end*
(Lab, Notebook, `nbclient`, VS Code, …) is Python, and it is not this kernel's
concern.

This is an interop on-ramp to the standard Jupyter ecosystem; it complements,
rather than replaces, the Tauri desktop notebook in [`../frontend`](../frontend).

## What it supports

| Jupyter request | Mathilda behaviour |
|---|---|
| `execute_request` | A cell may hold several statements (split like a Mathematica input cell). `Print` → stdout stream; `Head::tag` warnings → stderr; a `;`- or `Null`-valued statement shows nothing; results carry `text/latex` (KaTeX) + `text/plain`. The last result is the `execute_result` (`Out[n]`), earlier ones are `display_data`. `%`, `%%`, `Out[n]` work across cells. |
| `Plot[…]` / `Graphics[…]` | `application/vnd.plotly.v1+json` (renders in JupyterLab/nbviewer with the Plotly renderer). |
| `complete_request` | Tab completion over defined symbol names (prefix match). |
| `inspect_request` | `Shift-Tab` / `?name` docstring. |
| `is_complete_request` | Console continuation (balanced brackets / strings / comments). |

Follow-ups: `image/png` for `Image[…]` (currently `application/json` of raw
RGBA), a self-contained HTML+plotly.js fallback for front ends without the
Plotly renderer, and true mid-computation interrupt (needs an abort flag in the
evaluator loop).

## Build

The kernel needs the xeus toolchain (xeus ≥ 5, xeus-zmq ≥ 3, `nlohmann_json`).
The simplest source is conda-forge:

```bash
micromamba create -n mathilda-xeus -c conda-forge \
    xeus xeus-zmq cppzmq nlohmann_json xtl cmake jupyter_client
micromamba activate mathilda-xeus
```

Then, from the repository root:

```bash
cmake -S kernel -B kernel/build -DCMAKE_BUILD_TYPE=Release
cmake --build kernel/build -j
cmake --install kernel/build --prefix "$CONDA_PREFIX"   # installs the kernelspec
```

`libmathilda.a` is built automatically by the Mathilda makefile (with its own
GCC toolchain — it meets the C++ kernel only at the portable C ABI). For a
leaner, headless kernel, drop optional dependencies:

```bash
cmake -S kernel -B kernel/build \
    -DMATHILDA_MAKE_FLAGS="USE_GRAPHICS=0;USE_FFTW=0;USE_ECM=0"
```

Plotly output still works with `USE_GRAPHICS=0` (the `Graphics → Plotly JSON`
serializer is compiled regardless; only the on-screen/GL raster paths drop).

On Linux, if `libmathilda.a` was built against reference LAPACK, pass its link
libraries: `-DMATHILDA_EXTRA_LIBS="lapacke;lapack;blas;gfortran"`.

## Use

```bash
jupyter kernelspec list            # shows "mathilda"
jupyter console --kernel xmathilda # or open a notebook and pick the Mathilda kernel
```

```
In[1]:= Integrate[x^2, x]          (* typeset result *)
In[2]:= Table[Prime[k], {k, 5}]
In[3]:= Plot[Sin[x], {x, 0, 2 Pi}] (* Plotly figure *)
```
