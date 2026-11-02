# Mathilda Canvas Library — Full Implementation Plan

## What We're Building

An **infinite canvas** where every "library" (`.lb` file) contains multiple named
notebooks positioned in 2D space. Pan, pinch-zoom, two-finger scroll. Glass-dark
aesthetic. Works on macOS, Linux, and Windows as a single Tauri binary.

---

## 1. File Format — `.lb` (Library)

JSON, human-readable, Git-diffable. No stored outputs (ephemeral by design).

```jsonc
{
  "version": "1",
  "type":    "mathilda-library",
  "title":   "Calculus Toolkit",
  "notebooks": [
    {
      "id":        "nb-1",
      "title":     "Derivatives",
      "x":         60,
      "y":         60,
      "width":     640,
      "collapsed": false,
      "rows": [
        { "cells": [{ "type": "code", "source": "D[Sin[x]^2, x]" }] },
        { "cells": [{ "type": "code", "source": "D[x^x, x]" }] },
        { "cells": [{ "type": "section", "source": "Chain Rule" }] },
        { "cells": [{ "type": "code", "source": "D[Sin[Exp[x]], x]" }] }
      ]
    },
    {
      "id":    "nb-2",
      "title": "Integrals",
      "x":     740,
      "y":     60,
      ...
    }
  ]
}
```

**File extension:** `.lb`  
**Default save location:** `~/Documents/Mathilda/`  
**Auto-save:** debounced 2 s after last edit

---

## 2. Example Libraries (ship with the app)

| File | Contents | Notebooks |
|---|---|---|
| `examples/calculus.lb` | Derivatives, integrals, limits, series | 4 |
| `examples/number-theory.lb` | Primes, factoring, GCD, EulerPhi | 3 |
| `examples/linear-algebra.lb` | Matrices, eigenvalues, SVD | 3 |
| `examples/special-functions.lb` | Gamma, Zeta, Bessel, elliptic | 3 |
| `examples/plots.lb` | 2D plots, parametric, polar | 3 |

All examples use only functions already implemented in Mathilda.

---

## 3. Canvas UI

### Pan & Zoom (gestures)
```
macOS trackpad:
  Two-finger drag      → pan   (wheel event, ctrlKey=false)
  Pinch in/out         → zoom  (wheel event, ctrlKey=true, macOS sends fractional deltaY)

Mouse:
  Scroll wheel         → zoom  (centered on cursor)
  Click+drag on bg     → pan

Keyboard:
  Cmd+0                → fit all notebooks in viewport
  Cmd+= / Cmd+-        → zoom in / out  
  Cmd+N                → new notebook at viewport centre
```

### Spring Physics
Pan and zoom are spring-smoothed via `requestAnimationFrame`:
```
targetPanX, targetPanY, targetZoom  ← set instantly on input
panX, panY, zoom                    ← lerp toward target each frame (k=0.15)
```
Gives the fluid Figma/Motion feel without a physics library.

### Auto-collapse
- `zoom < 0.40` → all cards auto-collapse to compact preview
- `zoom >= 0.40` → cards expand back
- Individual collapse/expand always overrides

### Dot Grid
CSS `background-image: radial-gradient(circle, rgba(255,255,255,0.08) 1px, transparent 1px)`  
Grid spacing scales with zoom: `background-size: ${40*zoom}px ${40*zoom}px`  
Offset with pan: `background-position: ${panX}px ${panY}px`

### Card Design (glass-dark)
```css
background:       rgba(12, 15, 28, 0.85);
backdrop-filter:  blur(20px) saturate(140%);
border:           1px solid rgba(255,255,255,0.08);
box-shadow:       0 0 0 1px rgba(137,180,250,0.10),
                  0 24px 64px rgba(0,0,0,0.65),
                  0 1px 0 rgba(255,255,255,0.05) inset;
border-radius:    12px;
```

Title bar has a subtle gradient accent tinted to the notebook's colour.
Each notebook gets a unique hue from a 6-colour palette.

---

## 4. Implementation Tasks (Parallel)

### Track A — Canvas UI
Files: `Canvas.svelte`, `NotebookCard.svelte`, `CellShell.svelte` (store prop refactor)

- [ ] `Canvas.svelte`: pan/zoom/pinch wheel handler, dot grid, spring loop, Cmd+N/0
- [ ] `NotebookCard.svelte`: glass card, pointer-drag, collapse anim, inline title edit, colour accent, per-card store
- [ ] `CellShell.svelte`: replace singleton `notebook` import with `store` prop; fix all calls
- [ ] `App.svelte`: thin shell rendering `<Canvas />`

### Track B — Library Format + Save/Load
Files: `src-tauri/src/commands.rs`, `canvas.ts`, `ipc.ts`

- [ ] `commands.rs`: `save_library(path, json)`, `load_library(path) → json`
- [ ] `canvas.ts`: `saveLibrary()`, `loadLibrary()` serialise/deserialise `CanvasNotebook[]`
- [ ] `ipc.ts`: typed wrappers for new Tauri commands
- [ ] `App.svelte`: Cmd+S / Cmd+O wired to library I/O, title bar shows filename
- [ ] Auto-save store (debounced 2s)

### Track C — Example Libraries
Files: `examples/*.lb`

- [ ] `examples/calculus.lb` — 4 notebooks, 16+ cells
- [ ] `examples/number-theory.lb` — 3 notebooks, 12+ cells
- [ ] `examples/linear-algebra.lb` — 3 notebooks, 12+ cells
- [ ] `examples/special-functions.lb` — 3 notebooks, 10+ cells
- [ ] `examples/plots.lb` — 3 notebooks, 9+ cells

### Track D — Cross-Platform Build
Files: `src/repl.c`, `makefile`, `build-sidecar.sh`, `docs/building.md`

- [ ] **macOS**: Already working (aarch64, x86_64)
- [ ] **Linux x86_64**: `apt-get` deps, `pkg-config`, makefile autodetect
- [ ] **Windows x86_64**: `_isatty`/`_fileno` guards, no-readline fallback, MinGW/MSYS2 GMP
- [ ] `docs/building.md`: step-by-step for each platform
- [ ] Optional: GitHub Actions CI matrix (mac/linux/windows)

---

## 5. Cross-Platform C Binary Fixes

### Windows (repl.c)
```c
/* Portable isatty + fileno */
#ifdef _WIN32
  #include <io.h>
  #ifndef isatty
    #define isatty _isatty
  #endif
  #ifndef fileno
    #define fileno _fileno
  #endif
#else
  #include <unistd.h>
#endif
```

### Readline → fgets fallback
Readline is only used in `repl_loop()` (interactive mode).  
In pipe mode (`pipe_mode_loop`) we already use `fgets` — no change needed.  
For Windows interactive mode, compile with `USE_READLINE=0` and use `fgets`.

### makefile: detect platform
```make
UNAME := $(shell uname -s 2>/dev/null || echo Windows)
ifeq ($(UNAME),Windows_NT)
  LDFLAGS += -lgmp -lm
  CFLAGS  += -DUSE_READLINE=0
else ifeq ($(UNAME),Darwin)
  ...
else
  # Linux
  LDFLAGS += -lgmp -lm -lreadline
endif
```

### build-sidecar.sh → cross-platform
Add `build-sidecar.bat` for Windows (PowerShell).

---

## 6. Build Instructions (docs/building.md outline)

### macOS
```bash
brew install gmp readline
./frontend/build-sidecar.sh
cd frontend && cargo tauri build
```

### Linux (Ubuntu 22.04+)
```bash
sudo apt-get install libgmp-dev libreadline-dev \
  libwebkit2gtk-4.1-dev libgtk-3-dev
./frontend/build-sidecar.sh
cd frontend && cargo tauri build
```

### Windows (MSYS2 + MinGW64)
```powershell
# In MSYS2 MinGW64 shell:
pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-gcc
make USE_ECM=0 USE_READLINE=0
# Then in PowerShell:
cd frontend
.\build-sidecar.bat
cargo tauri build
```

---

## 7. Visual Design Spec (for Track A)

### Colour Palette (Catppuccin Mocha base)
```
Canvas bg:   #050810  (deep space)
Card bg:     rgba(12,15,28,0.85)
Border:      rgba(255,255,255,0.08)
Text:        #cdd6f4
Muted:       #585b70
Accent:      #89b4fa  (blue)

Notebook colours (6-hue palette):
  #89b4fa  blue     #a6e3a1  green
  #f38ba8  red      #fab387  peach
  #cba6f7  mauve    #94e2d5  teal
```

### Motion Spec
| Event | Duration | Easing |
|---|---|---|
| Card materialise | 180ms | spring (stiffness 400, damping 28) |
| Collapse/expand | 220ms | ease-out-quart |
| Pan/zoom spring | per-frame lerp k=0.14 | — |
| Hover glow | 120ms | ease-in-out |

### Typography
- Title bar: `Inter` or `SF Pro`, 13px, weight 500, tracking -0.01em
- Cell input: `JetBrains Mono` or `Fira Code`, 13px
- Canvas label (very zoomed out): 11px, weight 400, `opacity: 0.6`

---

## 8. Estimated Size

| Track | Files | LoC |
|---|---|---|
| A — Canvas UI | 3 new + 2 modified | ~800 |
| B — Library I/O | 2 modified | ~200 |
| C — Example libs | 5 new `.lb` | ~500 lines JSON |
| D — Build + docs | 1 modified C + docs | ~300 |
| **Total** | | **~1800** |

**Timeline:** Tracks A+B+C+D run in parallel. Estimated one session end-to-end.

---

## 9. What's Not In This Phase

- Edge wires between notebooks
- Cross-notebook value references (`Notebook2::Out[3]`)
- Real-time collaborative editing
- Minimap overlay
- Notebook templates/snippets
- Cloud sync
