# Mathilda Canvas — Design & Implementation Plan

## Vision

A single infinite canvas containing multiple named Mathilda notebooks.  
Pan, zoom, and pinch anywhere — every notebook visible simultaneously,  
no separate macOS windows. Next-gen aesthetic: dark glass, blur, glows,  
smooth spring physics.

---

## User Experience

### Canvas Navigation
| Gesture | Action |
|---|---|
| Two-finger drag (trackpad) | Pan |
| Pinch in/out (trackpad) | Zoom |
| Scroll wheel | Zoom (cursor-centred) |
| Click + drag on background | Pan (fallback) |
| Cmd+0 | Reset to fit-all |
| Cmd+= / Cmd+- | Zoom in / out |

### Notebook Cards
Each notebook lives as a **card** on the canvas:

```
╔════════════════════════════╗  ← drag handle
║  ◎ Notebook 1          ⊞ ✕ ║  ← title | collapse | close
╠════════════════════════════╣
║  In[1]:= Integrate[x^2,x] ║
║  Out[1]= x³/3             ║
║  In[2]:= Factor[x^4-1]    ║
║  Out[2]= (x-1)(x+1)(x²+1) ║
╚════════════════════════════╝
```

**Collapsed mode** (zoomed out or manually toggled):
```
╔══════════════════╗
║  ◎ Notebook 1 ⊞  ║
║  4 cells · ready ║
╚══════════════════╝
```

Auto-collapse when zoom < 0.45 (cards become too small to read).  
Auto-expand when zoom >= 0.45.

### Creating Notebooks
- `Cmd+N` → new notebook at centre of viewport
- Double-click empty canvas area → new notebook there
- `+` button in bottom-right corner

---

## Architecture

### Data Flow
```
canvas.ts  ←→  canvasState store  (pan, zoom, [CanvasNotebook])
               each CanvasNotebook has its own NotebookStore
               
Canvas.svelte      — the infinite pan/zoom stage
  └── NotebookCard.svelte  — draggable card (one per notebook)
        └── CellShell.svelte  — receives store as prop (no more singleton)
              └── Output.svelte, KernelStatus, etc.

App.svelte  — thin shell: renders <Canvas /> + global kernel status
```

### Key Design Decisions
1. **Per-notebook store** — `createNotebook()` called once per card. Keeps cell state, outputs, exec indices isolated.
2. **Shared kernel** — all notebooks talk to the single Mathilda process; correlation IDs keep responses routed correctly.
3. **CSS `transform` canvas** — `translate(panX, panY) scale(zoom)` on a single `<div>`. No WebGL, no canvas2D — just fast CSS GPU compositing.
4. **Pointer Events API** — `pointermove`, `pointerdown` for unified mouse/trackpad/touch drag.
5. **`wheel` event for pinch** — on macOS trackpad, pinch generates `wheel` events with `ctrlKey=true` and fractional `deltaY`. This is the standard way to detect pinch in browsers (no need for Touch API on desktop).
6. **Spring physics** — `panX/panY/zoom` smoothed with a lightweight spring (lerp each frame via `requestAnimationFrame`) for that fluid Figma feel.

### Pinch & Two-Finger Pan (macOS Trackpad)
```javascript
canvas.addEventListener('wheel', e => {
  e.preventDefault();
  if (e.ctrlKey) {
    // Pinch-to-zoom: deltaY encodes zoom delta on macOS
    const factor = 1 - e.deltaY * 0.01;
    setZoom(zoom * factor, e.clientX, e.clientY);
  } else {
    // Two-finger pan: deltaX/deltaY are pan deltas
    setPan(panX - e.deltaX, panY - e.deltaY);
  }
}, { passive: false });
```

---

## Visual Design — Next-Gen Dark Glass

### Colour System
```css
/* Canvas background */
--canvas-bg: radial-gradient(ellipse at 30% 30%, #0d1117 0%, #050810 100%);

/* Card chrome */
--card-bg:         rgba(15, 18, 30, 0.82);
--card-border:     rgba(255, 255, 255, 0.08);
--card-glow:       0 0 0 1px rgba(137,180,250,0.12),
                   0 20px 60px rgba(0,0,0,0.6),
                   0 1px 0 rgba(255,255,255,0.06) inset;
--card-blur:       backdrop-filter: blur(24px);

/* Title bar */
--titlebar-bg:     rgba(137,180,250,0.06);
--titlebar-border: rgba(255,255,255,0.06);

/* Accent */
--accent:          #89b4fa;   /* Catppuccin Mocha Blue */
--accent-warm:     #fab387;   /* Peach for running state */
--accent-green:    #a6e3a1;   /* Green for done */
```

### Typography
- Title: `SF Pro Display` / `Inter`, weight 500, tracking -0.02em
- Cell source: `JetBrains Mono` / `Fira Code`, 13px
- Output math: KaTeX rendered, 15px

### Micro-interactions
- Card drag: `box-shadow` deepens + subtle `scale(1.01)` on drag start
- Collapse: height animates with `grid-template-rows` trick (0fr → 1fr)
- New notebook: materialises with scale(0.92→1) + opacity(0→1), 200ms spring
- Zoom: spring-smoothed so it feels physical not instant

### Canvas Grid
Subtle dot grid at `opacity: 0.035` that scales with zoom.  
Dots at 40px spacing (canvas coords), rendered via CSS `background-image: radial-gradient(...)`.

---

## Implementation Order

### Step 1 — Canvas scaffold + pan/zoom ✓ (canvas.ts exists)
- `Canvas.svelte`: wheel handler (pinch + pan), background drag, grid
- Render notebook cards at their `(x, y)` positions
- Spring-smooth pan/zoom state

### Step 2 — NotebookCard.svelte
- Drag by title bar (pointer events)
- Collapse / expand toggle with animation
- Pass own store to child CellShell components
- Title inline-edit (double-click)

### Step 3 — CellShell.svelte → store prop
- Replace singleton `notebook` import with `store` prop
- All mutation calls go through `store.xxx()`

### Step 4 — App.svelte as thin shell
- Just renders `<Canvas />` + floating kernel status + Cmd+N handler

### Step 5 — Auto-collapse at zoom threshold
- Watch zoom; collapse all cards when `zoom < 0.45`
- Expand on zoom in

### Step 6 — Polish
- Fit-all (Cmd+0): compute bounding box of all cards, tween to fit
- Mini-map (optional): tiny overview in bottom-left
- Card color accents (each notebook gets a unique hue)

---

## Estimated Size
**L** — ~600 new lines across 3 new files + refactoring CellShell + App.  
Multi-PR if needed: canvas scaffold → card drag → store refactor → polish.

---

## Out of Scope (This Phase)
- Edges / wires between notebooks
- Cross-notebook value references
- Notebook saving to individual files
- Collaborative multiplayer
- Mobile / touch (pinch via Touch API)
