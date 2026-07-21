<!--
  Canvas.svelte
  Infinite pan/zoom stage. Renders NotebookCard for each notebook in canvasState.
  - Wheel: ctrlKey → zoom centred on cursor; else → two-finger pan
  - Pointer drag on background → pan
  - Spring-smooth animation (lerp k=0.14) via rAF
  - Dot-grid background that moves with pan/zoom
  - Auto-collapse all cards when zoom < 0.40, expand when >= 0.40
  - Cmd+0 → fit all notebooks in viewport
  - Cmd+N → add new notebook at viewport centre
-->
<script lang="ts">
  import { onMount, onDestroy } from 'svelte';
  import { get } from 'svelte/store';
  import NotebookCard from './NotebookCard.svelte';
  import Minimap from './Minimap.svelte';
  import { isTouchDevice } from './platform';
  import {
    canvasState,
    addNotebook,
    addNotebookAt,
    setFocused,
    loadStartupContent,
  } from './canvas';

  // ---------------------------------------------------------------------------
  // Active notebook — the most recently clicked card renders on top (z-index)
  let activeNbId: string | null = null;

  // ---------------------------------------------------------------------------
  // Rubber-band selection — drag on empty canvas to select multiple notebooks

  let selStart: { sx: number; sy: number } | null = null;  // screen coords
  let selCur:   { sx: number; sy: number } | null = null;

  // Compute the selection rect in screen coords for rendering
  $: selRect = selStart && selCur ? {
    x: Math.min(selStart.sx, selCur.sx),
    y: Math.min(selStart.sy, selCur.sy),
    w: Math.abs(selCur.sx - selStart.sx),
    h: Math.abs(selCur.sy - selStart.sy),
  } : null;

  function startSelection(e: PointerEvent) {
    selStart = { sx: e.clientX, sy: e.clientY };
    selCur   = { sx: e.clientX, sy: e.clientY };
    canvasState.update(s => ({ ...s, selectedIds: [] }));
  }

  function updateSelection(e: PointerEvent) {
    if (!selStart) return;
    selCur = { sx: e.clientX, sy: e.clientY };
  }

  function finishSelection() {
    if (!selStart || !selCur) { selStart = selCur = null; return; }
    const rect = canvasEl?.getBoundingClientRect() ?? { left: 0, top: 0 };
    // Convert screen rect to world rect
    const toWorld = (sx: number, sy: number) => ({
      wx: (sx - rect.left - panX) / zoom,
      wy: (sy - rect.top  - panY) / zoom,
    });
    const a = toWorld(selStart.sx, selStart.sy);
    const b = toWorld(selCur.sx,   selCur.sy);
    const minWx = Math.min(a.wx, b.wx), maxWx = Math.max(a.wx, b.wx);
    const minWy = Math.min(a.wy, b.wy), maxWy = Math.max(a.wy, b.wy);
    // Select notebooks whose bounding box overlaps the selection rect
    const ids = get(canvasState).notebooks
      .filter(nb => {
        const h = nb.height ?? 400;
        return nb.x < maxWx && nb.x + nb.width > minWx &&
               nb.y < maxWy && nb.y + h > minWy;
      })
      .map(nb => nb.id);
    canvasState.update(s => ({ ...s, selectedIds: ids }));
    selStart = selCur = null;
  }

  // ---------------------------------------------------------------------------
  // Group drag — when a selected notebook is dragged, move all selected ones

  let groupDragActive = false;
  let groupDragStart: { cx: number; cy: number } | null = null;
  let groupDragOrigins: Map<string, { x: number; y: number }> = new Map();

  function startGroupDrag(e: PointerEvent, nbId: string) {
    if (!$canvasState.selectedIds.includes(nbId)) return;
    groupDragActive = true;
    groupDragStart  = { cx: e.clientX, cy: e.clientY };
    groupDragOrigins = new Map(
      $canvasState.notebooks
        .filter(nb => $canvasState.selectedIds.includes(nb.id))
        .map(nb => [nb.id, { x: nb.x, y: nb.y }])
    );
    canvasEl?.setPointerCapture(e.pointerId);
  }

  function updateGroupDrag(e: PointerEvent) {
    if (!groupDragActive || !groupDragStart) return;
    const dx = (e.clientX - groupDragStart.cx) / zoom;
    const dy = (e.clientY - groupDragStart.cy) / zoom;
    canvasState.update(s => ({
      ...s,
      notebooks: s.notebooks.map(nb => {
        const origin = groupDragOrigins.get(nb.id);
        if (!origin) return nb;
        return { ...nb, x: origin.x + dx, y: origin.y + dy };
      }),
    }));
  }

  function endGroupDrag() {
    groupDragActive = false;
    groupDragStart  = null;
    groupDragOrigins.clear();
  }

  // ---------------------------------------------------------------------------
  // Animated display values — lerped towards store "targets" each rAF tick.

  let panX  = 0;
  let panY  = 0;
  let zoom  = 1.0;

  let targetPanX = 0;
  let targetPanY = 0;
  let targetZoom = 1.0;

  const unsub = canvasState.subscribe(s => {
    targetPanX = s.panX;
    targetPanY = s.panY;
    targetZoom = s.zoom;
  });

  // rAF spring loop
  let rafId: number;
  const K = 0.14;

  function lerp(a: number, b: number, k: number): number { return a + (b - a) * k; }

  function animate() {
    panX = lerp(panX, targetPanX, K);
    panY = lerp(panY, targetPanY, K);
    zoom = lerp(zoom, targetZoom, K);
    rafId = requestAnimationFrame(animate);
  }

  // ---------------------------------------------------------------------------
  // Canvas element ref

  let canvasEl: HTMLElement;
  let worldEl: HTMLElement;

  onMount(() => {
    // Seed display values immediately to avoid lerp-from-zero flash
    panX = targetPanX; panY = targetPanY; zoom = targetZoom;
    rafId = requestAnimationFrame(animate);
    // Load startup content after a tick to ensure all stores are ready
    setTimeout(() => {
      try { loadStartupContent(); } catch (e) { console.error('startup load failed:', e); }
    }, 100);
  });

  onDestroy(() => {
    cancelAnimationFrame(rafId);
    unsub();
  });

  // ---------------------------------------------------------------------------
  // Wheel handler

  function onWheel(e: WheelEvent) {
    if (e.ctrlKey) {
      // Pinch: always zoom the canvas (even over cards). Must preventDefault
      // to stop browser from zooming the page.
      e.preventDefault();
      const factor = 1 - e.deltaY * 0.008;
      if (factor < 1 && $canvasState.focusedId) { setFocused(null); return; }
      const rect = canvasEl?.getBoundingClientRect();
      if (!rect) return;
      const cx = e.clientX - rect.left;
      const cy = e.clientY - rect.top;
      canvasState.update(s => {
        const newZoom = Math.max(0.08, Math.min(3, s.zoom * factor));
        const zf = newZoom / s.zoom;
        return { ...s, zoom: newZoom, panX: cx - zf * (cx - s.panX), panY: cy - zf * (cy - s.panY) };
      });
    } else {
      // Two-finger scroll: scroll the notebook if a cell inside has focus;
      // pan the canvas otherwise (including after pressing Escape).
      const overCard = (e.target as HTMLElement).closest('.nb-card');
      const cellFocused = overCard && document.activeElement && overCard.contains(document.activeElement);
      if (cellFocused) return; // browser scrolls the card natively
      // No focused cell → pan canvas
      e.preventDefault();
      canvasState.update(s => ({ ...s, panX: s.panX - e.deltaX, panY: s.panY - e.deltaY }));
    }
  }

  // ---------------------------------------------------------------------------
  // Pointer drag on background → pan

  let dragging   = false;
  let dragStartX = 0;
  let dragStartY = 0;
  let dragPanX0  = 0;
  let dragPanY0  = 0;

  // Double-click on empty canvas → new notebook at cursor world position
  function onDblClick(e: MouseEvent) {
    if ((e.target as HTMLElement).closest('.nb-card-wrapper, .nb-card, button')) return;
    const rect   = canvasEl?.getBoundingClientRect() ?? { left: 0, top: 0 };
    const worldX = (e.clientX - rect.left - panX) / zoom - 320;
    const worldY = (e.clientY - rect.top  - panY) / zoom - 30;
    addNotebookAt(worldX, worldY);  // notebook appears at cursor, no zoom change
  }

  // ---------------------------------------------------------------------------
  // Touch (coarse pointer): 1-finger pan, 2-finger pinch-zoom.
  //
  // On mouse, panning rides on wheel events and dragging a card moves it. Touch
  // devices fire neither wheel nor a hoverable cursor, so here a single finger
  // on empty canvas pans, two fingers anywhere pan+zoom, and a single finger on
  // a card is left alone (tap to focus, scroll the card body) — card dragging
  // and rubber-band selection are disabled (see NotebookCard: isTouchDevice).
  //
  // We drive the store *and* seed the animated display values so touch tracks
  // the finger 1:1; the rAF lerp (tuned for the mouse-wheel spring) would
  // otherwise add visible lag to direct manipulation.

  const touchPts = new Map<number, { x: number; y: number }>();
  let touchPanning = false;
  let pinchPrevDist = 0;
  let gestPrevX = 0;   // previous anchor (finger, or 2-finger midpoint) in screen px
  let gestPrevY = 0;

  function twoFingerDist(): number {
    const p = [...touchPts.values()];
    return p.length < 2 ? 0 : Math.hypot(p[0].x - p[1].x, p[0].y - p[1].y);
  }
  function twoFingerMid(): { x: number; y: number } {
    const p = [...touchPts.values()];
    return { x: (p[0].x + p[1].x) / 2, y: (p[0].y + p[1].y) / 2 };
  }

  // Commit a new view and snap the animated display to it (no lerp lag).
  function applyView(zoomV: number, px: number, py: number) {
    canvasState.update(s => ({ ...s, zoom: zoomV, panX: px, panY: py }));
    zoom = targetZoom = zoomV;
    panX = targetPanX = px;
    panY = targetPanY = py;
  }

  function onTouchDown(e: PointerEvent) {
    const el = e.target as HTMLElement;
    const overCard = el.closest('.nb-card');
    const overInteractive = el.closest(
      'button, input, a, [role="button"], .cm-editor, [contenteditable="true"]'
    );
    touchPts.set(e.pointerId, { x: e.clientX, y: e.clientY });

    if (touchPts.size >= 2) {
      // Second finger down → begin pinch-zoom/pan (works even over a card).
      pinchPrevDist = twoFingerDist();
      const m = twoFingerMid();
      gestPrevX = m.x; gestPrevY = m.y;
      touchPanning = true;
      canvasEl.setPointerCapture(e.pointerId);
      e.preventDefault();
      return;
    }

    // Single finger over a card / control → let the card handle it (tap, scroll,
    // edit). No canvas pan, no card drag.
    if (overCard || overInteractive) return;

    // Single finger on empty canvas → pan.
    touchPanning = true;
    gestPrevX = e.clientX; gestPrevY = e.clientY;
    canvasEl.setPointerCapture(e.pointerId);
  }

  function onTouchMove(e: PointerEvent) {
    if (!touchPts.has(e.pointerId)) return;
    touchPts.set(e.pointerId, { x: e.clientX, y: e.clientY });
    if (!touchPanning) return;

    if (touchPts.size >= 2) {
      const dist = twoFingerDist();
      const m = twoFingerMid();
      const rect = canvasEl.getBoundingClientRect();
      const cx = m.x - rect.left, cy = m.y - rect.top;
      if (pinchPrevDist > 0 && dist > 0) {
        const factor = dist / pinchPrevDist;
        const newZoom = Math.max(0.08, Math.min(3, zoom * factor));
        const zf = newZoom / zoom;
        const mdx = m.x - gestPrevX, mdy = m.y - gestPrevY;   // midpoint travel → pan
        applyView(newZoom, cx - zf * (cx - panX) + mdx, cy - zf * (cy - panY) + mdy);
      }
      pinchPrevDist = dist;
      gestPrevX = m.x; gestPrevY = m.y;
      e.preventDefault();
      return;
    }

    // Single-finger pan.
    const dx = e.clientX - gestPrevX, dy = e.clientY - gestPrevY;
    gestPrevX = e.clientX; gestPrevY = e.clientY;
    applyView(zoom, panX + dx, panY + dy);
    e.preventDefault();
  }

  function onTouchUp(e: PointerEvent) {
    touchPts.delete(e.pointerId);
    if (touchPts.size < 2) pinchPrevDist = 0;
    if (touchPts.size === 1) {
      // Dropped from two fingers to one → re-anchor pan on the remaining finger
      // so the view doesn't jump.
      const p = [...touchPts.values()][0];
      gestPrevX = p.x; gestPrevY = p.y;
    } else if (touchPts.size === 0) {
      touchPanning = false;
    }
  }

  function onPointerDown(e: PointerEvent) {
    if (isTouchDevice) { onTouchDown(e); return; }
    if (e.button !== 0) return;
    const overCard = (e.target as HTMLElement).closest('.nb-card');
    // Interactive elements inside cards — don't interfere
    if ((e.target as HTMLElement).closest('button, input, a, [role="button"]')) return;

    if (overCard) {
      // Clicking a card: if it's already selected, start group drag
      const wrapper = (e.target as HTMLElement).closest('[data-nb-id]') as HTMLElement | null;
      if (wrapper) {
        const nbId = wrapper.getAttribute('data-nb-id');
        if (nbId) {
          activeNbId = nbId;
          if ($canvasState.selectedIds.includes(nbId)) {
            startGroupDrag(e, nbId);
            return;
          }
        }
      }
      return; // let NotebookCard handle normal single drag
    }

    // Empty canvas — plain drag draws a rubber-band selection rect.
    // Two-finger pan and pinch-zoom are handled by onWheel (no pointer drag needed for pan).
    startSelection(e);
    canvasEl.setPointerCapture(e.pointerId);
  }

  function onPointerMove(e: PointerEvent) {
    if (isTouchDevice) { onTouchMove(e); return; }
    if (groupDragActive) { updateGroupDrag(e); return; }
    if (selStart) { updateSelection(e); return; }
    if (!dragging) return;
    canvasState.update(s => ({
      ...s,
      panX: dragPanX0 + (e.clientX - dragStartX),
      panY: dragPanY0 + (e.clientY - dragStartY),
    }));
  }

  function onPointerUp(_e: PointerEvent) {
    if (isTouchDevice) { onTouchUp(_e); return; }
    if (groupDragActive) { endGroupDrag(); return; }
    if (selStart) { finishSelection(); return; }
    dragging = false;
  }

  // ---------------------------------------------------------------------------
  // Auto-collapse when zoomed out past threshold

  // Auto-collapse removed — notebooks stay expanded at all zoom levels.

  // ---------------------------------------------------------------------------
  // Keyboard: Cmd+0 fit-all (Cmd+N is macOS "new window" and can't be reliably intercepted)

  function onKeydown(e: KeyboardEvent) {
    const mod = e.metaKey || e.ctrlKey;
    if (mod && e.key === '0') { e.preventDefault(); fitAll(); return; }
    // 'N' key (no modifier) when not in an editor → add notebook at centre
    if (!mod && !e.shiftKey && e.key === 'n' && !$canvasState.focusedId) {
      const target = e.target as HTMLElement;
      if (!target.closest('.cm-editor') && !target.isContentEditable) {
        e.preventDefault();
        addAtCentre();
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Right-click = add notebook at cursor world position
  function onContextMenu(e: MouseEvent) {
    // On touch, contextmenu fires on long-press and would spawn a notebook
    // mid-pan. Suppress it there; the "＋ New Notebook" button and double-tap
    // remain as add paths.
    if (isTouchDevice) { e.preventDefault(); return; }
    if ((e.target as HTMLElement).closest('.nb-card-wrapper, .nb-card, button')) return;
    e.preventDefault();
    const rect   = canvasEl?.getBoundingClientRect() ?? { left: 0, top: 0 };
    const worldX = (e.clientX - rect.left - panX) / zoom - 320;
    const worldY = (e.clientY - rect.top  - panY) / zoom - 30;
    addNotebookAt(worldX, worldY);
  }

  function fitAll() {
    const nbs = get(canvasState).notebooks;
    if (!nbs.length) return;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const nb of nbs) {
      const h = nb.collapsed ? 52 : 520;
      minX = Math.min(minX, nb.x);
      minY = Math.min(minY, nb.y);
      maxX = Math.max(maxX, nb.x + nb.width);
      maxY = Math.max(maxY, nb.y + h);
    }
    const PAD = 60;
    const bbW = maxX - minX + PAD * 2;
    const bbH = maxY - minY + PAD * 2;
    const vw  = canvasEl?.clientWidth  || window.innerWidth;
    const vh  = canvasEl?.clientHeight || window.innerHeight;
    const newZoom = Math.max(0.08, Math.min(3, Math.min(vw / bbW, vh / bbH)));
    const newPanX = (vw - bbW * newZoom) / 2 - (minX - PAD) * newZoom;
    const newPanY = (vh - bbH * newZoom) / 2 - (minY - PAD) * newZoom;
    canvasState.update(s => ({ ...s, panX: newPanX, panY: newPanY, zoom: newZoom }));
  }

  function addAtCentre() {
    const vw = canvasEl?.clientWidth  || window.innerWidth;
    const vh = canvasEl?.clientHeight || window.innerHeight;
    const n  = get(canvasState).notebooks.length;
    // World coords for the current viewport centre
    const cx = (vw / 2 - targetPanX) / targetZoom - 320 + n * 20;
    const cy = (vh / 2 - targetPanY) / targetZoom - 100 + n * 20;
    addNotebook();
    canvasState.update(s => {
      const nbs  = s.notebooks.slice();
      nbs[nbs.length - 1] = { ...nbs[nbs.length - 1], x: cx, y: cy };
      return { ...s, notebooks: nbs };
    });
  }

  // ---------------------------------------------------------------------------
  // Minimap click: pan+zoom to the clicked notebook

  function onMinimapNotebookClick(nb: import('./canvas').CanvasNotebook) {
    const vw = canvasEl?.clientWidth  || window.innerWidth;
    const vh = canvasEl?.clientHeight || window.innerHeight;
    const nbH = nb.height ?? 420;
    const pad = 60;
    const newZoom = Math.min(1.2,
      (vw - pad * 2) / nb.width,
      (vh - pad * 2) / nbH
    );
    const newPanX = (vw - nb.width * newZoom) / 2 - nb.x * newZoom;
    const newPanY = (vh - nbH   * newZoom) / 2 - nb.y * newZoom;
    canvasState.update(s => ({ ...s, panX: newPanX, panY: newPanY, zoom: newZoom }));
  }

  // ---------------------------------------------------------------------------
  // Dot grid (reacts to animated pan/zoom)

  $: dotSpacing = 28 * zoom;
  $: gridBgStyle =
    `background-image: radial-gradient(circle, rgba(255,255,255,0.07) 1px, transparent 1px);` +
    `background-size: ${dotSpacing}px ${dotSpacing}px;` +
    `background-position: ${panX % dotSpacing}px ${panY % dotSpacing}px;`;

</script>

<svelte:window on:keydown={onKeydown} />

{#if $canvasState.focusedId}
  <!-- ── Focused (full-screen) mode — pinch out to return ── -->
  {@const fnb = $canvasState.notebooks.find(n => n.id === $canvasState.focusedId)}
  {#if fnb}
    <div
      class="focused-view"
    >
      <div class="focused-view-inner">
        <NotebookCard nb={fnb} currentZoom={1} focused={true} />
      </div>
    </div>
  {/if}
{:else}
  <!-- ── Canvas mode ── -->
  <!-- svelte-ignore a11y-no-static-element-interactions -->
  <div
    class="canvas-stage"
    style={gridBgStyle}
    bind:this={canvasEl}
    on:wheel|nonpassive={onWheel}
    on:pointerdown={onPointerDown}
    on:pointermove={onPointerMove}
    on:pointerup={onPointerUp}
    on:pointercancel={onPointerUp}
    on:dblclick={onDblClick}
    on:contextmenu={onContextMenu}
  >
    <div
      class="canvas-world"
      bind:this={worldEl}
      style="transform: translate({panX}px, {panY}px) scale({zoom}); transform-origin: 0 0;"
    >
      {#each $canvasState.notebooks as nb (nb.id)}
        <!-- svelte-ignore a11y-no-static-element-interactions -->
        <div
          class="nb-card-wrapper"
          data-nb-id={nb.id}
          style="left:{nb.x}px;top:{nb.y}px;width:{nb.width}px;z-index:{activeNbId===nb.id?10:1};"
          on:pointerdown|capture={() => activeNbId = nb.id}
        >
          <NotebookCard
            {nb}
            currentZoom={zoom}
            isSelected={$canvasState.selectedIds.includes(nb.id)}
            on:focusNotebook={(e) => setFocused(e.detail.id)}
            on:groupMoveEnd={() => groupDragOrigins.clear()}
            on:groupMove={(e) => {
              const { dx, dy, originX, originY, id } = e.detail;
              // Lazy-init origins for all other selected notebooks on first move
              if (groupDragOrigins.size === 0) {
                const s = get(canvasState);
                s.notebooks
                  .filter(n => s.selectedIds.includes(n.id) && n.id !== id)
                  .forEach(n => groupDragOrigins.set(n.id, { x: n.x, y: n.y }));
              }
              canvasState.update(s => ({
                ...s,
                notebooks: s.notebooks.map(n => {
                  if (!s.selectedIds.includes(n.id)) return n;
                  if (n.id === id) return { ...n, x: originX + dx, y: originY + dy };
                  const origin = groupDragOrigins.get(n.id);
                  if (!origin) return n;
                  return { ...n, x: origin.x + dx, y: origin.y + dy };
                }),
              }));
            }}
          />
        </div>
      {/each}
    </div>

    <div class="canvas-hints">
      <button class="hint-new-btn" on:click={addAtCentre}>＋ New Notebook</button>
      <span class="hint-sep">·</span>
      <span>dbl-click or right-click canvas</span>
      <span class="hint-sep">·</span>
      <span>⌘0 fit · N key</span>
      <span class="hint-sep">·</span>
      <span>scroll pan · pinch zoom</span>
    </div>

    <!-- Rubber-band selection rectangle -->
    {#if selRect && selRect.w > 4 && selRect.h > 4}
      <div
        class="sel-rect"
        style="left:{selRect.x}px;top:{selRect.y}px;width:{selRect.w}px;height:{selRect.h}px;"
      ></div>
    {/if}
  </div>

  <!-- Minimap: click a notebook rect to jump to it -->
  <Minimap
    notebooks={$canvasState.notebooks}
    {panX} {panY} {zoom}
    viewportW={window.innerWidth}
    viewportH={window.innerHeight}
    onNotebookClick={onMinimapNotebookClick}
  />
{/if}

<style>
  .canvas-stage {
    position: fixed;
    inset: 0;
    width: 100vw;
    height: 100vh;
    background-color: var(--bg, #050810);
    overflow: hidden;
    cursor: default;
    /* prevent text selection while dragging */
    user-select: none;
    -webkit-user-select: none;
    /* Take over touch gestures (finger pan / pinch-zoom) instead of letting the
       browser scroll or page-zoom. Card bodies re-enable vertical scrolling via
       their own touch-action (see NotebookCard .card-body). */
    touch-action: none;
  }

  .canvas-world {
    position: absolute;
    top: 0;
    left: 0;
    width: 0;
    height: 0;
    will-change: transform;
  }

  .nb-card-wrapper {
    position: absolute;
  }

  /* Rubber-band selection rectangle */
  .sel-rect {
    position: absolute;
    pointer-events: none;
    border: 1.5px dashed var(--accent, #89b4fa);
    background: rgba(137, 180, 250, 0.08);
    border-radius: 4px;
    z-index: 9999;
  }

  /* Canvas hints must NOT scale with Cmd+/- — use px not rem */
  .canvas-hints {
    position: fixed;
    /* Sit above the device's bottom gesture/nav bar on mobile (0 on desktop). */
    bottom: calc(14px + env(safe-area-inset-bottom, 0px));
    left: 50%;
    transform: translateX(-50%);
    display: flex;
    align-items: center;
    gap: 10px;
    font-size: 11px;   /* fixed px — immune to root font-size changes */
    color: rgba(255,255,255,0.25);
    pointer-events: none;
    letter-spacing: 0.02em;
    white-space: nowrap;
  }
  .hint-new-btn {
    background: rgba(137,180,250,0.15);
    border: 1px solid rgba(137,180,250,0.35);
    color: rgba(137,180,250,0.9);
    border-radius: 5px;
    padding: 4px 10px;
    font-size: 12px;   /* fixed px */
    cursor: pointer;
    pointer-events: auto;
    transition: background 0.12s;
    letter-spacing: 0.01em;
    white-space: nowrap;
  }
  .hint-new-btn:hover { background: rgba(137,180,250,0.25); }
  .hint-sep { opacity: 0.4; }

  /* ---- Floating add button ---- */

  /* ---- Focused (full-screen) view — truly edge to edge ---- */
  .focused-view {
    position: fixed;
    inset: 0;
    /* Use card-bg so light mode doesn't show dark canvas edges */
    background: var(--card-bg, #050810);
    overflow-y: auto;
    z-index: 50;
  }
  .focused-view-inner {
    width: 100%;
    /* No max-width, no side padding — notebook card fills the window */
  }
  /* Override card styles when in focused view so it has no border-radius or side margins */
  .focused-view-inner :global(.nb-card) {
    border-radius: 0;
    border-left: none;
    border-right: none;
    min-height: 100vh;
  }

</style>
