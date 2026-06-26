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
  import {
    canvasState,
    addNotebook,
  } from './canvas';

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
  });

  onDestroy(() => {
    cancelAnimationFrame(rafId);
    unsub();
  });

  // ---------------------------------------------------------------------------
  // Wheel handler

  function onWheel(e: WheelEvent) {
    e.preventDefault();
    if (e.ctrlKey) {
      // Pinch-to-zoom centred on cursor
      const factor  = 1 - e.deltaY * 0.008;
      const rect    = canvasEl.getBoundingClientRect();
      const cx      = e.clientX - rect.left;
      const cy      = e.clientY - rect.top;
      canvasState.update(s => {
        const newZoom  = Math.max(0.08, Math.min(3, s.zoom * factor));
        const zFactor  = newZoom / s.zoom;
        return {
          ...s,
          zoom:  newZoom,
          panX:  cx - zFactor * (cx - s.panX),
          panY:  cy - zFactor * (cy - s.panY),
        };
      });
    } else {
      // Two-finger scroll → translate
      canvasState.update(s => ({
        ...s,
        panX: s.panX - e.deltaX,
        panY: s.panY - e.deltaY,
      }));
    }
    maybeAutoCollapse();
  }

  // ---------------------------------------------------------------------------
  // Pointer drag on background → pan

  let dragging   = false;
  let dragStartX = 0;
  let dragStartY = 0;
  let dragPanX0  = 0;
  let dragPanY0  = 0;

  function onPointerDown(e: PointerEvent) {
    if ((e.target as HTMLElement).closest('.nb-card')) return;
    dragging   = true;
    dragStartX = e.clientX;
    dragStartY = e.clientY;
    dragPanX0  = targetPanX;
    dragPanY0  = targetPanY;
    canvasEl.setPointerCapture(e.pointerId);
  }

  function onPointerMove(e: PointerEvent) {
    if (!dragging) return;
    canvasState.update(s => ({
      ...s,
      panX: dragPanX0 + (e.clientX - dragStartX),
      panY: dragPanY0 + (e.clientY - dragStartY),
    }));
  }

  function onPointerUp(_e: PointerEvent) {
    dragging = false;
  }

  // ---------------------------------------------------------------------------
  // Auto-collapse when zoomed out past threshold

  const COLLAPSE_THRESHOLD = 0.40;
  let prevAutoCollapsed: boolean | null = null;

  function maybeAutoCollapse() {
    const curZoom = get(canvasState).zoom;
    const shouldCollapse = curZoom < COLLAPSE_THRESHOLD;
    if (shouldCollapse === prevAutoCollapsed) return;
    prevAutoCollapsed = shouldCollapse;
    canvasState.update(s => ({
      ...s,
      notebooks: s.notebooks.map(nb => ({ ...nb, collapsed: shouldCollapse })),
    }));
  }

  // ---------------------------------------------------------------------------
  // Cmd+0 — fit all / Cmd+N — new notebook at centre

  function onKeydown(e: KeyboardEvent) {
    const mod = e.metaKey || e.ctrlKey;
    if (!mod) return;
    if (e.key === '0') { e.preventDefault(); fitAll(); }
    else if (e.key === 'n' || e.key === 'N') { e.preventDefault(); addAtCentre(); }
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
  // Dot grid (reacts to animated pan/zoom)

  $: dotSpacing = 28 * zoom;
  $: gridBgStyle =
    `background-image: radial-gradient(circle, rgba(255,255,255,0.07) 1px, transparent 1px);` +
    `background-size: ${dotSpacing}px ${dotSpacing}px;` +
    `background-position: ${panX % dotSpacing}px ${panY % dotSpacing}px;`;

</script>

<svelte:window on:keydown={onKeydown} />

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
>
  <div
    class="canvas-world"
    bind:this={worldEl}
    style="transform: translate({panX}px, {panY}px) scale({zoom}); transform-origin: 0 0;"
  >
    {#each $canvasState.notebooks as nb (nb.id)}
      <div
        class="nb-card-wrapper"
        style="left: {nb.x}px; top: {nb.y}px; width: {nb.width}px;"
      >
        <NotebookCard {nb} currentZoom={zoom} />
      </div>
    {/each}
  </div>

  <!-- Keyboard hints overlay -->
  <div class="canvas-hints">
    <span>⌘N new</span>
    <span>⌘0 fit</span>
    <span>scroll pan · ctrl+scroll zoom</span>
  </div>
</div>

<style>
  .canvas-stage {
    position: fixed;
    inset: 0;
    width: 100vw;
    height: 100vh;
    /* deep space base; dot grid applied via inline style */
    background-color: #050810;
    overflow: hidden;
    cursor: default;
    /* prevent text selection while dragging */
    user-select: none;
    -webkit-user-select: none;
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

  .canvas-hints {
    position: fixed;
    bottom: 14px;
    left: 50%;
    transform: translateX(-50%);
    display: flex;
    gap: 1.2rem;
    font-size: 0.72rem;
    color: rgba(255,255,255,0.22);
    pointer-events: none;
    letter-spacing: 0.02em;
  }
</style>
