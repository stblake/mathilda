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
  import type { FocusLayout } from './canvas';
  import NotebookCard from './NotebookCard.svelte';
  import Minimap from './Minimap.svelte';
  import StatusBar from './StatusBar.svelte';
  import { showStatusBar } from './status';
  import { isTouchDevice } from './platform';
  import {
    canvasState,
    addNotebook,
    addNotebookAt,
    setFocused,
    loadStartupContent,
    setStageOrigin,
    removePane,
    addPaneToFocus,
    setFocusLayout,
    splitWithNext,
    openPanes,
    MAX_PANES,
    setFocusedActive,
    setFocusedSizes,
    setFocusedGrid,
  } from './canvas';

  // ---------------------------------------------------------------------------
  // Active notebook — the most recently clicked card renders on top (z-index).
  // Held in canvasState so code outside this component can raise a card;
  // openRefpage does, so a documentation page opens in front of the notebook it
  // was asked from rather than behind it.
  function setActive(id: string) {
    canvasState.update(s => (s.activeId === id ? s : { ...s, activeId: id }));
  }

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

  /* Publish where the stage sits, so a viewport point (a clicked symbol) can be
     turned into a canvas coordinate. Re-read on resize because the app bar and
     window chrome move it. */
  function publishStageOrigin() {
    const r = canvasEl?.getBoundingClientRect();
    if (r) setStageOrigin(r.left, r.top);
  }

  onMount(() => {
    // Seed display values immediately to avoid lerp-from-zero flash
    panX = targetPanX; panY = targetPanY; zoom = targetZoom;
    publishStageOrigin();
    window.addEventListener('resize', publishStageOrigin);
    rafId = requestAnimationFrame(animate);
    // Load startup content after a tick to ensure all stores are ready
    setTimeout(() => {
      try { loadStartupContent(); } catch (e) { console.error('startup load failed:', e); }
    }, 100);
  });

  onDestroy(() => {
    cancelAnimationFrame(rafId);
    window.removeEventListener('resize', publishStageOrigin);
    unsub();
  });

  // ---------------------------------------------------------------------------
  // Wheel handler

  /* Is there an element between `from` and `stop` (inclusive) that can still
     scroll in the direction of this gesture?
     
     "Can still scroll" matters as much as "is scrollable": a card scrolled to
     its bottom should hand the rest of the gesture to the canvas rather than
     swallowing it, which is how nested scrolling behaves everywhere else. The
     1px tolerance absorbs fractional scroll positions at fractional zoom. */
  function canScroll(from: HTMLElement, stop: Element, dx: number, dy: number): boolean {
    let el: HTMLElement | null = from;
    while (el) {
      const style = getComputedStyle(el);
      const scrollableY = /(auto|scroll)/.test(style.overflowY);
      const scrollableX = /(auto|scroll)/.test(style.overflowX);
      if (scrollableY && el.scrollHeight > el.clientHeight + 1) {
        const room = dy > 0
          ? el.scrollHeight - el.clientHeight - el.scrollTop > 1
          : el.scrollTop > 1;
        if (dy !== 0 && room) return true;
      }
      if (scrollableX && el.scrollWidth > el.clientWidth + 1) {
        const room = dx > 0
          ? el.scrollWidth - el.clientWidth - el.scrollLeft > 1
          : el.scrollLeft > 1;
        if (dx !== 0 && room) return true;
      }
      if (el === stop) break;
      el = el.parentElement;
    }
    return false;
  }

  /* Scroll chaining, gesture-aware.
   *
   * A trackpad flick keeps delivering wheel events long after the fingers lift.
   * Without this, reaching the bottom of a notebook mid-flick handed the rest of
   * the momentum to the canvas and the view shot away. Native scrolling avoids
   * that because a scroll gesture is bound to the element it started in.
   *
   * Same rule here: events closer together than GESTURE_GAP_MS are one gesture,
   * and a gesture that has scrolled a card can never pan the canvas -- it stops
   * dead at the boundary. Pausing ends the gesture, so scrolling again from a
   * stop does pan, which is how you get out of a card deliberately. */
  const GESTURE_GAP_MS = 160;
  let lastWheelAt = 0;
  let gestureScrolledCard = false;

  function onWheel(e: WheelEvent) {
    if (e.ctrlKey) {
      // Pinch: always zoom the canvas (even over cards). Must preventDefault
      // to stop browser from zooming the page.
      e.preventDefault();
      const factor = 1 - e.deltaY * 0.008;
      /* Pinch out leaves focused mode entirely, however many panes are open —
         it is the one gesture that means "zoom out of this whole thing". */
      if (factor < 1 && $canvasState.focusedIds.length) { setFocused(null); return; }
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
      // Two-finger scroll: let a card consume it only if the card can actually
      // scroll the way the gesture is going; otherwise pan the canvas.
      //
      // This used to defer to the card whenever a cell inside it held focus,
      // which meant that clicking into a notebook killed panning over that card
      // entirely -- on a card with nothing to scroll the gesture did nothing at
      // all. Focus is the wrong question; scrollability is the right one, and it
      // is what the browser itself uses to chain a scroll to an ancestor.
      const now = e.timeStamp || performance.now();
      if (now - lastWheelAt > GESTURE_GAP_MS) gestureScrolledCard = false;  // new gesture
      lastWheelAt = now;

      const overCard = (e.target as HTMLElement).closest('.nb-card');
      if (overCard && canScroll(e.target as HTMLElement, overCard, e.deltaX, e.deltaY)) {
        gestureScrolledCard = true;
        return;                       // the browser scrolls the card natively
      }
      if (gestureScrolledCard) {
        // Boundary reached mid-flick: absorb the momentum rather than pan.
        e.preventDefault();
        return;
      }
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

  // Route by the actual pointer that fired, not a device-wide flag: a mouse on
  // a touchscreen laptop still gets the mouse model (drag / rubber-band), while
  // a finger anywhere gets pan / pinch. On phones every event is `touch`.
  const isTouchEvent = (e: PointerEvent) => e.pointerType === 'touch' || e.pointerType === 'pen';

  function onPointerDown(e: PointerEvent) {
    if (isTouchEvent(e)) { onTouchDown(e); return; }
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
          setActive(nbId);
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
    if (isTouchEvent(e)) { onTouchMove(e); return; }
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
    if (isTouchEvent(_e)) { onTouchUp(_e); return; }
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

    /* Splitting from the keyboard. Works from inside a cell -- these are the
       shortcuts you reach for mid-edit, so requiring the caret to be outside an
       editor first would defeat them. Backslash rather than a letter because
       every plain Cmd+letter is either taken here or by the native menu. */
    if (mod && (e.key === '\\' || e.code === 'Backslash')) {
      e.preventDefault();
      const layout: FocusLayout = e.shiftKey ? 'v' : 'h';
      if ($canvasState.focusedIds.length === 0) {
        /* From the canvas: open the selection if there is one, else the card on
           top, and split it with the next notebook. */
        const sel = $canvasState.selectedIds;
        if (sel.length >= 2) { openPanes(sel, layout); return; }
        const seed = sel[0] ?? $canvasState.activeId ?? $canvasState.notebooks[0]?.id;
        if (!seed) return;
        setFocused(seed);
      }
      splitWithNext(layout);
      return;
    }

    /* Open the rubber-band selection side by side. Enter is safe here only
       because this branch requires the caret to be outside any editor. */
    if (!mod && e.key === 'Enter' && !$canvasState.focusedIds.length
        && $canvasState.selectedIds.length >= 2) {
      const target = e.target as HTMLElement;
      if (!target.closest('.cm-editor') && !target.isContentEditable) {
        e.preventDefault();
        openPanes($canvasState.selectedIds, 'h');
        return;
      }
    }

    // 'N' key (no modifier) when not in an editor → add notebook at centre
    if (!mod && !e.shiftKey && e.key === 'n' && !$canvasState.focusedIds.length) {
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
  // Tiled focused view

  let focusRootEl: HTMLElement;

  $: fIds    = $canvasState.focusedIds;
  $: fLayout = $canvasState.focusedLayout;
  /* filter(Boolean) guards the frame between a notebook being removed and
     normalizeFocus dropping its id. */
  $: fPanes  = fIds
       .map(id => $canvasState.notebooks.find(n => n.id === id))
       .filter(Boolean) as import('./canvas').CanvasNotebook[];

  /** Divider thickness in px. Its own grid track, so it takes real space rather
   *  than overlapping a pane. */
  const DIV = 10;
  /** A pane narrower than this is not usefully editable. */
  const MIN_PANE_PX = 240;

  /* Sizes come from drag-local state while a drag is in flight, so a 60Hz
     pointermove does not write the canvas store 60 times a second and
     invalidate every card, cell and the minimap along with it. */
  let dragKind: number | 'gx' | 'gy' | null = null;
  let dragSizes: number[] | null = null;
  let dragGrid: { x: number; y: number } | null = null;

  $: sizes   = dragSizes ?? $canvasState.focusedSizes;
  $: gridPos = dragGrid  ?? $canvasState.focusedGrid;

  /* Panes and dividers are both grid items: N sizes interleaved with N-1 fixed
     tracks. Using fr means the fixed px are subtracted automatically instead of
     pushing the total past 100%. */
  $: flowTemplate = sizes.length
    ? sizes.map(p => `minmax(0, ${Math.max(p, 0.001)}fr)`).join(` ${DIV}px `)
    : 'minmax(0, 1fr)';

  $: containerStyle =
      fLayout === 'grid'
        ? `grid-template-columns: minmax(0, ${gridPos.x}fr) ${DIV}px minmax(0, ${100 - gridPos.x}fr);` +
          `grid-template-rows: minmax(0, ${gridPos.y}fr) ${DIV}px minmax(0, ${100 - gridPos.y}fr);`
      : fLayout === 'h'
        ? `grid-template-columns: ${flowTemplate}; grid-template-rows: minmax(0, 1fr);`
        : `grid-template-rows: ${flowTemplate}; grid-template-columns: minmax(0, 1fr);`;

  /** Which grid cell a pane occupies. With three panes the third spans the
   *  bottom row, so the layout has no empty quadrant. */
  function gridArea(i: number, n: number): string {
    if (i === 0) return 'grid-column: 1; grid-row: 1;';
    if (i === 1) return 'grid-column: 3; grid-row: 1;';
    if (i === 2) return n === 3 ? 'grid-column: 1 / -1; grid-row: 3;' : 'grid-column: 1; grid-row: 3;';
    return 'grid-column: 3; grid-row: 3;';
  }

  /* Activation is pointerdown in the CAPTURE phase, mirroring how canvas cards
     are raised: capture runs before CodeMirror's own mouse handling and before
     the several |stopPropagation handlers inside a card, so it cannot be
     swallowed. `focusin` covers the keyboard and programmatic paths (Tab, or the
     card moving the caret itself), since focus does not bubble but focusin does.
     Deliberately not hover: with four panes, a mouse travelling to the toolbar
     would re-target it en route, and Run would fire on the wrong notebook. */
  function activatePane(id: string) {
    if (get(canvasState).focusedActiveId === id) return;   // no store write per click
    setFocusedActive(id);
  }

  function clamp(v: number, lo: number, hi: number) { return Math.max(lo, Math.min(hi, v)); }

  function onDividerDown(e: PointerEvent, kind: number | 'gx' | 'gy') {
    if (e.button !== 0) return;
    e.preventDefault();
    e.stopPropagation();
    dragKind  = kind;
    dragSizes = [...get(canvasState).focusedSizes];
    dragGrid  = { ...get(canvasState).focusedGrid };
    (e.currentTarget as HTMLElement).setPointerCapture(e.pointerId);
  }

  function onDividerMove(e: PointerEvent) {
    if (dragKind === null || !focusRootEl) return;
    e.stopPropagation();
    const r = focusRootEl.getBoundingClientRect();
    const alongX = dragKind === 'gx' || (typeof dragKind === 'number' && fLayout === 'h');
    const extent = alongX ? r.width : r.height;
    if (extent <= 0) return;
    const pos    = alongX ? e.clientX - r.left : e.clientY - r.top;
    const pct    = (pos / extent) * 100;
    const minPct = (MIN_PANE_PX / extent) * 100;

    if (dragKind === 'gx')      dragGrid = { ...gridPos, x: clamp(pct, minPct, 100 - minPct) };
    else if (dragKind === 'gy') dragGrid = { ...gridPos, y: clamp(pct, minPct, 100 - minPct) };
    else {
      /* Divider i moves only the boundary between pane i and i+1: the two trade
         against each other, so the total stays 100 with no renormalise pass and
         the panes you are not touching hold still. */
      const i = dragKind;
      const base = dragSizes ?? [];
      const before = base.slice(0, i).reduce((a, b) => a + b, 0);
      const pair   = (base[i] ?? 0) + (base[i + 1] ?? 0);
      const want   = clamp(pct - before, Math.min(minPct, pair / 2), Math.max(pair - minPct, pair / 2));
      const next = [...base];
      next[i] = want;
      next[i + 1] = pair - want;
      dragSizes = next;
    }
  }

  function onDividerUp() {
    if (dragKind === null) return;
    if (dragKind === 'gx' || dragKind === 'gy') setFocusedGrid(dragGrid!);
    else if (dragSizes) setFocusedSizes(dragSizes);
    dragKind = null; dragSizes = null; dragGrid = null;
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

{#if $canvasState.focusedIds.length}
  <!-- ── Focused mode — 1 to 4 notebooks tiled. Pinch out to return ──
       A CSS grid, not nested flex: flex handles side-by-side and stacked, but the
       2x2 would need two levels of nesting, two kinds of divider, and a
       special case for the 3-pane row that spans. One grid template does all
       three. Tracks are minmax(0, Nfr) so a pane can shrink below its content's
       intrinsic width and scroll, instead of pushing the grid past the window. -->
  <div
    class="focused-view"
    class:with-status={$showStatusBar}
    class:multi={fPanes.length > 1}
    class:dragging={dragKind !== null}
    style={containerStyle}
    bind:this={focusRootEl}
  >
    {#each fPanes as fnb, i (fnb.id)}
      <!-- svelte-ignore a11y-no-static-element-interactions -->
      <div
        class="focused-pane"
        class:active={fPanes.length > 1 && fnb.id === $canvasState.focusedActiveId}
        style={fLayout === 'grid' ? gridArea(i, fPanes.length) : ''}
        data-pane-id={fnb.id}
        on:pointerdown|capture={() => activatePane(fnb.id)}
        on:focusin={() => activatePane(fnb.id)}
      >
        {#if fPanes.length > 1}
          <div class="pane-header">
            <span class="pane-title">{fnb.title}</span>
            <button
              class="pane-close"
              title="Remove this pane (the notebook stays on the canvas)"
              tabindex="-1"
              on:pointerdown|stopPropagation|preventDefault
              on:click|stopPropagation={() => removePane(fnb.id)}
            >✕</button>
          </div>
        {/if}
        <div class="focused-view-inner">
          <NotebookCard nb={fnb} currentZoom={1} focused={true} />
        </div>
      </div>

      <!-- One divider after every pane but the last, along the flow axis. -->
      {#if fLayout !== 'grid' && i < fPanes.length - 1}
        <!-- svelte-ignore a11y-no-static-element-interactions -->
        <div
          class="pane-divider"
          class:vertical={fLayout === 'v'}
          on:pointerdown={(e) => onDividerDown(e, i)}
          on:pointermove={onDividerMove}
          on:pointerup={onDividerUp}
          on:pointercancel={onDividerUp}
        ></div>
      {/if}
    {/each}

    {#if fLayout === 'grid'}
      <!-- The two crossing dividers. The column one is drawn second and takes
           the intersection, so a drag there moves one axis rather than both. -->
      <!-- svelte-ignore a11y-no-static-element-interactions -->
      <div
        class="pane-divider vertical grid-row"
        on:pointerdown={(e) => onDividerDown(e, 'gy')}
        on:pointermove={onDividerMove}
        on:pointerup={onDividerUp}
        on:pointercancel={onDividerUp}
      ></div>
      <!-- svelte-ignore a11y-no-static-element-interactions -->
      <div
        class="pane-divider grid-col"
        on:pointerdown={(e) => onDividerDown(e, 'gx')}
        on:pointermove={onDividerMove}
        on:pointerup={onDividerUp}
        on:pointercancel={onDividerUp}
      ></div>
    {/if}
  </div>

  <!-- A sibling of the scroller, not a child: a status bar that scrolls away
       with the notebook is not a status bar. Fixed to the bottom, with
       .focused-view inset above it by the same height. -->
  {#if $showStatusBar}
    <div class="status-dock"><StatusBar /></div>
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
          style="left:{nb.x}px;top:{nb.y}px;width:{nb.width}px;z-index:{$canvasState.activeId===nb.id?10:1};"
          on:pointerdown|capture={() => setActive(nb.id)}
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

    <!-- With a selection, the hint strip becomes the route into the split view.
         Rubber-band a few cards and open them together, which is the answer to
         "how do I get these two side by side" from the canvas. -->
    {#if $canvasState.selectedIds.length >= 2}
      <div class="canvas-hints selection-actions">
        <button class="hint-new-btn" on:click={() => openPanes($canvasState.selectedIds, 'h')}>
          Open {Math.min($canvasState.selectedIds.length, MAX_PANES)} side by side
        </button>
        <button class="hint-new-btn" on:click={() => openPanes($canvasState.selectedIds, 'v')}>
          Stacked
        </button>
        {#if $canvasState.selectedIds.length >= 3}
          <button class="hint-new-btn" on:click={() => openPanes($canvasState.selectedIds, 'grid')}>
            2×2
          </button>
        {/if}
        <span class="hint-sep">·</span>
        <span>⏎ or ⌘\</span>
        {#if $canvasState.selectedIds.length > MAX_PANES}
          <span class="hint-sep">·</span>
          <span>first {MAX_PANES} of {$canvasState.selectedIds.length}</span>
        {/if}
      </div>
    {:else}
      <div class="canvas-hints">
        <button class="hint-new-btn" on:click={addAtCentre}>＋ New Notebook</button>
        <span class="hint-sep">·</span>
        <span>drag-select 2+ to open side by side</span>
        <span class="hint-sep">·</span>
        <span>⌘\ split · ⌘0 fit · N key</span>
        <span class="hint-sep">·</span>
        <span>scroll pan · pinch zoom</span>
      </div>
    {/if}

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
    inset: var(--appbar-h, 34px) 0 0 0;   /* clear the app bar */
    width: 100vw;
    height: calc(100vh - var(--appbar-h, 34px));
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
  /* With a selection this strip is offering actions, not reciting hints, so it
     gets a surface and full opacity instead of the near-invisible grey. */
  .canvas-hints.selection-actions {
    pointer-events: auto;
    color: var(--text);
    background: var(--menu-bg);
    border: 1px solid var(--menu-border);
    border-radius: 9px;
    padding: 5px 9px;
    box-shadow: var(--menu-shadow);
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
    /* The toolbar, not the canvas bar: focused mode swaps the app bar for the
       taller grouped toolbar, and these two numbers must agree or the view
       either overlaps the bar or leaves a gap under it. */
    /* The toolbar, not the canvas bar: focused mode swaps the app bar for the taller grouped
       toolbar, and these two numbers must agree with App.svelte's .app-bar.toolbar-mode height
       or the view either overlaps the bar or leaves a gap under it. */
    inset: var(--toolbar-h, 46px) 0 0 0;

    /* Use card-bg so light mode doesn't show dark canvas edges */
    background: var(--card-bg, #050810);
    z-index: 50;
    /* A grid of panes. overflow:hidden is load-bearing, not tidiness: while the
       WINDOW was the scroller, dragging a divider scrolled the page and no pane
       could scroll on its own. Each .focused-pane is now its own scroller, which
       is what "independent panes" actually requires. */
    display: grid;
    overflow: hidden;
  }

  .focused-pane {
    /* Both non-optional. A grid item's automatic minimum size is min-content,
       and a notebook full of CodeMirror never reports a small one -- so without
       these the tracks blow past the container and the fr maths silently stops
       meaning anything. */
    min-width: 0;
    min-height: 0;
    overflow-y: auto;
    overflow-x: hidden;
    position: relative;
    display: flex;
    flex-direction: column;
    /* The canvas sets touch-action:none to own pan/zoom; re-enable finger scroll
       inside a pane. */
    touch-action: pan-y;
  }
  /* Seams only when tiled: with one pane there is no header, no border and no
     accent, so the single-notebook view is unchanged. */
  .focused-view.multi .focused-pane { border: 1px solid var(--border); }

  /* The active pane gets 2px on one edge and nothing else. The rubber-band
     selection language (outline + double box-shadow) belongs to a different
     concept and would shout here. */
  .focused-view.multi .focused-pane.active { border-color: var(--accent); }
  .focused-view.multi .focused-pane.active::before {
    content: '';
    position: absolute;
    inset: 0 0 auto 0;
    height: 2px;
    background: var(--accent);
    z-index: 11;
  }

  .pane-header {
    position: sticky;
    top: 0;
    z-index: 10;
    display: flex;
    align-items: center;
    gap: 6px;
    height: 24px;
    flex-shrink: 0;
    padding: 0 6px 0 9px;
    background: var(--gutter-bg);
    border-bottom: 1px solid var(--border);
    /* Deliberately slimmer than the card's own 36px titlebar, so a pane header
       never reads as a second app bar. */
    font: 600 11px/1 var(--sans);
    color: var(--text);
    user-select: none;
    -webkit-user-select: none;
    /* iOS: a top-row pane's header would otherwise sit under the notch. */
    padding-top: env(safe-area-inset-top, 0px);
  }
  /* Dim the CHROME of an inactive pane, never its content: a notebook you can
     still read and scroll must not look disabled. */
  .focused-pane:not(.active) .pane-header { opacity: 0.55; }

  .pane-title { flex: 1; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }

  .pane-close {
    background: none;
    border: none;
    color: var(--text-muted);
    font-size: 11px;
    line-height: 1;
    padding: 2px 4px;
    border-radius: 4px;
    cursor: pointer;
    flex-shrink: 0;
  }
  .pane-close:hover { background: var(--surface-2); color: var(--err); }

  /* ---- Dividers ---- */
  .pane-divider {
    grid-column: auto;
    cursor: col-resize;
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 12;
    /* The pane's pan-y would otherwise steal the drag on a touchscreen. */
    touch-action: none;
  }
  .pane-divider.vertical { cursor: row-resize; }

  .pane-divider::before {
    content: '';
    display: block;
    width: 2px;
    height: 46px;
    max-height: 70%;
    border-radius: 2px;
    background: var(--tb-rule);
    transition: background 0.15s, height 0.15s;
  }
  .pane-divider.vertical::before { width: 46px; max-width: 70%; height: 2px; }
  .pane-divider:hover::before, .pane-divider:active::before { background: var(--accent); height: 64px; }
  .pane-divider.vertical:hover::before, .pane-divider.vertical:active::before { width: 64px; height: 2px; }

  /* The two crossing dividers of the 2x2. The column one is painted above so it
     owns the 10x10 intersection and a drag there moves one axis, not both. */
  .pane-divider.grid-col { grid-column: 2; grid-row: 1 / -1; z-index: 14; }
  .pane-divider.grid-row { grid-column: 1 / -1; grid-row: 2; z-index: 13; }

  /* During a drag, panes stop taking pointer events so a fast movement across
     one cannot start a CodeMirror text selection behind the divider. */
  .focused-view.dragging .focused-pane { pointer-events: none; user-select: none; }

  /* Leave room for the status bar when it is shown, so the notebook's last cell
     is never hidden behind it. */
  .focused-view.with-status { bottom: var(--statusbar-h, 22px); }

  .status-dock {
    position: fixed;
    left: 0;
    right: 0;
    bottom: 0;
    z-index: 60;   /* above .focused-view (50), below the app bar (200) */
  }

  .focused-view-inner {
    width: 100%;
    /* Grow to fill, but never shrink below content. */
    flex: 1 0 auto;
    display: flex;
    flex-direction: column;
    /* No max-width, no side padding — notebook card fills the window */
  }
  /* Override card styles when in focused view so it has no border-radius or side margins */
  .focused-view-inner :global(.nb-card) {
    border-radius: 0;
    border-left: none;
    border-right: none;
    flex: 1 0 auto;
  }

</style>
