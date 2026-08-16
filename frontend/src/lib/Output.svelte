<script lang="ts">
  import type { OutputItem } from './notebook';
  import katex from 'katex';
  import 'katex/dist/katex.min.css';

  /* Opens a symbol's own reference page. Passed in so this stays a renderer. */
  export let onOpenDoc: ((name: string) => void) | null = null;
  export let items: OutputItem[] = [];

  // Max height before output is collapsed with a "Show more" toggle
  const MAX_HEIGHT = 180; // px

  // Per-item expanded/overflow state
  let expanded: Record<number, boolean> = {};
  let overflows: Record<number, boolean> = {};

  // Svelte action: measures actual scrollHeight vs offsetHeight.
  // Triggers reactivity only when overflow state changes.
  function measureOverflow(node: HTMLElement, idx: number) {
    function check() {
      // While expanded the clamp is lifted, so scrollHeight == offsetHeight and
      // the node reports "no overflow" -- which would hide the collapse toggle
      // and strand the user in the expanded state. Measure against the clamp
      // instead, so the answer means the same thing in both states.
      const does = expanded[idx]
        ? node.scrollHeight > MAX_HEIGHT + 4
        : node.scrollHeight > node.offsetHeight + 4;
      if (overflows[idx] !== does) {
        overflows[idx] = does;
        overflows = { ...overflows };
      }
    }
    // First check after mount, then again after async content (KaTeX, Plotly).
    requestAnimationFrame(check);
    const t = setTimeout(check, 150);

    // Re-measure on resize. Without this the overflow flag is decided once, at
    // the width the card happened to have on mount: widening the card can make
    // text that needed the clamp fit, and narrowing it can make text that fit
    // need the clamp, but neither re-ran the check -- so the toggle went stale
    // and stopped matching what was on screen.
    let ro: ResizeObserver | undefined;
    if (typeof ResizeObserver !== 'undefined') {
      ro = new ResizeObserver(() => check());
      ro.observe(node);
    }
    return {
      update() { check(); },
      destroy() { clearTimeout(t); ro?.disconnect(); }
    };
  }

  /* A usage message is structured, not a blob: alternating signature lines
     (flush left) and their descriptions (indented, hard-wrapped at ~70 columns
     for the terminal REPL). Rendering it as one <pre> loses all of that -- the
     hard wraps rewrap again at the card edge, a wrapped continuation drops back
     to column 0 so the indent stops meaning anything, and prose set in a
     monospace column is simply harder to read than it needs to be.

     Parse it back into blocks instead: each flush-left line is a signature,
     each run of indented lines is one description paragraph. The renderer then
     sets signatures in mono and descriptions as ordinary text, which is how
     Mathematica displays ?sym too. */
  type UsageBlock = { kind: 'sig' | 'body'; text: string };

  function parseUsage(text: string): UsageBlock[] {
    const out: UsageBlock[] = [];
    for (const line of text.split('\n')) {
      if (!line.trim()) continue;
      const indented = /^\s+\S/.test(line);
      const prev = out.length ? out[out.length - 1] : null;
      if (indented && prev && prev.kind === 'body') {
        prev.text += ' ' + line.trim();          // continuation of the paragraph
      } else if (indented) {
        out.push({ kind: 'body', text: line.trim() });
      } else {
        out.push({ kind: 'sig', text: line.trim() });
      }
    }
    return out;
  }

  // Heuristic: expressions that are long lists of numbers/symbols don't
  // benefit from KaTeX (no fractions/superscripts) and KaTeX can't wrap them.
  // Render as code with word-break so they don't overflow the card.
  // Heuristic: if output has >4 commas or is long, it's a list/sequence.
  // KaTeX can't wrap math spans so we use plain code with word-break.
  function isListOutput(text: string): boolean {
    const commas = (text.match(/,/g) ?? []).length;
    return commas > 4 || text.length > 200;
  }

  function renderOutput(text: string, latex?: string): string {
    // Long lists: always use wrapping code regardless of latex field.
    // KaTeX renders math spans without line-breaking, so even \{1,2,...\}
    // produces a single wide unbreakable line.
    if (isListOutput(text)) {
      const wrapped = text.replace(/,\s+/g, ', ');
      return `<code class="out-code-wrap">${wrapped}</code>`;
    }
    // Short expressions: prefer LaTeX from the kernel (StandardForm)
    if (latex && latex.length > 0) {
      try {
        return katex.renderToString(latex, { throwOnError: false, displayMode: false });
      } catch { /* fall through */ }
    }
    try {
      return katex.renderToString(text, { throwOnError: false, displayMode: false });
    } catch {
      return `<code>${text}</code>`;
    }
  }

  /* Draw a base64-RGBA payload onto a canvas.
   *
   * putImageData rather than an <img src="data:image/png"> because the kernel sends raw samples, not
   * an encoded file: encoding a PNG in C to have the browser decode it again would be two conversions
   * to display what is already a pixel buffer.
   *
   * The canvas is its natural pixel size and CSS scales it, with image-rendering: pixelated, so a 3x3
   * image is a visible 3x3 grid of squares rather than a 3-pixel dot or a blurred smear. */
  /* Per-output display width, in CSS pixels, once the reader has dragged the corner.
     Keyed by output index and deliberately NOT written back into the notebook: a display
     size is a way of looking at a result, not part of it, and persisting it would make an
     unedited notebook dirty. */
  let imgWidth: Record<number, number> = {};

  /* The default is a compromise the reader can override: a tiny image is useless at its
     natural size (an 8x8 result is a speck) and a large one must not push the cell wider
     than the pane, so small images are magnified and everything is capped. */
  function defaultImgWidth(it: { w: number }) {
    return Math.min(Math.max(it.w * (it.w < 64 ? 12 : 1), 48), 512);
  }

  function shownWidth(idx: number, it: { w: number }) {
    return imgWidth[idx] ?? defaultImgWidth(it);
  }

  /* Corner drag. Pointer events with setPointerCapture rather than mousemove on window:
     capture keeps the stream coming when the cursor leaves the handle (which it does
     immediately, since the handle moves with the image), and one release ends it. */
  function startResize(ev: PointerEvent, idx: number, it: { w: number }) {
    ev.preventDefault();
    ev.stopPropagation();
    const handle = ev.currentTarget as HTMLElement;
    const frame = handle.parentElement as HTMLElement | null;
    const startX = ev.clientX;
    const startW = shownWidth(idx, it);
    handle.setPointerCapture(ev.pointerId);

    /* THE DRAG DOES NOT TOUCH SVELTE STATE. Writing `imgWidth` per pointermove invalidated the
       component, which re-rendered the whole output list and re-laid-out the page around it --
       on a reference page of sixty examples that is a full document reflow per mouse move, and
       the resize visibly lagged the cursor. During the drag the width is written straight onto
       the one element that changes, and the state is committed once on release so the size
       survives a later re-render.

       Coalesced to one write per animation frame as well: a trackpad emits pointermove far faster
       than the compositor paints, and every extra write was a layout nobody saw. */
    let pending = 0;
    let latest = startW;
    const flush = () => {
      pending = 0;
      if (frame) frame.style.width = `${latest}px`;
    };

    const move = (e: PointerEvent) => {
      /* Clamped at both ends: below ~24px the handle would be most of the image and the
         drag could not be undone by dragging back. */
      latest = Math.max(24, Math.min(4096, startW + (e.clientX - startX)));
      if (!pending) pending = requestAnimationFrame(flush);
    };
    const up = (e: PointerEvent) => {
      handle.releasePointerCapture(e.pointerId);
      handle.removeEventListener('pointermove', move);
      handle.removeEventListener('pointerup', up);
      handle.removeEventListener('pointercancel', up);
      if (pending) { cancelAnimationFrame(pending); flush(); }
      imgWidth[idx] = latest;
      imgWidth = imgWidth;                    /* Svelte 4: reassign to trigger */
    };
    handle.addEventListener('pointermove', move);
    handle.addEventListener('pointerup', up);
    handle.addEventListener('pointercancel', up);
  }

  /* Double-click the handle to go back to the default, so a drag is never a one-way door. */
  function resetResize(idx: number, ev?: Event) {
    /* The drag writes an inline width directly onto the frame, so clearing the state is not
       enough -- the inline style would win over the re-rendered one. */
    const handle = ev?.currentTarget as HTMLElement | undefined;
    const frame = handle?.parentElement as HTMLElement | undefined;
    if (frame) frame.style.width = '';
    delete imgWidth[idx];
    imgWidth = imgWidth;
  }

  /* ---------------------------------------------------------------- volumes
   *
   * An Image3D drawn as one slice is a volume the reader cannot see -- the Image3D page showed a
   * 4x3 rectangle and nothing about it said "block". Mathematica draws a shaded box carrying the
   * voxels on its visible faces, and lets you turn it; this does the same.
   *
   * ORTHOGRAPHIC PROJECTION MAKES THIS EXACT, NOT APPROXIMATE. Under an orthographic camera the
   * projection of a rectangular face is an AFFINE map of that face's texture, so each face is one
   * `setTransform` + `drawImage` -- no triangles, no WebGL, no per-pixel work in JavaScript. The
   * six faces arrive from the kernel; the three facing the camera are drawn, back to front.
   *
   * Shading is per-face by the angle between its normal and the light, which is what keeps the
   * three visible faces distinguishable when they carry similar voxels: without it a cube of
   * uniform grey reads as a flat hexagon.
   */
  type Face = { w: number; h: number; data: string };
  type FaceMap = Record<string, Face>;

  /* Rotation per output index, in radians. Yaw around the vertical, pitch around the horizontal.
     The default is the three-quarter view every 3-D thumbnail uses: enough of the top and side to
     read as a box, with the front face still square-on enough to read as an image. */
  let volYaw: Record<number, number> = {};
  let volPitch: Record<number, number> = {};
  const YAW0 = -0.62, PITCH0 = 0.42;

  const yawOf = (i: number) => volYaw[i] ?? YAW0;
  const pitchOf = (i: number) => volPitch[i] ?? PITCH0;

  /* One decoded face, cached: decoding six base64 faces on every animation frame of a drag would
     be six atob() calls per frame for no gain, since the pixels never change. */
  const faceCache = new Map<string, HTMLCanvasElement>();

  function faceCanvas(key: string, f: Face): HTMLCanvasElement | null {
    const hit = faceCache.get(key);
    if (hit) return hit;
    const need = f.w * f.h * 4;
    const bin = atob(f.data);
    if (bin.length < need) return null;
    const buf = new Uint8ClampedArray(need);
    for (let i = 0; i < need; i++) buf[i] = bin.charCodeAt(i);
    const c = document.createElement('canvas');
    c.width = f.w; c.height = f.h;
    const cx = c.getContext('2d');
    if (!cx) return null;
    cx.putImageData(new ImageData(buf, f.w, f.h), 0, 0);
    faceCache.set(key, c);
    return c;
  }

  /* The unit cube's eight corners, scaled so the longest axis is 1 and the others keep their
     proportion -- a 4 x 4 x 32 volume has to look like a column, not a cube. */
  function cubeCorners(w: number, h: number, d: number) {
    const m = Math.max(w, h, d);
    const sx = w / m, sy = h / m, sz = d / m;
    const pts: [number, number, number][] = [];
    for (const z of [-sz / 2, sz / 2])
      for (const y of [-sy / 2, sy / 2])
        for (const x of [-sx / 2, sx / 2]) pts.push([x, y, z]);
    return pts;   /* index = x + 2*y + 4*z, each 0 or 1 */
  }

  function rotate(p: [number, number, number], yaw: number, pitch: number) {
    const [x, y, z] = p;
    /* Yaw about the y axis, then pitch about the x axis. */
    const cx = Math.cos(yaw), sx = Math.sin(yaw);
    const x1 = x * cx + z * sx, z1 = -x * sx + z * cx;
    const cy = Math.cos(pitch), sy = Math.sin(pitch);
    const y2 = y * cy - z1 * sy, z2 = y * sy + z1 * cy;
    return [x1, y2, z2] as [number, number, number];
  }

  /* Each face as the three corners the affine map needs (origin, +u, +v) plus its outward normal.
     Corner indices are into cubeCorners' x + 2y + 4z ordering, and the u/v directions match the
     kernel's face export: front/back are (x, y), left/right are (z, y), top/bottom are (x, z). */
  const FACE_GEOM: Record<string, { o: number; u: number; v: number; n: [number, number, number] }> = {
    front:  { o: 0, u: 1, v: 2, n: [0, 0, -1] },
    back:   { o: 5, u: 4, v: 7, n: [0, 0, 1] },
    left:   { o: 4, u: 0, v: 6, n: [-1, 0, 0] },
    right:  { o: 1, u: 5, v: 3, n: [1, 0, 0] },
    top:    { o: 0, u: 1, v: 4, n: [0, -1, 0] },
    bottom: { o: 2, u: 3, v: 6, n: [0, 1, 0] },
  };

  function mountVolume(node: HTMLCanvasElement,
                       arg: { item: { w: number; h: number; depth?: number; faces?: FaceMap }, idx: number }) {
    let cur = arg;

    const draw = (a: typeof arg) => {
      cur = a;
      const { item, idx } = a;
      const faces = item.faces;
      const ctx = node.getContext('2d');
      if (!ctx || !faces) return;

      /* A fixed backing size with CSS scaling on top: the drag rotates the box, the corner handle
         resizes the element, and the two must not fight over the pixel buffer. */
      const S = 320;
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      node.width = S * dpr; node.height = S * dpr;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, S, S);
      ctx.imageSmoothingEnabled = false;

      const yaw = yawOf(idx), pitch = pitchOf(idx);
      const pts = cubeCorners(item.w, item.h, item.depth ?? 1).map(p => rotate(p, yaw, pitch));
      const scale = S * 0.62;
      const proj = (p: [number, number, number]) => [S / 2 + p[0] * scale, S / 2 + p[1] * scale];

      /* Light from the upper left, in camera space. */
      const L = [-0.45, -0.75, -0.5];
      const order = Object.keys(FACE_GEOM)
        .map(k => {
          const g = FACE_GEOM[k];
          const n = rotate(g.n, yaw, pitch);
          /* Depth of the face centre, for back-to-front ordering. */
          const c = [(pts[g.o][0] + pts[g.u][0] + pts[g.v][0]) / 3,
                     (pts[g.o][1] + pts[g.u][1] + pts[g.v][1]) / 3,
                     (pts[g.o][2] + pts[g.u][2] + pts[g.v][2]) / 3];
          return { k, g, n, z: c[2], visible: n[2] < 0 };
        })
        .filter(f => f.visible && faces[f.k])
        .sort((a, b) => b.z - a.z);

      for (const f of order) {
        const cvs = faceCanvas(`${idx}:${f.k}`, faces[f.k]);
        if (!cvs) continue;
        const [ox, oy] = proj(pts[f.g.o]);
        const [ux, uy] = proj(pts[f.g.u]);
        const [vx, vy] = proj(pts[f.g.v]);
        /* The affine map taking the texture's (0,0)-(w,0)-(0,h) to those three points. */
        const a = (ux - ox) / cvs.width,  b = (uy - oy) / cvs.width;
        const c = (vx - ox) / cvs.height, d = (vy - oy) / cvs.height;
        ctx.save();
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        ctx.transform(a, b, c, d, ox, oy);
        ctx.drawImage(cvs, 0, 0);
        ctx.restore();

        /* Per-face shading, so three faces of similar voxels still read as three faces. */
        const nl = Math.abs(f.n[0] * L[0] + f.n[1] * L[1] + f.n[2] * L[2]);
        const shade = 0.38 * (1 - Math.min(1, nl));
        if (shade > 0.01) {
          ctx.save();
          ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
          ctx.beginPath();
          const p4 = [pts[f.g.o], pts[f.g.u],
                      /* the fourth corner: o + (u - o) + (v - o) */
                      [pts[f.g.u][0] + pts[f.g.v][0] - pts[f.g.o][0],
                       pts[f.g.u][1] + pts[f.g.v][1] - pts[f.g.o][1],
                       pts[f.g.u][2] + pts[f.g.v][2] - pts[f.g.o][2]] as [number, number, number],
                      pts[f.g.v]];
          p4.forEach((q, i2) => {
            const [X, Y] = proj(q as [number, number, number]);
            if (i2 === 0) ctx.moveTo(X, Y); else ctx.lineTo(X, Y);
          });
          ctx.closePath();
          ctx.fillStyle = `rgba(0,0,0,${shade})`;
          ctx.fill();
          ctx.restore();
        }
      }
    };

    draw(arg);
    /* Redraw on rotation without remounting: the action's update only fires when `arg` changes, so
       the drag handler calls this directly. */
    (node as any).__redraw = () => draw(cur);
    return { update: draw };
  }

  /* Drag to rotate. Same pointer-capture shape as the resize handle, and the same reason: the
     cursor leaves the element immediately once the box turns under it. */
  function startVolDrag(ev: PointerEvent, idx: number) {
    ev.preventDefault();
    const node = ev.currentTarget as HTMLCanvasElement;
    const x0 = ev.clientX, y0 = ev.clientY;
    const yaw0 = yawOf(idx), pitch0 = pitchOf(idx);
    node.setPointerCapture(ev.pointerId);

    /* One redraw per animation frame, for the reason the resize handle documents: pointermove
       outruns the compositor, and three drawImage calls per event that never gets painted is
       work thrown away. Rotation is local to this canvas, so no Svelte state is touched until
       the pointer is released. */
    let pending = 0;
    let lastEv: PointerEvent | null = null;
    const flush = () => {
      pending = 0;
      const e = lastEv;
      if (!e) return;
      volYaw[idx] = yaw0 + (e.clientX - x0) * 0.011;
      /* Pitch is clamped just short of a pole: at exactly +-pi/2 the box degenerates to a line and
         there is no way to tell which way a further drag should go. */
      const lim = Math.PI / 2 - 0.05;
      volPitch[idx] = Math.max(-lim, Math.min(lim, pitch0 + (e.clientY - y0) * 0.011));
      (node as any).__redraw?.();
    };

    const move = (e: PointerEvent) => {
      lastEv = e;
      if (!pending) pending = requestAnimationFrame(flush);
    };
    const up = (e: PointerEvent) => {
      node.releasePointerCapture(e.pointerId);
      node.removeEventListener('pointermove', move);
      node.removeEventListener('pointerup', up);
      node.removeEventListener('pointercancel', up);
      if (pending) { cancelAnimationFrame(pending); flush(); }
    };
    node.addEventListener('pointermove', move);
    node.addEventListener('pointerup', up);
    node.addEventListener('pointercancel', up);
  }

  function resetVol(idx: number) {
    delete volYaw[idx];
    delete volPitch[idx];
    volYaw = volYaw; volPitch = volPitch;
  }

  function mountImage(node: HTMLCanvasElement,
                      item: { w: number; h: number; data: string }) {
    const draw = (it: { w: number; h: number; data: string }) => {
      if (!it.w || !it.h || !it.data) return;
      const ctx = node.getContext('2d');
      if (!ctx) return;
      node.width = it.w;
      node.height = it.h;
      const bin = atob(it.data);
      const need = it.w * it.h * 4;
      if (bin.length < need) return;          /* truncated payload: draw nothing rather than garbage */
      const buf = new Uint8ClampedArray(need);
      for (let i = 0; i < need; i++) buf[i] = bin.charCodeAt(i);
      ctx.putImageData(new ImageData(buf, it.w, it.h), 0, 0);
    };
    draw(item);
    return { update: draw };
  }

  function mountPlot(node: HTMLElement, data: object) {
    import('plotly.js-dist-min').then((Plotly: any) => {
      const spec = data as any;
      const dark = !document.documentElement.classList.contains('light');
      const sceneAxisDark  = { gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' }, backgroundcolor: '#181825', showbackground: true };
      const sceneAxisLight = { gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' }, backgroundcolor: '#f0f0f8', showbackground: true };
      const layoutOverride = dark ? {
        plot_bgcolor:  '#181825', paper_bgcolor: '#181825',
        font: { color: '#cdd6f4' },
        xaxis: { ...(spec.layout?.xaxis ?? {}), gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' } },
        yaxis: { ...(spec.layout?.yaxis ?? {}), gridcolor: '#313244', zerolinecolor: '#585b70', tickfont: { color: '#cdd6f4' } },
        scene: { ...(spec.layout?.scene ?? {}), bgcolor: '#181825', xaxis: sceneAxisDark, yaxis: sceneAxisDark, zaxis: sceneAxisDark },
      } : {
        plot_bgcolor:  '#ffffff', paper_bgcolor: '#f5f5fa',
        font: { color: '#1c1c2e' },
        xaxis: { ...(spec.layout?.xaxis ?? {}), gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' } },
        yaxis: { ...(spec.layout?.yaxis ?? {}), gridcolor: '#d8d9e8', zerolinecolor: '#9999bb', tickfont: { color: '#1c1c2e' } },
        scene: { ...(spec.layout?.scene ?? {}), bgcolor: '#f5f5fa', xaxis: sceneAxisLight, yaxis: sceneAxisLight, zaxis: sceneAxisLight },
      };
      Plotly.react(node, spec.data ?? [spec], { ...(spec.layout ?? {}), ...layoutOverride }, {
        responsive: true, displayModeBar: true,
      });
    });
  }

  // Measure height of a rendered output element to decide if it needs collapse
  function checkOverflow(node: HTMLElement, idx: number) {
    requestAnimationFrame(() => {
      if (node.scrollHeight > MAX_HEIGHT + 20) {
        // tall enough to warrant collapsing by default — nothing to do,
        // the CSS max-height handles it; the button appears via CSS
      }
    });
    return {};
  }
</script>

<div class="output">
  {#each items as item, idx (idx)}
    <div class="out-item" class:expanded={expanded[idx]} class:overflowing={overflows[idx]}>
      {#if item.kind === 'expr'}
        <div class="out-collapsible" use:measureOverflow={idx}>
          <div class="out-expr">{@html renderOutput(item.text, item.latex)}</div>
        </div>
      {:else if item.kind === 'expected'}
        <!-- A reference-page example that has not been run yet. Shown as plain
             text, never typeset: the recorded value is Mathilda syntax, and
             handing it to KaTeX (which is what an `expr` with no kernel LaTeX
             does) renders Derivative[1][g][x] as italic mathematics. -->
        <div class="out-collapsible" use:measureOverflow={idx}>
          <pre class="out-expected">{item.text}</pre>
        </div>
      {:else if item.kind === 'usage'}
        <!-- A usage message is documentation, not an expression: render it
             verbatim. It must not reach KaTeX, which cannot typeset it and
             would fall back to showing the InputForm escapes. -->
        <div class="out-collapsible" use:measureOverflow={idx}>
          <div class="out-usage">
            {#each parseUsage(item.text) as blk, bi (bi)}
              {#if blk.kind === 'sig'}
                <div class="usage-sig">{blk.text}</div>
              {:else}
                <p class="usage-body">{blk.text}</p>
              {/if}
            {/each}
          </div>
        </div>
          {#if item.symbol}
            <!-- The docstring answers "what does this do"; the page answers "show me it working".
                 Linked to Mathilda's OWN documentation, never an external site. -->
            <button class="usage-doc" on:click={() => { const sy = item.symbol; if (onOpenDoc && sy) onOpenDoc(sy); }}>
              Documentation: {item.symbol} &rarr;
            </button>
          {/if}
      {:else if item.kind === 'names'}
        <!-- `?pat*` is a symbol search. A grid reads far better than one long
             braced line, and each name is a discrete thing to scan for. -->
        <div class="out-collapsible" use:measureOverflow={idx}>
          {#if item.names.length === 0}
            <div class="names-empty">No symbol matches that pattern.</div>
          {:else}
            <div class="names-grid">
              {#each item.names as nm (nm)}
                <button
                  class="name-chip"
                  title={`Open the reference page for ${nm}`}
                  on:click={(e) => e.currentTarget.dispatchEvent(
                    new CustomEvent('mathilda-refpage',
                      { detail: { name: nm }, bubbles: true }))}
                >{nm}</button>
              {/each}
            </div>
          {/if}
        </div>
      {:else if item.kind === 'error'}
        <div class="out-error">{item.text}</div>
      {:else if item.kind === 'stream'}
        <div class="out-collapsible" use:measureOverflow={idx}>
          <pre class="out-stream">{item.text}</pre>
        </div>
      {:else if item.kind === 'image'}
        <div class="out-image">
          <!-- svelte-ignore a11y-no-static-element-interactions -->
          <div class="img-frame" style={`width: ${shownWidth(idx, item)}px`}>
            {#if item.faces && item.depth}
              <!-- A VOLUME, drawn as the block it is. Three of its six faces face the camera at
                   any angle, so three are drawn, back to front, each an affine map of its own
                   texture -- exact under an orthographic camera. Drag to turn it. -->
              <canvas
                class="out-canvas vol-canvas"
                title="Drag to rotate · double-click to reset the view"
                use:mountVolume={{ item, idx }}
                on:pointerdown={(e) => startVolDrag(e, idx)}
                on:dblclick={() => resetVol(idx)}
              ></canvas>
            {:else}
              <canvas class="out-canvas" use:mountImage={item}></canvas>
            {/if}
            <!-- The grab corner. Aspect ratio is preserved because only the width is set and
                 the canvas keeps `height: auto`, so an image cannot be squashed by accident. -->
            <div
              class="img-handle"
              title="Drag to resize · double-click to reset"
              on:pointerdown={(e) => startResize(e, idx, item)}
              on:dblclick={(e) => resetResize(idx, e)}
            ></div>
          </div>
          <span class="out-image-note">
            {item.w}&times;{item.h}{#if item.depth}&times;{item.depth}{/if}{item.channels > 1 ? `\u00d7${item.channels}` : ''}{#if item.depth && item.faces}
              &nbsp;· volume{:else if item.depth}
              &nbsp;· slice {item.slice} of {item.depth}{/if}
            {#if imgWidth[idx]}&nbsp;· shown at {Math.round(imgWidth[idx])}px{/if}
          </span>
        </div>
      {:else if item.kind === 'plot'}
        <div class="out-plot" use:mountPlot={item.data}></div>
      {:else if item.kind === 'html'}
        <div class="out-html">{@html item.html}</div>
      {/if}

      <!-- Always show toggle for collapsible output so user can expand/collapse -->
      {#if item.kind !== 'plot' && item.kind !== 'error' && item.kind !== 'image'}
        <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
        <div class="out-toggle" class:hidden={!overflows[idx] && !expanded[idx]} on:click={() => expanded[idx] = !expanded[idx]}>
          {expanded[idx] ? '▲ collapse' : '▼ show all'}
        </div>
      {/if}
    </div>
  {/each}
</div>

<style>
  .output {
    padding: 0.3rem 0.75rem 0.5rem;
    min-height: 1px;
    text-align: left;
    min-width: 0;    /* prevent output from pushing cell-content wider */
    overflow: hidden; /* clip anything that escapes a collapsible */
  }

  .out-item {
    position: relative;
    margin-bottom: 0.2rem;
    min-width: 0;
  }

  /* Collapsible wrapper: clips vertically, scrolls horizontally */
  .out-collapsible {
    width: 100%;
    min-width: 0;
    max-height: 180px;
    overflow-x: auto;
    overflow-y: hidden;
    /* Fade applied only when content actually overflows (via .overflowing class) */
  }

  /* Only fade when content genuinely overflows the cap */
  .overflowing .out-collapsible {
    -webkit-mask-image: linear-gradient(to bottom, black 55%, transparent 100%);
    mask-image:         linear-gradient(to bottom, black 55%, transparent 100%);
  }

  .expanded .out-collapsible {
    max-height: none;
    overflow-y: visible;
    -webkit-mask-image: none;
    mask-image: none;
  }

  /* Show-more button: hidden by default, shown only when content overflows */
  .out-toggle {
    font-size: 0.68rem;
    color: var(--accent, #89b4fa);
    cursor: pointer;
    padding: 2px 0 0;
    user-select: none;
    transition: opacity 0.1s;
    text-align: left;
  }
  .out-toggle:hover { opacity: 0.7; }
  /* Hide when content fits and not yet expanded */
  .out-toggle.hidden { display: none; }

  /* Expression output — overflow handled by parent .out-collapsible */
  .out-expr {
    font-size: 1.05em;
    padding: 0.25rem 0;
    color: var(--out-text, #222);
    text-align: left;
  }

  /* Long list/sequence outputs rendered as wrapping code */
  :global(.out-code-wrap) {
    display: block;
    font-family: 'SF Mono', 'Fira Code', monospace;
    font-size: 0.95em;
    color: var(--out-text, #cdd6f4);
    background: transparent;  /* override browser default <code> background */
    white-space: normal;
    word-break: break-word;
    overflow-wrap: break-word;
    line-height: 1.7;
    padding: 0.15rem 0;
  }

  /* Error output */
  .out-error {
    color: #e74c3c;
    font-family: 'SF Mono', monospace;
    font-size: 0.88em;
    background: rgba(231,76,60,0.08);
    border-left: 3px solid #e74c3c;
    padding: 0.4rem 0.8rem;
    border-radius: 3px;
    text-align: left;
    overflow-x: auto;
  }

  /* Stream (print) output */
  .out-stream {
    color: var(--text-muted);
    font-size: 0.84em;
    margin: 0;
    white-space: pre-wrap;
    word-break: break-all;
    font-family: 'SF Mono', 'Fira Code', monospace;
    text-align: left;
    overflow-x: hidden;
  }

  /* Usage message from `?sym`, rendered as structure rather than a blob:
     signatures in mono so they read as code, descriptions as ordinary text so
     a paragraph reads as a paragraph. The body is indented as a block, so a
     wrapped line stays aligned under the first instead of falling back to
     column 0 the way a <pre> would. */
  .out-usage {
    text-align: left;
    margin: 0.1rem 0 0.2rem;
  }

  .usage-sig {
    font-family: 'SF Mono', 'Fira Code', monospace;
    font-size: 0.84em;
    color: var(--text);
    margin-top: 0.55rem;
    word-break: break-word;
  }
  .usage-sig:first-child { margin-top: 0; }

  .usage-body {
    margin: 0.15rem 0 0 1.6em;
    font-size: 0.88em;
    line-height: 1.5;
    color: var(--text-muted);
    max-width: 78ch;          /* prose stops being readable much past this */
    word-break: break-word;
  }

  /* `?pat*` symbol search: an auto-fitting grid, so a wide card shows more
     columns instead of one tall column. minmax keeps a long name like
     NeighborhoodContraction from being clipped while still packing short
     ones. */
  .names-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(15rem, 1fr));
    gap: 0.15rem 0.9rem;
    text-align: left;
  }

  /* A <button>, not a <code>: the global `code` rule paints a themed
     background and padding, which turned each name into a filled box whose
     text stopped contrasting with it.

     Colour is `inherit`, deliberately, not a theme token. app.css defines its
     dark palette only under @media (prefers-color-scheme: dark) and has no
     .light override, but the app also has its own light/dark toggle that sets
     a class on <html>. Toggle the app to light while the OS is dark and every
     var() still holds its DARK value on a light surface -- which is why
     --accent read as washed out here and --text-h read as white-on-white.
     Inheriting from the surrounding output text sidesteps that entirely: these
     names are exactly as legible as the text beside them, in every combination
     of OS and app theme. The underline, not colour, carries the affordance. */
  .out-expected {
    font: 0.82rem/1.6 var(--mono);
    color: var(--text-dim, var(--text));
    margin: 0;
    padding: 0;
    white-space: pre-wrap;
    word-break: break-word;
    opacity: 0.85;
  }

  .names-grid .name-chip {
    font-family: 'SF Mono', 'Fira Code', monospace;
    font-size: 0.84em;
    color: inherit;
    background: none;
    border: none;
    padding: 0.1rem 0;
    margin: 0;
    text-align: left;
    cursor: pointer;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    text-decoration: underline;
    text-decoration-color: color-mix(in srgb, currentColor 40%, transparent);
    text-underline-offset: 2px;
  }
  .names-grid .name-chip:hover,
  .names-grid .name-chip:focus-visible {
    text-decoration-color: currentColor;
    text-decoration-thickness: 2px;
  }

  .names-empty {
    color: var(--text-muted);
    font-size: 0.86em;
    text-align: left;
  }

  /* Plot output */
  .out-plot {
    width: 100%;
    min-height: 320px;
  }

  /* HTML output */
  .out-html {
    font-size: 0.95em;
    text-align: left;
    color: var(--text);
    overflow-x: auto;
  }
  /* A recorded image result from a documentation page: nearest-neighbour like the live
     canvas, so a small image reads as a grid of squares rather than a blurred smear, and
     capped so a large one cannot push the cell wider than its pane. */
  :global(.out-html img.ref-shot) {
    image-rendering: pixelated;
    max-width: min(100%, 512px);
    height: auto;
    display: block;
  }
  .out-image { display: flex; flex-direction: column; gap: 4px; align-items: flex-start; }
  .img-frame {
    position: relative;
    display: inline-block;
    line-height: 0;          /* no descender gap under the canvas */
    max-width: 100%;
  }
  .img-frame .out-canvas { width: 100%; display: block; }
  /* A volume is square whatever its voxel dimensions -- the box is drawn inside a square field so
     that rotating it never changes the element's size and never reflows the cell. */
  .vol-canvas { aspect-ratio: 1 / 1; cursor: grab; touch-action: none; }
  .vol-canvas:active { cursor: grabbing; }
  .img-handle {
    position: absolute;
    right: -3px;
    bottom: -3px;
    width: 14px;
    height: 14px;
    cursor: nwse-resize;
    /* Two hairlines in the corner: the conventional resize affordance, and quiet enough
       not to compete with the pixels it sits on. */
    background:
      linear-gradient(135deg, transparent 0 55%, var(--text-muted) 55% 65%, transparent 65% 78%,
                      var(--text-muted) 78% 88%, transparent 88%);
    opacity: 0;
    transition: opacity 120ms ease;
    touch-action: none;      /* so a touch drag resizes instead of scrolling */
  }
  .img-frame:hover .img-handle { opacity: 0.9; }
  .out-canvas {
    /* Nearest-neighbour, so small images read as grids of squares rather than blurred smears. */
    image-rendering: pixelated;
    height: auto;
    max-width: 100%;
    border: 1px solid var(--border);
    background: var(--cell-bg);
  }
  .out-image-note { font-size: 10px; color: var(--text-muted); font-variant-numeric: tabular-nums; }
  .usage-doc {
    align-self: flex-start;
    margin-top: 4px;
    padding: 2px 6px;
    font: inherit;
    font-size: 11px;
    color: var(--accent);
    background: none;
    border: 1px solid var(--border);
    border-radius: 3px;
    cursor: pointer;
  }
  .usage-doc:hover { border-color: var(--accent); }
</style>
