<!--
  Minimap.svelte — bird's-eye overview of the canvas.
  Shows all notebooks proportionally sized/positioned.
  Visible when zoom < 0.7. pointer-events auto so clicking works.
  Clicking a notebook rect pans+zooms the canvas to that notebook.
-->
<script lang="ts">
  import type { CanvasNotebook } from './canvas';

  export let notebooks: CanvasNotebook[] = [];
  export let panX: number = 0;
  export let panY: number = 0;
  export let zoom: number = 1.0;
  export let viewportW: number = 1280;
  export let viewportH: number = 800;
  // Callback when user clicks a notebook in the minimap
  export let onNotebookClick: ((nb: CanvasNotebook) => void) | null = null;

  const MAP_W = 280;
  const MAP_H = 200;
  const PAD   = 80;

  $: show = zoom < 0.7;

  // Compute world-space bounding box of all notebooks
  $: bbox = (() => {
    if (!notebooks.length) return { x: 0, y: 0, w: 1200, h: 800 };
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const nb of notebooks) {
      const h = nb.collapsed ? 52 : (nb.height ?? 400);
      minX = Math.min(minX, nb.x);
      minY = Math.min(minY, nb.y);
      maxX = Math.max(maxX, nb.x + nb.width);
      maxY = Math.max(maxY, nb.y + h);
    }
    return {
      x: minX - PAD,
      y: minY - PAD,
      w: Math.max(1, maxX - minX + PAD * 2),
      h: Math.max(1, maxY - minY + PAD * 2),
    };
  })();

  // Uniform scale to fit bounding box into minimap, preserving aspect ratio
  $: scale = Math.min(MAP_W / bbox.w, MAP_H / bbox.h);
  $: offsetX = (MAP_W - bbox.w * scale) / 2;
  $: offsetY = (MAP_H - bbox.h * scale) / 2;

  function toMapX(wx: number) { return offsetX + (wx - bbox.x) * scale; }
  function toMapY(wy: number) { return offsetY + (wy - bbox.y) * scale; }

  const PALETTE = ['#89b4fa','#a6e3a1','#f38ba8','#fab387','#cba6f7','#94e2d5',
                   '#89dceb','#f9e2af','#cba6f7','#b4befe'];

  $: nbRects = notebooks.map((nb, i) => {
    const h = nb.collapsed ? 52 : (nb.height ?? 400);
    return {
      nb,
      x: toMapX(nb.x),
      y: toMapY(nb.y),
      w: Math.max(6, nb.width * scale),
      h: Math.max(4, h * scale),
      color: PALETTE[i % PALETTE.length],
      label: nb.title.slice(0, 18),
    };
  });

  // Viewport rectangle in minimap coords
  $: vpRect = (() => {
    const wx = -panX / zoom;
    const wy = -panY / zoom;
    return {
      x: toMapX(wx),
      y: toMapY(wy),
      w: Math.max(8, (viewportW / zoom) * scale),
      h: Math.max(8, (viewportH / zoom) * scale),
    };
  })();

  // Click: find which notebook was clicked and call onNotebookClick
  function handleClick(e: MouseEvent) {
    e.stopPropagation();
    if (!onNotebookClick) return;
    const rect = (e.currentTarget as SVGSVGElement).getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    for (const r of nbRects) {
      if (mx >= r.x && mx <= r.x + r.w && my >= r.y && my <= r.y + r.h) {
        onNotebookClick(r.nb);
        return;
      }
    }
  }
</script>

{#if show}
  <div class="minimap-wrap">
    <!-- svelte-ignore a11y-click-events-have-key-events a11y-no-static-element-interactions -->
    <svg
      width={MAP_W}
      height={MAP_H}
      class="minimap-svg"
      on:click={handleClick}
    >
      <!-- Notebook rectangles -->
      {#each nbRects as r}
        <rect
          x={r.x} y={r.y}
          width={r.w} height={r.h}
          fill={r.color} fill-opacity="0.28"
          rx="2"
          class="nb-rect"
        />
        <rect
          x={r.x} y={r.y}
          width={r.w} height={r.h}
          fill="none" stroke={r.color} stroke-opacity="0.5" stroke-width="0.8"
          rx="2"
        />
        {#if r.h > 10}
          <text
            x={r.x + 3} y={r.y + Math.min(10, r.h - 2)}
            font-size="6.5" fill={r.color} fill-opacity="0.9"
            font-family="SF Mono, monospace"
            style="pointer-events:none"
          >{r.label}</text>
        {/if}
      {/each}

      <!-- Viewport indicator -->
      <rect
        x={vpRect.x} y={vpRect.y}
        width={vpRect.w} height={vpRect.h}
        fill="rgba(255,255,255,0.06)"
        stroke="rgba(255,255,255,0.55)"
        stroke-width="1.2"
        rx="1"
        style="pointer-events:none"
      />
    </svg>
    <div class="minimap-label">click to jump</div>
  </div>
{/if}

<style>
  .minimap-wrap {
    position: fixed;
    bottom: 42px;
    left: 14px;
    background: rgba(8, 10, 22, 0.88);
    border: 1px solid rgba(255,255,255,0.12);
    border-radius: 8px;
    overflow: hidden;
    z-index: 90;
    box-shadow: 0 4px 20px rgba(0,0,0,0.45);
    /* pointer-events: auto so clicks reach the SVG */
    pointer-events: auto;
    cursor: pointer;
    transition: box-shadow 0.15s;
  }
  .minimap-wrap:hover {
    box-shadow: 0 4px 28px rgba(0,0,0,0.6), 0 0 0 1px rgba(137,180,250,0.2);
  }

  .minimap-svg { display: block; }

  .nb-rect { cursor: pointer; }
  .nb-rect:hover { fill-opacity: 0.5 !important; }

  .minimap-label {
    text-align: right;
    padding: 1px 6px 3px;
    font-size: 0.52rem;
    color: rgba(255,255,255,0.2);
    letter-spacing: 0.04em;
    font-family: 'SF Mono', monospace;
    line-height: 1;
  }
</style>
