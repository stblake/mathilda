<!--
  Minimap.svelte
  Bird's-eye overview of the canvas. Shows notebook positions and the
  current viewport. Appears when zoom < 0.7. pointer-events: none so it
  never blocks canvas interactions.
-->
<script lang="ts">
  import type { CanvasNotebook } from './canvas';

  export let notebooks: CanvasNotebook[] = [];
  export let panX: number = 0;
  export let panY: number = 0;
  export let zoom: number = 1.0;
  export let viewportW: number = 1280;
  export let viewportH: number = 800;

  const MAP_W = 160;
  const MAP_H = 110;
  const PAD   = 60;   // world-space padding around bounding box

  const PALETTE = ['#89b4fa','#a6e3a1','#f38ba8','#fab387','#cba6f7','#94e2d5'];

  $: show = zoom < 0.7;

  // Compute world-space bounding box of all notebooks
  $: bbox = (() => {
    if (!notebooks.length) return { x: 0, y: 0, w: 1000, h: 700 };
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

  // Scale factors: world → minimap pixels
  $: scaleX = MAP_W / bbox.w;
  $: scaleY = MAP_H / bbox.h;

  function toMapX(wx: number) { return (wx - bbox.x) * scaleX; }
  function toMapY(wy: number) { return (wy - bbox.y) * scaleY; }

  // Notebook rects in minimap space
  $: nbRects = notebooks.map((nb, i) => ({
    x: toMapX(nb.x),
    y: toMapY(nb.y),
    w: nb.width * scaleX,
    h: (nb.collapsed ? 52 : (nb.height ?? 400)) * scaleY,
    color: PALETTE[i % PALETTE.length],
    title: nb.title,
  }));

  // Viewport rect in minimap space
  // The canvas-world is transformed by translate(panX, panY) scale(zoom).
  // The top-left world coordinate visible = -panX/zoom, -panY/zoom
  // The visible world size = viewportW/zoom × viewportH/zoom
  $: vpRect = (() => {
    const wx = -panX / zoom;
    const wy = -panY / zoom;
    const ww = viewportW / zoom;
    const wh = viewportH / zoom;
    return {
      x: toMapX(wx),
      y: toMapY(wy),
      w: ww * scaleX,
      h: wh * scaleY,
    };
  })();
</script>

{#if show}
  <div class="minimap">
    <svg width={MAP_W} height={MAP_H} xmlns="http://www.w3.org/2000/svg">
      <!-- Notebook rectangles -->
      {#each nbRects as r}
        <rect
          x={Math.max(0, r.x)}
          y={Math.max(0, r.y)}
          width={Math.min(MAP_W, r.w)}
          height={Math.min(MAP_H, r.h)}
          fill={r.color}
          fill-opacity="0.35"
          rx="2"
        />
        <!-- Notebook title label -->
        <text
          x={Math.max(2, r.x + 3)}
          y={Math.max(8, r.y + 9)}
          font-size="6"
          fill={r.color}
          fill-opacity="0.9"
          font-family="SF Mono, monospace"
        >{r.title.slice(0, 14)}</text>
      {/each}

      <!-- Viewport indicator -->
      <rect
        x={vpRect.x}
        y={vpRect.y}
        width={Math.max(4, vpRect.w)}
        height={Math.max(4, vpRect.h)}
        fill="none"
        stroke="rgba(255,255,255,0.5)"
        stroke-width="1"
        rx="1"
      />
    </svg>
    <div class="minimap-label">overview</div>
  </div>
{/if}

<style>
  .minimap {
    position: fixed;
    bottom: 42px;
    left: 14px;
    width: 160px;
    height: 110px;
    background: rgba(8, 10, 20, 0.82);
    border: 1px solid rgba(255,255,255,0.12);
    border-radius: 7px;
    overflow: hidden;
    pointer-events: none;   /* CRITICAL: never block canvas clicks */
    z-index: 90;
    box-shadow: 0 4px 16px rgba(0,0,0,0.4);
  }
  .minimap-label {
    position: absolute;
    bottom: 3px;
    right: 6px;
    font-size: 0.55rem;
    color: rgba(255,255,255,0.25);
    letter-spacing: 0.05em;
    font-family: 'SF Mono', monospace;
  }
  svg {
    display: block;
  }
</style>
