<!--
  Minimap.svelte — disabled (was intercepting canvas clicks and had wrong math)
  Re-enable once pointer-events and coordinate math are fixed.
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
  const CARD_H_APPROX = 400; // rough height estimate when expanded

  // Bounding box of all notebooks in world space
  $: bounds = (() => {
    if (!notebooks.length) return { minX: 0, minY: 0, maxX: 1000, maxY: 600 };
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const nb of notebooks) {
      const h = nb.collapsed ? 52 : CARD_H_APPROX;
      minX = Math.min(minX, nb.x);
      minY = Math.min(minY, nb.y);
      maxX = Math.max(maxX, nb.x + nb.width);
      maxY = Math.max(maxY, nb.y + h);
    }
    // Add some padding
    const pad = 100;
    return { minX: minX - pad, minY: minY - pad, maxX: maxX + pad, maxY: maxY + pad };
  })();

  $: bbW = bounds.maxX - bounds.minX;
  $: bbH = bounds.maxY - bounds.minY;

  // Scale factor: world → minimap pixels
  $: scale = Math.min(MAP_W / bbW, MAP_H / bbH);

  function worldToMap(wx: number, wy: number): { x: number; y: number } {
    return {
      x: (wx - bounds.minX) * scale,
      y: (wy - bounds.minY) * scale,
    };
  }

  // Viewport rect in world space
  $: vpRect = (() => {
    const vpW = viewportW / zoom;
    const vpH = viewportH / zoom;
    const vpX = -panX / zoom;
    const vpY = -panY / zoom;
    const a = worldToMap(vpX, vpY);
    return {
      x: a.x,
      y: a.y,
      w: vpW * scale,
      h: vpH * scale,
    };
  })();

  const PALETTE = ['#89b4fa','#a6e3a1','#f38ba8','#fab387','#cba6f7','#94e2d5'];
  function nbColor(nb: CanvasNotebook) {
    return PALETTE[parseInt(nb.id.replace('nb-', ''), 10) % PALETTE.length] ?? '#89b4fa';
  }
</script>

{#if zoom < 0.7}
  <div class="minimap">
    <svg width={MAP_W} height={MAP_H}>
      <!-- Notebook rects -->
      {#each notebooks as nb}
        {@const pos = worldToMap(nb.x, nb.y)}
        {@const w = nb.width * scale}
        {@const h = (nb.collapsed ? 52 : CARD_H_APPROX) * scale}
        <rect
          x={pos.x}
          y={pos.y}
          width={Math.max(2, w)}
          height={Math.max(2, h)}
          fill={nbColor(nb)}
          fill-opacity="0.45"
          rx="1"
        />
      {/each}

      <!-- Viewport rect -->
      <rect
        x={vpRect.x}
        y={vpRect.y}
        width={Math.max(4, vpRect.w)}
        height={Math.max(4, vpRect.h)}
        fill="none"
        stroke="rgba(255,255,255,0.6)"
        stroke-width="1"
        rx="1"
      />
    </svg>
  </div>
{/if}

<style>
  .minimap {
    position: fixed;
    left: 12px;
    bottom: 40px;
    width: 160px;
    height: 110px;
    background: rgba(8,10,20,0.85);
    border: 1px solid rgba(255,255,255,0.1);
    border-radius: 6px;
    overflow: hidden;
    pointer-events: none;
    z-index: 90;
  }
</style>
