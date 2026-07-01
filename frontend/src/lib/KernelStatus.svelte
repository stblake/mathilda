<script lang="ts">
  import { kernelStatus } from './notebook';

  const labels: Record<string, string> = {
    starting:    '⬤ Starting',
    ready:       '⬤ Ready',
    busy:        '⬤ Busy',
    restarting:  '⬤ Restarting',
    dead:        '⬤ Dead',
  };
</script>

<span class="kernel-status" data-status={$kernelStatus}>
  {labels[$kernelStatus] ?? $kernelStatus}
</span>

<style>
  .kernel-status {
    font-size: 0.78rem;
    font-family: 'SF Mono', monospace;
    padding: 2px 10px;
    border-radius: 12px;
    letter-spacing: 0.02em;
    background: #f0f0f0;
    color: #555;
  }
  [data-status='ready']      { background: #e6f9f0; color: #27ae60; }
  [data-status='busy']       { background: #fff8e6; color: #e67e22; }
  [data-status='starting']   { background: #eaf3ff; color: #2980b9; }
  [data-status='restarting'] { background: #fff3cd; color: #856404; }
  [data-status='dead']       { background: #fdf0f0; color: #c0392b; }
</style>
