"""Decorative background: small ball-and-stick molecules drifting behind the page.

They drift on their own, are pushed away by the mouse pointer, and faint bonds reach from
the pointer to nearby molecules. It is one canvas under all content, drawn by a short script
that Streamlit runs in the page itself (``unsafe_allow_javascript``). The script is a fixed
string with no user input, which is what makes that flag acceptable here.
"""

from __future__ import annotations

import streamlit as st

# Idempotent: Streamlit re-sends the element on every rerun, but only one canvas and one
# animation loop may exist. Honours "reduce motion" by drawing a single still frame.
SCRIPT = """
<script>
(() => {
  if (window.__moleculeBackground) return;
  const app = document.querySelector('[data-testid="stApp"]');
  if (!app) return;

  const canvas = document.createElement("canvas");
  canvas.id = "molecule-background";
  canvas.setAttribute("aria-hidden", "true");
  app.prepend(canvas);
  window.__moleculeBackground = canvas;
  const ctx = canvas.getContext("2d");

  const reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  const COLOURS = ["108,92,231", "166,108,255", "255,126,182", "56,189,248"];
  const TEMPLATES = [
    {atoms: [[0, 0, 7], [-15, 11, 4.5], [15, 11, 4.5]], bonds: [[0, 1], [0, 2]]},  // water
    {atoms: [[-18, 0, 5], [0, 0, 6], [18, 0, 5]], bonds: [[0, 1], [1, 2]]},  // CO2
    {atoms: Array.from({length: 6}, (_, i) => [
      17 * Math.cos(i * Math.PI / 3), 17 * Math.sin(i * Math.PI / 3), 4.5]),
     bonds: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]},  // ring
    {atoms: [[-27, 6, 4.5], [-9, -6, 4.5], [9, 6, 4.5], [27, -6, 4.5]],
     bonds: [[0, 1], [1, 2], [2, 3]]},  // chain
    {atoms: [[0, 0, 6.5], [-17, -11, 4.5], [17, -11, 4.5], [0, 20, 5], [26, 18, 4]],
     bonds: [[0, 1], [0, 2], [0, 3], [3, 4]]},  // branched
  ];
  const PUSH_RADIUS = 170;
  const LINK_RADIUS = 230;

  const mouse = {x: -9999, y: -9999};
  let width = 0, height = 0, molecules = [], frame = 0, dark = false;

  const rand = (a, b) => a + Math.random() * (b - a);

  function resize() {
    const dpr = window.devicePixelRatio || 1;
    width = window.innerWidth;
    height = window.innerHeight;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const wanted = Math.max(18, Math.min(55, Math.round(width * height / 30000)));
    while (molecules.length < wanted) molecules.push(spawn());
    molecules.length = wanted;
  }

  function spawn() {
    const depth = rand(0.6, 1.5);
    const bx = rand(-0.25, 0.25), by = rand(-0.2, 0.2);
    return {
      template: TEMPLATES[Math.floor(Math.random() * TEMPLATES.length)],
      colour: COLOURS[Math.floor(Math.random() * COLOURS.length)],
      x: rand(0, width), y: rand(0, height), vx: bx, vy: by, bx, by,
      angle: rand(0, Math.PI * 2), spin: rand(-0.004, 0.004), depth,
    };
  }

  function detectTheme() {
    // Streamlit paints the theme colour on the app container; read it so the molecules
    // stay readable after the user switches between light and dark.
    const colour = getComputedStyle(app).backgroundColor;
    // Cards use it as a mostly opaque base so molecules do not show through their text.
    document.documentElement.style.setProperty("--page-colour", colour);
    const m = colour.match(/\\d+/g);
    dark = m ? (0.299 * m[0] + 0.587 * m[1] + 0.114 * m[2]) < 128 : false;
  }

  function step(m) {
    const dx = m.x - mouse.x, dy = m.y - mouse.y, dist = Math.hypot(dx, dy);
    if (dist < PUSH_RADIUS && dist > 0.01) {
      const push = Math.pow(1 - dist / PUSH_RADIUS, 2) * 1.1;
      m.vx += (dx / dist) * push;
      m.vy += (dy / dist) * push;
      m.spin += (dx > 0 ? 1 : -1) * push * 0.002;
    }
    m.vx += (m.bx - m.vx) * 0.02;
    m.vy += (m.by - m.vy) * 0.02;
    m.spin *= 0.995;
    m.x += m.vx;
    m.y += m.vy;
    m.angle += m.spin;
    const pad = 60;
    if (m.x < -pad) m.x = width + pad; else if (m.x > width + pad) m.x = -pad;
    if (m.y < -pad) m.y = height + pad; else if (m.y > height + pad) m.y = -pad;
  }

  function draw(m) {
    const scale = m.depth, cos = Math.cos(m.angle), sin = Math.sin(m.angle);
    const bondAlpha = dark ? 0.42 : 0.3, atomAlpha = dark ? 0.6 : 0.42;
    const points = m.template.atoms.map(([ax, ay, r]) => [
      m.x + (ax * cos - ay * sin) * scale, m.y + (ax * sin + ay * cos) * scale, r * scale]);
    ctx.lineWidth = 1.8 * scale;
    ctx.strokeStyle = `rgba(${m.colour},${bondAlpha})`;
    for (const [a, b] of m.template.bonds) {
      ctx.beginPath();
      ctx.moveTo(points[a][0], points[a][1]);
      ctx.lineTo(points[b][0], points[b][1]);
      ctx.stroke();
    }
    ctx.fillStyle = `rgba(${m.colour},${atomAlpha})`;
    for (const [px, py, r] of points) {
      ctx.beginPath();
      ctx.arc(px, py, r, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  function drawPointer() {
    if (mouse.x < 0) return;
    const glow = ctx.createRadialGradient(mouse.x, mouse.y, 0, mouse.x, mouse.y, 220);
    glow.addColorStop(0, `rgba(166,108,255,${dark ? 0.16 : 0.12})`);
    glow.addColorStop(1, "rgba(166,108,255,0)");
    ctx.fillStyle = glow;
    ctx.fillRect(mouse.x - 220, mouse.y - 220, 440, 440);
    ctx.lineWidth = 1.2;
    for (const m of molecules) {
      const d = Math.hypot(m.x - mouse.x, m.y - mouse.y);
      if (d >= LINK_RADIUS) continue;
      ctx.strokeStyle = `rgba(${m.colour},${(1 - d / LINK_RADIUS) * 0.5})`;
      ctx.beginPath();
      ctx.moveTo(mouse.x, mouse.y);
      ctx.lineTo(m.x, m.y);
      ctx.stroke();
    }
  }

  function render() {
    ctx.clearRect(0, 0, width, height);
    for (const m of molecules) { if (!reduceMotion) step(m); draw(m); }
    drawPointer();
  }

  function loop() {
    if (!document.hidden) {
      if (frame++ % 30 === 0) detectTheme();
      render();
    }
    requestAnimationFrame(loop);
  }

  document.addEventListener("mousemove", (e) => { mouse.x = e.clientX; mouse.y = e.clientY; });
  document.documentElement.addEventListener("mouseleave", () => { mouse.x = mouse.y = -9999; });
  window.addEventListener("resize", resize);
  resize();
  detectTheme();
  if (reduceMotion) render(); else loop();
})();
</script>
"""


def inject() -> None:
    st.html(SCRIPT, unsafe_allow_javascript=True)
