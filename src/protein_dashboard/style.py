"""Visual polish kept apart from app logic: one stylesheet and one heading helper.

Colours that must work in both light and dark mode are either the accent (readable on
either background) or translucent greys, so no per-theme branching is needed here.
"""

from __future__ import annotations

from html import escape

import streamlit as st

CSS = """
<style>
:root {
  --accent: #6c5ce7;
  --accent-2: #a66cff;
  --card-border: rgba(128, 128, 128, 0.28);
  --card-fill: rgba(128, 128, 128, 0.07);
  --card-base: color-mix(in srgb, var(--page-colour, transparent) 88%, transparent);
}

/* Animated molecule canvas (see background.py): sits above the page colour, below all content */
[data-testid="stApp"] { isolation: isolate; }
#molecule-background {
  position: fixed;
  inset: 0;
  width: 100vw;
  height: 100vh;
  z-index: -1;
  pointer-events: none;
}

/* Page title: gradient text */
[data-testid="stMain"] h1 {
  background: linear-gradient(90deg, var(--accent), var(--accent-2) 70%, #ff7eb6);
  -webkit-background-clip: text;
  background-clip: text;
  -webkit-text-fill-color: transparent;
  font-weight: 800;
  letter-spacing: -0.5px;
}

/* Section headings as accent pills */
.section-pill {
  display: inline-block;
  margin: 1.1rem 0 0.7rem;
  padding: 0.45rem 1.15rem;
  border-radius: 999px;
  background: linear-gradient(135deg, var(--accent), var(--accent-2));
  box-shadow: 0 4px 14px rgba(108, 92, 231, 0.35);
  color: #fff;
  font-size: 1.2rem;
  font-weight: 700;
  line-height: 1.3;
  letter-spacing: 0.2px;
}
[data-testid="stSidebar"] .section-pill {
  font-size: 1.05rem;
  padding: 0.4rem 1rem;
}

/* Buttons */
.stButton button, .stDownloadButton button, [data-testid="stFormSubmitButton"] button {
  border-radius: 12px;
  font-weight: 600;
  transition: transform 0.12s ease, box-shadow 0.12s ease, border-color 0.12s ease;
}
.stButton button:hover, .stDownloadButton button:hover,
[data-testid="stFormSubmitButton"] button:hover {
  transform: translateY(-1px);
  box-shadow: 0 4px 12px rgba(108, 92, 231, 0.25);
  border-color: var(--accent);
}

/* Metric cards */
[data-testid="stMetric"] {
  background: linear-gradient(var(--card-fill), var(--card-fill)), var(--card-base);
  border: 1px solid var(--card-border);
  border-radius: 16px;
  padding: 0.9rem 1.15rem;
  backdrop-filter: blur(8px);  /* frosted glass over the molecules */
}
[data-testid="stMetricValue"] { font-weight: 700; }

/* Example cards */
[class*="st-key-card_"] {
  border-radius: 18px;
  background: linear-gradient(var(--card-fill), var(--card-fill)), var(--card-base);
  backdrop-filter: blur(8px);
  transition: transform 0.15s ease, box-shadow 0.15s ease, border-color 0.15s ease;
}
[class*="st-key-card_"]:hover {
  transform: translateY(-3px);
  border-color: var(--accent);
  box-shadow: 0 10px 26px rgba(108, 92, 231, 0.22);
}

/* 3D viewer: always a white canvas, framed like a card */
[data-testid="stMain"] iframe {
  background: #fff;
  border: 1px solid var(--card-border);
  border-radius: 18px;
  box-shadow: 0 8px 28px rgba(0, 0, 0, 0.14);
}

/* Sidebar search form */
[data-testid="stSidebar"] [data-testid="stForm"] { border-radius: 16px; }
</style>
"""


def inject() -> None:
    # st.html adds no layout gap, unlike an st.markdown element holding only a <style>.
    st.html(CSS)


def section(title: str) -> None:
    """A large accent pill used as a section heading (in the main area or the sidebar)."""
    st.html(f'<div class="section-pill">{escape(title)}</div>')
