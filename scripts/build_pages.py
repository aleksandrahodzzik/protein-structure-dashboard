"""Build the static GitHub Pages version of the app, which runs in the browser with stlite.

uv run python scripts/build_pages.py            # writes _site/index.html
python -m http.server -d _site 8000             # then open http://localhost:8000

stlite (https://github.com/whitphx/stlite) runs Streamlit on Pyodide inside a web worker,
so the page needs no server: RCSB's APIs answer the browser directly (they send
Access-Control-Allow-Origin: *). The app's sources, its Streamlit config and its
requirements are embedded in one index.html, so the page cannot drift from the code it
was built from. Streamlit itself comes with stlite, at the version that stlite bundles.
"""

from __future__ import annotations

import argparse
import json
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "src" / "protein_dashboard"
OUT_DIR = ROOT / "_site"
STLITE_VERSION = "1.12.0"  # bundles Streamlit 1.65
STLITE_CDN = f"https://cdn.jsdelivr.net/npm/@stlite/browser@{STLITE_VERSION}/build"
TITLE = "Protein Structure Dashboard"
DESCRIPTION = (
    "Explore any structure from the RCSB Protein Data Bank in 3D and analyse its chains "
    "(composition, pI, GRAVY, Kyte-Doolittle profile). Runs entirely in your browser."
)
# Built for Pyodide (biopython has C extensions, requests is patched by stlite), so these
# come from Pyodide's own distribution at whatever version it ships. The pure-Python rest
# is pinned to uv.lock, so the page runs the same plotly and py3Dmol as the server does.
FROM_PYODIDE = frozenset({"biopython", "requests"})
LIVE_DEMO = "https://protein-structure-dashboard.streamlit.app"
REPO = "https://github.com/aleksandrahodzzik/protein-structure-dashboard"

PAGE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <meta name="description" content="__DESCRIPTION__">
  <link rel="stylesheet" href="__CDN__/stlite.css">
  <style>
    html, body { margin: 0; height: 100%; background: #fff; }
    @media (prefers-color-scheme: dark) { html, body { background: #0e1117; } }
    #boot {
      position: fixed; inset: 0; display: grid; place-content: center; gap: 0.75rem;
      font: 16px/1.5 system-ui, sans-serif; color: #31333f; text-align: center;
      padding: 0 16px; pointer-events: none;
    }
    @media (prefers-color-scheme: dark) { #boot { color: #fafafa; } }
    #boot small { opacity: 0.7; }
    #boot a { pointer-events: auto; color: inherit; }
  </style>
</head>
<body>
  <div id="boot" aria-live="polite">
    <strong>Loading __TITLE__&hellip;</strong>
    <small>Python and the app's packages (about 40&nbsp;MB on the first visit) run
      inside your browser, so the app takes 30&ndash;45&nbsp;s to start.</small>
    <small><a href="__REPO__">Source on GitHub</a> &middot;
      <a href="__LIVE_DEMO__">Server-hosted version</a></small>
  </div>
  <div id="root"></div>
  <noscript>This app runs Python in your browser and needs JavaScript.</noscript>
  <script type="application/json" id="app-config">__CONFIG__</script>
  <script type="module">
    import { mount } from "__CDN__/stlite.js";
    const config = JSON.parse(document.getElementById("app-config").textContent);
    mount(config, document.getElementById("root"));
    // stlite draws its own progress toasts and the app over #root; the boot text only
    // covers the blank moment before that, so it goes once the app has rendered.
    const boot = document.getElementById("boot");
    const gone = new MutationObserver(() => {
      if (document.querySelector('[data-testid="stApp"], .stApp')) {
        boot.remove();
        gone.disconnect();
      }
    });
    gone.observe(document.body, { childList: true, subtree: true });
  </script>
</body>
</html>
"""


def canonical(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def requirements(root: Path = ROOT) -> list[str]:
    """The app's runtime dependencies, minus Streamlit, which stlite provides."""
    pyproject = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    lock = tomllib.loads((root / "uv.lock").read_text(encoding="utf-8"))
    locked = {canonical(p["name"]): p["version"] for p in lock["package"]}
    result = []
    for dep in pyproject["project"]["dependencies"]:
        match = re.match(r"[A-Za-z0-9._-]+", dep)
        assert match, dep
        name, key = match.group(0), canonical(match.group(0))
        if key == "streamlit":
            continue
        result.append(dep if key in FROM_PYODIDE else f"{name}=={locked[key]}")
    return result


def app_files() -> dict[str, str]:
    """Virtual file system of the app: the entry point, its config and the package."""
    files = {
        "streamlit_app.py": (ROOT / "streamlit_app.py").read_text(encoding="utf-8"),
        ".streamlit/config.toml": (ROOT / ".streamlit" / "config.toml").read_text(encoding="utf-8"),
    }
    for source in sorted(PACKAGE.glob("*.py")):
        files[f"protein_dashboard/{source.name}"] = source.read_text(encoding="utf-8")
    return files


def mount_config() -> dict[str, object]:
    return {
        "entrypoint": "streamlit_app.py",
        "files": app_files(),
        "requirements": requirements(),
        # The same menu as the deployed app's viewers get, wherever the page is served.
        "streamlitConfig": {"client.toolbarMode": "viewer"},
    }


def render(config: dict[str, object]) -> str:
    # "</" is escaped so no embedded source (some hold <script> strings) can end the
    # JSON <script> block early; "<\\/" is still valid JSON for "</".
    payload = json.dumps(config, ensure_ascii=False).replace("</", "<\\/")
    page = PAGE
    for key, value in {
        "__TITLE__": TITLE,
        "__DESCRIPTION__": DESCRIPTION,
        "__CDN__": STLITE_CDN,
        "__REPO__": REPO,
        "__LIVE_DEMO__": LIVE_DEMO,
    }.items():
        page = page.replace(key, value)
    return page.replace("__CONFIG__", payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=OUT_DIR, help="output directory")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    index = args.out / "index.html"
    index.write_text(render(mount_config()), encoding="utf-8")
    # Serve files as they are: Pages would otherwise run Jekyll over the site.
    (args.out / ".nojekyll").write_text("", encoding="utf-8")
    print(f"Wrote {index} ({index.stat().st_size / 1024:.0f} KiB)")


if __name__ == "__main__":
    main()
