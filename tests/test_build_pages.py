"""The GitHub Pages build: one index.html that mounts the app in stlite."""

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
_spec = importlib.util.spec_from_file_location("build_pages", ROOT / "scripts" / "build_pages.py")
assert _spec and _spec.loader
build_pages = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build_pages)

CONFIG_BLOCK = re.compile(
    r'<script type="application/json" id="app-config">(.*?)</script>', re.DOTALL
)


def embedded_config(page: str) -> dict:
    match = CONFIG_BLOCK.search(page)
    assert match
    return json.loads(match.group(1))


def test_every_module_of_the_package_is_mounted() -> None:
    config = build_pages.mount_config()
    modules = {p.name for p in (ROOT / "src" / "protein_dashboard").glob("*.py")}
    mounted = {
        path.removeprefix("protein_dashboard/")
        for path in config["files"]
        if path.startswith("protein_dashboard/")
    }
    assert mounted == modules
    assert config["entrypoint"] in config["files"]
    assert ".streamlit/config.toml" in config["files"]


def test_requirements_leave_streamlit_to_stlite_and_pin_pure_python() -> None:
    reqs = build_pages.requirements()
    names = [re.match(r"[A-Za-z0-9._-]+", r).group(0).lower() for r in reqs]  # type: ignore[union-attr]
    assert "streamlit" not in names
    assert set(names) == {"biopython", "plotly", "py3dmol", "requests"}
    for req in reqs:
        name = req.split("=")[0].split(">")[0].lower()
        pinned = "==" in req
        assert pinned == (name not in build_pages.FROM_PYODIDE), req


def test_embedded_sources_cannot_end_the_config_block_early() -> None:
    config = {"files": {"x.py": "html = '</script><script>alert(1)</script>'"}}
    assert embedded_config(build_pages.render(config)) == config


def test_page_round_trips_the_real_app() -> None:
    config = build_pages.mount_config()
    assert embedded_config(build_pages.render(config)) == config


def test_main_writes_the_site(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["build_pages.py", "--out", str(tmp_path)])
    build_pages.main()
    assert (tmp_path / ".nojekyll").exists()
    page = (tmp_path / "index.html").read_text(encoding="utf-8")
    assert f"@stlite/browser@{build_pages.STLITE_VERSION}/build/stlite.js" in page
