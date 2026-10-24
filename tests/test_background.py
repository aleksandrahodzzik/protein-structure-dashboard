"""The decorative layer: guard against duplicate loops and against any network access."""

import pytest
import streamlit as st

from protein_dashboard import background, style


def test_script_is_idempotent_and_wrapped_in_one_script_tag() -> None:
    assert background.SCRIPT.count("<script>") == background.SCRIPT.count("</script>") == 1
    assert "if (window.__moleculeBackground) return;" in background.SCRIPT


@pytest.mark.parametrize("api", ["fetch(", "XMLHttpRequest", "WebSocket", "localStorage", "cookie"])
def test_script_touches_no_network_or_storage(api: str) -> None:
    assert api not in background.SCRIPT


def test_script_honours_reduced_motion() -> None:
    assert "prefers-reduced-motion" in background.SCRIPT


def test_inject_runs_the_script_in_the_page(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str, bool]] = []
    monkeypatch.setattr(
        st,
        "html",
        lambda body, **kw: calls.append((body, kw.get("unsafe_allow_javascript", False))),
    )

    background.inject()

    assert calls == [(background.SCRIPT, True)]


def test_section_title_escapes_html(monkeypatch: pytest.MonkeyPatch) -> None:
    rendered: list[str] = []
    monkeypatch.setattr(st, "html", lambda body, **kw: rendered.append(body))

    style.section("<b>Chain</b> & co")

    assert rendered == ['<div class="section-pill">&lt;b&gt;Chain&lt;/b&gt; &amp; co</div>']
