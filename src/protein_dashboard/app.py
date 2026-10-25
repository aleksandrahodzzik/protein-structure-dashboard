"""Streamlit UI. Deliberately thin: widgets and layout live here, logic in sibling modules."""

from __future__ import annotations

import textwrap
from datetime import timedelta

import streamlit as st
import streamlit.components.v1 as components

from . import analysis, background, charts, examples, rcsb_client, style, viewer
from .rcsb_client import EntryMetadata, RcsbClient, RcsbError
from .structure import UNKNOWN_RESIDUE, ParsedStructure, StructureParseError, parse_structure

CACHE_TTL = timedelta(hours=24)
VIEWER_HEIGHT = 600

EXAMPLES_PER_ROW = 4

METRIC_HELP = {
    "length": "Modelled residues of the selected chain that are standard amino acids.",
    "mw": "Average molecular weight of the observed sequence, in daltons.",
    "pi": "Isoelectric point: the pH at which the protein carries no net charge.",
    "gravy": "Grand average of hydropathy (Kyte-Doolittle). Positive is hydrophobic, "
    "negative is hydrophilic.",
}


def get_client() -> RcsbClient:
    """Single construction point so tests can swap in a fake client."""
    return RcsbClient()


# --- cached network access ---------------------------------------------------------
# RCSB data changes rarely, so a day-long TTL keeps repeat views instant and keeps us
# well below any rate limit. Exceptions are never cached by Streamlit.


@st.cache_data(ttl=CACHE_TTL, max_entries=200, show_spinner=False)
def load_entry(pdb_id: str, with_organisms: bool = True) -> EntryMetadata:
    return get_client().fetch_entry(pdb_id, with_organisms=with_organisms)


# Structure text is megabytes for large entries, so the number kept in memory is capped.
@st.cache_data(ttl=CACHE_TTL, max_entries=8, show_spinner="Downloading structure from RCSB...")
def load_structure_text(pdb_id: str) -> str:
    return get_client().fetch_structure(pdb_id)


@st.cache_data(ttl=CACHE_TTL, max_entries=8, show_spinner="Parsing structure...")
def load_parsed(pdb_id: str) -> ParsedStructure:
    return parse_structure(load_structure_text(pdb_id), pdb_id)


@st.cache_data(ttl=CACHE_TTL, max_entries=50, show_spinner="Searching RCSB...")
def search_entries(query: str) -> list[str]:
    return get_client().search(query, limit=10)


# --- navigation: the selected entry lives in the URL (?pdb=1UBQ) so it is shareable ----


def select_entry(pdb_id: str | None) -> None:
    if pdb_id:
        st.query_params["pdb"] = pdb_id
    else:
        st.query_params.clear()


def set_kd_window(window: int) -> None:
    st.session_state["kd_window"] = window


def go_home() -> None:
    st.session_state.pop("search_results", None)
    select_entry(None)


def is_dark_theme() -> bool:
    return st.context.theme.get("type") == "dark"


def request_theme_switch() -> None:
    st.session_state["theme_target"] = "Light" if is_dark_theme() else "Dark"


# Streamlit has no API to change the theme from Python. Its own menu does it natively
# (every widget and chart restyles, nothing reloads), so the button presses that menu item.
# st.context.theme only catches up on the next script run, and the 3D viewer (an iframe)
# is coloured from it, so the script also presses the menu's own "Rerun".
THEME_SWITCH_SCRIPT = """
<script>
const doc = window.parent.document;
const menu = doc.querySelector('[data-testid="stMainMenu"]');
const opener = menu && (menu.querySelector("button") || menu);
const press = (id) => {
  const item = doc.querySelector('[data-testid="' + id + '"]');
  if (item) item.click();
};
if (opener) {
  opener.click();
  setTimeout(() => {
    press("stMainMenuItem-theme-%s");
    setTimeout(() => {
      press("stMainMenuItem-rerun");
      setTimeout(() => doc.body.dispatchEvent(
        new KeyboardEvent("keydown", {key: "Escape", bubbles: true})), 150);
    }, 300);
  }, 150);
}
</script>
"""


def render_menu() -> None:
    """Top of the sidebar: back to the example gallery, and the light/dark switch."""
    with st.sidebar:
        home_column, theme_column = st.columns([3, 1], vertical_alignment="center")
        home_column.button(
            "Home",
            icon=":material/home:",
            key="home",
            on_click=go_home,
            width="stretch",
            type="primary",
        )
        dark = is_dark_theme()
        theme_column.button(
            ":material/light_mode:" if dark else ":material/dark_mode:",  # icon-only label
            key="theme_toggle",
            help="Switch to light mode" if dark else "Switch to dark mode",
            on_click=request_theme_switch,
            width="stretch",
        )
        target = st.session_state.pop("theme_target", None)
        if target:
            components.html(THEME_SWITCH_SCRIPT % target, height=0)


def render_search() -> None:
    with st.sidebar:
        style.section("Find a structure")
        with st.form("search"):
            query = st.text_input("PDB ID or keyword", placeholder="1UBQ, or e.g. hemoglobin")
            submitted = st.form_submit_button("Search", width="stretch")
        if submitted and query.strip():
            handle_query(query.strip())
        results = st.session_state.get("search_results")
        if results is not None:
            render_search_results(*results)


def handle_query(query: str) -> None:
    if rcsb_client.is_valid_pdb_id(query):
        st.session_state.pop("search_results", None)
        select_entry(query.upper())
        return
    try:
        ids = search_entries(query)
    except RcsbError as exc:
        st.sidebar.error(str(exc))
        return
    st.session_state["search_results"] = (query, ids)


def render_search_results(query: str, ids: list[str]) -> None:
    if not ids:
        st.warning(f"No structures found for '{query}'.")
        return
    st.caption(f"Top {len(ids)} results for '{query}'")
    for pdb_id in ids:
        try:
            title = load_entry(pdb_id, with_organisms=False).title
        except RcsbError:
            title = ""
        label = (
            f"{pdb_id} · {textwrap.shorten(title.title(), 60, placeholder='…')}"
            if title
            else pdb_id
        )
        st.button(label, key=f"hit_{pdb_id}", on_click=select_entry, args=(pdb_id,))


# --- pages ----------------------------------------------------------------------------


def render_landing() -> None:
    st.markdown(
        "Fetch any entry from the [RCSB Protein Data Bank](https://www.rcsb.org), explore it in "
        "3D, and inspect its sequence composition and hydrophobicity. Search in the sidebar or "
        f"start with one of {len(examples.ALL)} examples:"
    )
    for tab, category in zip(
        st.tabs([c.name for c in examples.CATEGORIES]), examples.CATEGORIES, strict=True
    ):
        with tab:
            render_example_grid(category.items)


def render_example_grid(items: tuple[examples.Example, ...]) -> None:
    for start in range(0, len(items), EXAMPLES_PER_ROW):
        row = items[start : start + EXAMPLES_PER_ROW]
        for column, item in zip(st.columns(EXAMPLES_PER_ROW), row, strict=False):
            with column, st.container(border=True, height=220, key=f"card_{item.pdb_id}"):
                st.markdown(f"**{item.name}**")
                st.caption(item.blurb)
                st.button(
                    item.pdb_id,
                    key=f"example_{item.pdb_id}",
                    on_click=select_entry,
                    args=(item.pdb_id,),
                )


def render_header(entry: EntryMetadata) -> None:
    st.subheader(entry.title or entry.pdb_id)
    resolution = f"{entry.resolution:.2f} Å" if entry.resolution else "resolution n/a"
    organism = ", ".join(entry.organisms) or "organism n/a"
    st.caption(
        f"**{entry.pdb_id}** · {entry.method.title()} · {resolution} · *{organism}* · "
        f"deposited {entry.deposition_date or 'n/a'} · [View on RCSB]({entry.url})"
    )


def render_controls(
    pdb_id: str, parsed: ParsedStructure
) -> tuple[str, int | None, viewer.ViewerOptions]:
    """Sidebar controls that depend on the loaded structure."""
    with st.sidebar:
        st.divider()
        style.section("Analysis")
        chain_id = st.selectbox(
            "Chain",
            [c.chain_id for c in parsed.chains],
            format_func=lambda cid: f"{cid} ({parsed.chain(cid).residue_count} residues)",
            key=f"chain_{pdb_id}",
        )
        usable = parsed.chain(chain_id).residue_count - parsed.chain(chain_id).unknown_count
        max_window = analysis.max_window(usable)
        window = render_window_controls(max_window)

        style.section("3D view")
        mode = st.selectbox("Colour by", list(viewer.ColorMode))
        surface = st.checkbox("Show surface", disabled=viewer.is_heavy(parsed))
        ligands = st.checkbox(
            "Show ligands as sticks", value=True, disabled=viewer.is_heavy(parsed)
        )
    options = viewer.ViewerOptions(
        color_mode=mode or viewer.ColorMode.CHAIN,  # selectbox is typed Optional
        show_surface=surface,
        show_ligands=ligands,
        focus_chain=chain_id if len(parsed.chains) > 1 else None,
    )
    return chain_id, window, options


def render_window_controls(max_window: int) -> int | None:
    if max_window == 0:
        st.caption("Chain too short for a hydrophobicity profile.")
        return None
    if max_window == analysis.MIN_WINDOW:
        return max_window
    st.session_state.setdefault("kd_window", analysis.DEFAULT_WINDOW)
    st.session_state["kd_window"] = min(st.session_state["kd_window"], max_window)
    window = st.slider(
        "Hydrophobicity window (odd)",
        analysis.MIN_WINDOW,
        max_window,
        step=2,
        key="kd_window",
    )
    preset = min(analysis.TM_WINDOW, max_window)
    st.button(
        f"Transmembrane preset (window {preset})",
        on_click=set_kd_window,
        args=(preset,),
        help="Window 19 with a 1.6 cut-off is the classic Kyte-Doolittle test for "
        "transmembrane helices.",
    )
    return int(window)


def render_metrics(metrics: analysis.ProteinMetrics) -> None:
    cols = st.columns(4)
    cols[0].metric("Length", f"{metrics.length} aa", help=METRIC_HELP["length"])
    cols[1].metric(
        "Molecular weight", f"{metrics.molecular_weight:,.1f} Da", help=METRIC_HELP["mw"]
    )
    cols[2].metric(
        "Isoelectric point (pI)", f"{metrics.isoelectric_point:.2f}", help=METRIC_HELP["pi"]
    )
    cols[3].metric("GRAVY", f"{metrics.gravy:+.3f}", help=METRIC_HELP["gravy"])
    st.caption(
        f"Instability index {metrics.instability_index:.1f} "
        f"({'unstable' if metrics.instability_index > 40 else 'stable'} in vitro, cut-off 40) · "
        f"aromaticity {metrics.aromaticity:.3f}"
    )


def render_entry(pdb_id: str, entry: EntryMetadata, parsed: ParsedStructure) -> None:
    render_header(entry)
    chain_id, window, options = render_controls(pdb_id, parsed)
    chain = parsed.chain(chain_id)

    if viewer.is_heavy(parsed):
        st.warning(
            f"Large structure ({parsed.atom_count:,} atoms). To keep the browser responsive "
            "only the backbone trace is drawn; surface and ligands are disabled."
        )
    style.section("3D structure")
    view = viewer.build_view(load_structure_text(pdb_id), parsed, options, VIEWER_HEIGHT)
    viewer.show_view(view, VIEWER_HEIGHT)
    if options.color_mode is viewer.ColorMode.HYDROPHOBICITY:
        st.caption("Per-residue Kyte-Doolittle value: blue is hydrophilic, orange is hydrophobic.")
    if chain.ligand_names:
        st.caption(f"Ligands on chain {chain_id}: {', '.join(chain.ligand_names)}")

    style.section(f"Chain {chain_id} analytics")
    st.caption(
        f"Observed sequence: {chain.residue_count} modelled residues in {chain.segment_count} "
        "contiguous segment(s). Residues missing from the density are not included, so this "
        "can be shorter than the full (SEQRES) sequence."
    )
    if chain.unknown_count:
        st.info(f"{chain.unknown_count} non-standard residue(s) (X) are excluded from the metrics.")
    try:
        metrics = analysis.compute_metrics(chain.sequence)
    except ValueError as exc:
        st.warning(str(exc))
        return

    render_metrics(metrics)
    st.plotly_chart(charts.composition_chart(metrics.composition), key="composition")
    if window is not None:
        usable = [r for r in chain.residues if r.letter != UNKNOWN_RESIDUE]
        sequence = "".join(r.letter for r in usable)
        profile = analysis.kd_profile(sequence, window)
        st.plotly_chart(
            charts.kd_profile_chart(profile, [r.number for r in usable], sequence),
            key="kd_profile",
        )
    st.download_button(
        "Download metrics (CSV)",
        data=analysis.metrics_to_csv(f"{pdb_id}_{chain_id}", metrics),
        file_name=f"{pdb_id}_{chain_id}_metrics.csv",
        mime="text/csv",
    )


def main() -> None:
    st.set_page_config(page_title="Protein Structure Dashboard", layout="wide")
    style.inject()
    background.inject()
    st.title("Protein Structure Dashboard")
    render_menu()
    render_search()

    requested = st.query_params.get("pdb")
    if not requested:
        render_landing()
        return
    try:
        pdb_id = rcsb_client.normalise_pdb_id(requested)
        entry = load_entry(pdb_id)
        parsed = load_parsed(pdb_id)
    except (RcsbError, StructureParseError) as exc:
        st.error(str(exc))
        st.button("Back to examples", on_click=select_entry, args=(None,))
        return
    render_entry(pdb_id, entry, parsed)
