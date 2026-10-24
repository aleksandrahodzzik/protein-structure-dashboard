"""Capture the README screenshots from a running app (live RCSB data, headless Chromium).

uv sync --group screenshots && uv run playwright install chromium
uv run streamlit run streamlit_app.py &
uv run python scripts/screenshots.py --url http://localhost:8501
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

OUT_DIR = Path(__file__).resolve().parents[1] / "docs" / "img"
# Software GL so the 3D view renders on machines (and CI runners) without a GPU.
CHROMIUM_ARGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]
HYDROPHOBICITY = "Hydrophobicity (Kyte-Doolittle)"


def open_entry(page: Page, base_url: str, pdb_id: str) -> None:
    page.goto(f"{base_url}/?pdb={pdb_id}", wait_until="networkidle")
    page.get_by_text("analytics").first.wait_for(timeout=120_000)


def settle_3d(page: Page) -> None:
    """Headless Chromium only produces WebGL frames once a screenshot asks for one,
    so take a throwaway shot and give 3Dmol time to draw before the real capture."""
    time.sleep(5)
    page.locator("iframe").first.screenshot()
    time.sleep(3)


def shoot_3d_view(page: Page, base_url: str) -> None:
    open_entry(page, base_url, "4HHB")
    page.locator("[data-testid=stSidebar] [data-testid=stSelectbox]").nth(1).click()
    page.get_by_role("option", name=HYDROPHOBICITY).click()
    settle_3d(page)
    page.screenshot(
        path=OUT_DIR / "3d-view.png", clip={"x": 0, "y": 0, "width": 1440, "height": 960}
    )


def shoot_charts(page: Page, base_url: str) -> None:
    open_entry(page, base_url, "1BL8")
    page.get_by_role("button", name="Transmembrane preset").click()
    page.get_by_text("Kyte-Doolittle hydrophobicity (window 19)").wait_for(timeout=30_000)
    time.sleep(3)
    top = page.get_by_text("Chain A analytics").bounding_box()
    bottom = page.get_by_role("button", name="Download metrics (CSV)").bounding_box()
    assert top and bottom
    clip = {
        "x": 180,
        "y": top["y"] - 16,
        "width": 1260,
        "height": bottom["y"] + bottom["height"] - top["y"] + 32,
    }
    page.screenshot(path=OUT_DIR / "charts.png", clip=clip)


def shoot_search(page: Page, base_url: str) -> None:
    page.goto(base_url, wait_until="networkidle")
    page.get_by_label("PDB ID or keyword").fill("hemoglobin")
    page.get_by_role("button", name="Search").click()
    page.get_by_text("Top 10 results").wait_for(timeout=60_000)
    # Titles load one request at a time; wait for the last result button.
    page.locator("[data-testid=stSidebar] [data-testid=stButton]").nth(9).wait_for(timeout=60_000)
    time.sleep(2)
    page.screenshot(
        path=OUT_DIR / "search.png", clip={"x": 0, "y": 0, "width": 1440, "height": 900}
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://localhost:8501")
    base_url = parser.parse_args().url.rstrip("/")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(args=CHROMIUM_ARGS)
        for shoot, height in ((shoot_3d_view, 1000), (shoot_charts, 2400), (shoot_search, 1000)):
            page = browser.new_page(viewport={"width": 1440, "height": height})
            shoot(page, base_url)
            page.close()
        browser.close()
    print(f"Wrote screenshots to {OUT_DIR}")


if __name__ == "__main__":
    main()
