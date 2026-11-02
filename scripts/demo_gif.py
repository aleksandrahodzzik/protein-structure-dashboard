"""Record the README demo GIF from a running app: search -> 3D view -> hydrophobicity plot.

uv sync --group screenshots && uv run playwright install chromium
uv run python scripts/demo_gif.py --url https://protein-structure-dashboard.streamlit.app/~/+/

On Streamlit Community Cloud the app itself is served under /~/+/ (the bare URL wraps it
in an iframe). A first, unrecorded pass wakes the app and fills its caches, so the
recording shows the app's speed rather than RCSB's. Playwright records a video, and the
ffmpeg bundled with imageio-ffmpeg turns it into a palette-optimised GIF.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
import time
from pathlib import Path

import imageio_ffmpeg
from playwright.sync_api import Page, sync_playwright

OUT = Path(__file__).resolve().parents[1] / "docs" / "img" / "demo.gif"
# Software GL so the 3D view renders on machines (and CI runners) without a GPU.
CHROMIUM_ARGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"]
VIEWPORT = {"width": 1280, "height": 720}
QUERY = "bacteriorhodopsin"
PDB_ID = "1KG8"  # a hit for the query, with a disordered loop that breaks the profile
HYDROPHOBICITY = "Hydrophobicity (Kyte-Doolittle)"
# Viewport positions of points along the profile's SVG line(s).
PROFILE_POINTS_JS = """el => {
  const points = [];
  for (const path of el.querySelectorAll(".scatterlayer path.js-line")) {
    const total = path.getTotalLength(), m = path.getScreenCTM();
    for (let d = 0; d <= total; d += 6) {
      const p = path.getPointAtLength(d);
      points.push([p.x * m.a + p.y * m.c + m.e, p.x * m.b + p.y * m.d + m.f]);
    }
  }
  return points;
}"""
SECONDS = 14.0  # plus a one-second hold on the last frame: 15 s in all
HOLD_SECONDS = 1.0
GIF_WIDTH = 960
FPS = 12


def wait_idle(page: Page) -> None:
    """Streamlit is done when its 'running' indicator has gone."""
    time.sleep(0.3)
    page.locator("[data-testid=stStatusWidget]").wait_for(state="hidden", timeout=180_000)


def open_app(page: Page, url: str) -> None:
    page.goto(url, wait_until="domcontentloaded")
    wake = page.get_by_role("button", name="Yes, get this app back up!")
    title = page.get_by_text("Protein Structure Dashboard").first
    title.or_(wake).wait_for(timeout=120_000)
    if wake.is_visible():
        wake.click()
    title.wait_for(timeout=300_000)
    page.get_by_role("button", name="1UBQ").wait_for(timeout=120_000)  # the gallery is drawn
    wait_idle(page)


def search_and_open(page: Page, *, typing_delay: float) -> None:
    page.get_by_label("PDB ID or keyword").click()
    page.get_by_label("PDB ID or keyword").type(QUERY, delay=typing_delay)
    page.get_by_role("button", name="Search").click()
    page.get_by_text("Top 10 results").wait_for(timeout=120_000)
    wait_idle(page)
    time.sleep(0.6)
    page.get_by_role("button", name=PDB_ID).click()
    page.get_by_text("analytics").first.wait_for(timeout=180_000)
    wait_idle(page)


def colour_by_hydrophobicity(page: Page) -> None:
    page.locator("[data-testid=stSidebar] [data-testid=stSelectbox]").nth(1).click()
    page.get_by_role("option", name=HYDROPHOBICITY).click()
    wait_idle(page)


def rotate_molecule(page: Page) -> None:
    box = page.locator("[data-testid=stMain] iframe").first.bounding_box()
    assert box
    x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
    page.mouse.move(x - 120, y)
    page.mouse.down()
    for step in range(40):
        page.mouse.move(x - 120 + step * 7, y + (step % 20 - 10) * 1.5)
        time.sleep(0.03)
    page.mouse.up()


def show_profile(page: Page) -> None:
    page.get_by_role("button", name="Transmembrane preset").click()
    # Off the button and out of focus, so its help tooltip closes.
    page.mouse.move(VIEWPORT["width"] - 40, 80)
    page.evaluate("document.activeElement && document.activeElement.blur()")
    page.get_by_text("Kyte-Doolittle hydrophobicity (window 19)").wait_for(timeout=60_000)
    wait_idle(page)
    # A smooth scroll so the GIF shows the page moving; a mouse wheel over the 3D view
    # would zoom the molecule instead.
    chart = page.locator(".stPlotlyChart").nth(1)
    chart.evaluate("el => el.scrollIntoView({behavior: 'smooth', block: 'center'})")
    time.sleep(1.5)
    # Follow the drawn line up to its highest peak, so Plotly's hover label shows.
    points = chart.evaluate(PROFILE_POINTS_JS)
    peak = min(range(len(points)), key=lambda i: points[i][1])
    for x, y in points[max(0, peak - 14) : peak + 1]:
        page.mouse.move(x, y)
        time.sleep(0.05)
    time.sleep(1.5)


def flash_start_marker(page: Page) -> None:
    """Black out the page for a moment. ffmpeg finds that flash in the video and cuts just
    after it: wall-clock times drift from the video's own timeline by a second or two."""
    page.evaluate(
        """() => {
          const cover = document.createElement("div");
          cover.id = "gif-start-marker";
          cover.style.cssText = "position:fixed;inset:0;background:#000;z-index:2147483647";
          document.body.append(cover);
        }"""
    )
    time.sleep(1.0)  # at 0.6 s the recording sometimes missed it
    page.evaluate("document.getElementById('gif-start-marker').remove()")


def record(url: str, video_dir: Path) -> Path:
    """Warm the app up, then record the demo; returns the video file."""
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(args=CHROMIUM_ARGS)

        warm = browser.new_page(viewport=VIEWPORT)
        open_app(warm, url)
        search_and_open(warm, typing_delay=0)
        show_profile(warm)
        warm.close()

        context = browser.new_context(
            viewport=VIEWPORT, record_video_dir=str(video_dir), record_video_size=VIEWPORT
        )
        page = context.new_page()
        open_app(page, url)
        flash_start_marker(page)
        time.sleep(0.8)
        search_and_open(page, typing_delay=70)
        colour_by_hydrophobicity(page)
        # Wake the WebGL frames. A whole-page screenshot, because an element screenshot
        # leaks into the video as a few frames of that element alone on grey.
        page.screenshot()
        time.sleep(1.0)
        rotate_molecule(page)
        time.sleep(0.5)
        show_profile(page)
        video = page.video
        context.close()  # the video ends here
        browser.close()
        assert video
        return Path(video.path())


def video_seconds(ffmpeg: str, video: Path) -> float:
    probe = subprocess.run([ffmpeg, "-i", str(video)], capture_output=True, text=True)
    match = re.search(r"Duration: (\d+):(\d+):([\d.]+)", probe.stderr)
    assert match, probe.stderr
    hours, minutes, seconds = match.groups()
    return int(hours) * 3600 + int(minutes) * 60 + float(seconds)


def after_start_marker(ffmpeg: str, video: Path) -> float:
    """Seconds into the video where the black start marker ends."""
    probe = subprocess.run(
        [ffmpeg, "-i", str(video), "-vf", "blackdetect=d=0.2:pix_th=0.1", "-an", "-f", "null", "-"],
        capture_output=True,
        text=True,
    )
    ends = re.findall(r"black_end:([\d.]+)", probe.stderr)
    assert ends, "the black start marker was not found in the video"
    return float(ends[-1]) + 0.1  # nothing in the app itself is ever all black


def to_gif(video: Path, out: Path) -> None:
    """Cut to the demo, fit it into SECONDS, hold the last frame, and encode with one
    shared palette."""
    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    begin, end = after_start_marker(ffmpeg, video), video_seconds(ffmpeg, video)
    speed = max(1.0, (end - begin) / SECONDS)
    filters = (
        f"setpts=PTS/{speed:.3f},fps={FPS},scale={GIF_WIDTH}:-1:flags=lanczos,"
        f"tpad=stop_mode=clone:stop_duration={HOLD_SECONDS},"
        "split[a][b];[a]palettegen=max_colors=128:stats_mode=diff[p];"
        "[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            ffmpeg,
            "-y",
            "-ss",
            f"{begin:.2f}",
            "-to",
            f"{end:.2f}",
            "-i",
            str(video),
            "-vf",
            filters,
            "-loop",
            "0",
            str(out),
        ],
        check=True,
        capture_output=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://localhost:8501")
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        to_gif(record(args.url, Path(tmp)), args.out)
    print(f"Wrote {args.out} ({args.out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
