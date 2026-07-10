#!/usr/bin/env python
"""Capture screenshots of the QuickLook web GUI for the docs.

Internal dev tool -- NOT part of the installed package (no console entry
point in pyproject.toml). It boots the Flask app headless, drives a
headless Chromium via Playwright to screenshot each route, then downscales
and palette-quantizes the PNGs with Pillow so they are small enough to
commit under docs/img/.

Requirements (install into your dev venv, not the package):

    pip install playwright pillow
    python -m playwright install chromium

Usage:

    python scripts/screenshot_gui.py                  # reuse or boot a server
    python scripts/screenshot_gui.py --url http://127.0.0.1:5000  # your ql-gui
    python scripts/screenshot_gui.py --port 5070      # boot on a custom port
    python scripts/screenshot_gui.py --no-optimize    # keep full-res PNGs

If a GUI is already serving on the target port (or --url is given), the
script reuses it -- instant, no ~40s import. Otherwise it boots a throwaway
server and tears it down afterwards, so the script stays reproducible in CI.

The default routes overwrite docs/img/ql-gui.png (landing form) and
docs/img/ql-gui-gallery.png (gallery), which the README and docs/gui.md
already reference.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
IMG_DIR = REPO_ROOT / "docs" / "img"

# Importing quicklook.app.app pulls in heavy astro stacks (~40s) and may make
# network calls, so give the server plenty of head start before giving up.
SERVER_BOOT_TIMEOUT = 180

# route -> output filename under docs/img/
DEFAULT_SHOTS = {
    "/": "ql-gui.png",
    "/gallery": "ql-gui-gallery.png",
}


def start_server(port: int, log_path: Path) -> subprocess.Popen:
    """Launch the Flask app in a subprocess, tee-ing its output to log_path."""
    env = dict(os.environ)
    # numexpr logs a non-fatal error if its cap is below the detected core
    # count; set it high so the message is suppressed on most machines.
    env.setdefault("NUMEXPR_MAX_THREADS", "64")
    code = (
        "from quicklook.app.app import app; "
        f"app.run(host='127.0.0.1', port={port}, debug=False, threaded=True)"
    )
    log = open(log_path, "w")
    return subprocess.Popen(
        [sys.executable, "-c", code],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=log,
        stderr=subprocess.STDOUT,
    )


def is_up(base: str) -> bool:
    """True if a GUI already answers HTTP 200 on the landing route."""
    try:
        with urllib.request.urlopen(base + "/", timeout=2) as resp:
            return resp.status == 200
    except Exception:
        return False


def wait_for_server(server: subprocess.Popen, port: int, log_path: Path) -> None:
    """Poll the landing route until HTTP 200, failing fast if the server dies."""
    url = f"http://127.0.0.1:{port}/"
    deadline = time.time() + SERVER_BOOT_TIMEOUT
    while time.time() < deadline:
        if server.poll() is not None:
            tail = log_path.read_text(errors="replace").strip().splitlines()[-15:]
            raise RuntimeError("GUI server process exited early:\n" + "\n".join(tail))
        try:
            with urllib.request.urlopen(url, timeout=2) as resp:
                if resp.status == 200:
                    return
        except Exception:
            time.sleep(1)
    raise RuntimeError(
        f"GUI server did not come up on port {port} within {SERVER_BOOT_TIMEOUT}s (see {log_path})"
    )


def capture(base: str, shots: dict[str, str]) -> list[Path]:
    """Screenshot each route full-page at 2x and return the saved paths."""
    from playwright.sync_api import sync_playwright

    IMG_DIR.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []
    base = base.rstrip("/")
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 900}, device_scale_factor=2)
        for route, name in shots.items():
            dest = IMG_DIR / name
            page.goto(base + route, wait_until="networkidle", timeout=60000)
            page.wait_for_timeout(1200)  # let fonts/layout settle
            page.screenshot(path=str(dest), full_page=True)
            print(f"captured {route} -> {dest.relative_to(REPO_ROOT)}")
            saved.append(dest)
        browser.close()
    return saved


def optimize(paths: list[Path]) -> None:
    """Halve the 2x capture and palette-quantize so PNGs are commit-sized."""
    from PIL import Image

    for path in paths:
        before = path.stat().st_size
        img = Image.open(path).convert("RGB")
        w, h = img.size
        img = img.resize((w // 2, h // 2), Image.LANCZOS)
        pal = img.quantize(colors=256, method=Image.FASTOCTREE, dither=Image.Dither.NONE)
        pal.save(path, optimize=True)
        after = path.stat().st_size
        print(
            f"optimized {path.name}: {before // 1024} KB -> {after // 1024} KB "
            f"({img.size[0]}x{img.size[1]})"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Avoid Chromium's ERR_UNSAFE_PORT list (e.g. 5060/SIP, 6000/X11); 5050 is fine.
    parser.add_argument("--port", type=int, default=5050, help="port to boot the GUI on")
    parser.add_argument(
        "--url",
        default=None,
        help="screenshot an already-running GUI at this base URL "
        "(e.g. http://127.0.0.1:5000); skips booting a server",
    )
    parser.add_argument(
        "--no-optimize",
        action="store_true",
        help="skip the downscale/quantize step (keeps full-res PNGs)",
    )
    args = parser.parse_args()

    # Decide whether to reuse a running GUI or boot a throwaway one.
    if args.url:
        base = args.url.rstrip("/")
        if not is_up(base):
            print(f"ERROR: no GUI responding at {base}", file=sys.stderr)
            return 1
        print(f"reusing running GUI at {base}")
        server = None
    else:
        base = f"http://127.0.0.1:{args.port}"
        if is_up(base):
            print(f"reusing running GUI at {base}")
            server = None
        else:
            log_path = Path(tempfile.gettempdir()) / f"ql_gui_screenshot_{args.port}.log"
            print(f"booting throwaway GUI on port {args.port} (~40s of imports)...")
            server = start_server(args.port, log_path)
            wait_for_server(server, args.port, log_path)

    try:
        saved = capture(base, DEFAULT_SHOTS)
    finally:
        if server is not None:
            server.terminate()
            try:
                server.wait(timeout=10)
            except subprocess.TimeoutExpired:
                server.kill()

    if not args.no_optimize:
        optimize(saved)
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
