"""
PDB Pasqua — Easter egg art from Protein Data Bank structures.
"""

from __future__ import annotations

import io
import os
import re
from dataclasses import dataclass

# Writable cache (avoids issues when ~/.matplotlib is not writable)
_mpl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".mplconfig")
os.makedirs(_mpl_dir, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", _mpl_dir)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import requests
import streamlit as st
from matplotlib import colormaps
from matplotlib.collections import LineCollection

RCSB_DOWNLOAD = "https://files.rcsb.org/download/{id}.pdb"

# Pastel egg shells + accent palettes (one set per slot)
EGG_FILLS = ["#fff5f7", "#f2fdf5", "#f7f3ff"]
EGG_EDGES = ["#e8a0b4", "#7ec9a0", "#b8a0e8"]
RIBBON_COLORS = ["#d4849a", "#5cb88a", "#8b7fd4"]
SCATTER_CMAPS = ["spring", "summer", "cool"]


@dataclass
class ParsedPDB:
    pdb_id: str
    ca_coords: np.ndarray  # (N, 3), sequence order


def normalize_pdb_id(raw: str) -> str | None:
    s = raw.strip().upper()
    if not s:
        return None
    if not re.fullmatch(r"[0-9A-Z]{4}", s):
        return None
    return s


def fetch_pdb_text(pdb_id: str) -> str:
    url = RCSB_DOWNLOAD.format(id=pdb_id)
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text


def parse_ca_trace(pdb_text: str, pdb_id: str) -> ParsedPDB:
    """Extract CA atoms in file order (typical N→C per chain)."""
    coords: list[list[float]] = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        coords.append([x, y, z])
    arr = np.array(coords, dtype=np.float64)
    return ParsedPDB(pdb_id=pdb_id, ca_coords=arr)


def pca_project_xy(coords: np.ndarray) -> np.ndarray:
    """First two principal components (N, 2)."""
    if coords.shape[0] < 2:
        return np.zeros((coords.shape[0], 2))
    X = coords - coords.mean(axis=0)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    xy = X @ vt[:2].T
    m = np.abs(xy).max()
    if m < 1e-9:
        return xy
    return xy / m


def egg_curve(n: int = 320) -> tuple[np.ndarray, np.ndarray]:
    t = np.linspace(0, 2 * np.pi, n)
    # Slightly asymmetric egg silhouette
    x = 0.95 * np.cos(t)
    y = 1.18 * np.sin(t) * (1.0 + 0.22 * np.cos(t))
    return x, y


def render_egg_figure(parsed: ParsedPDB, style_index: int) -> io.BytesIO:
    """Render one egg PNG into a BytesIO buffer."""
    ca = parsed.ca_coords
    pdb_id = parsed.pdb_id
    idx = style_index % 3

    fig, ax = plt.subplots(figsize=(4.2, 5.2), dpi=160)
    bg = "#fff8fb" if idx == 0 else "#f8fffa" if idx == 1 else "#faf8ff"
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)

    ex, ey = egg_curve()
    ax.fill(ex, ey, color=EGG_FILLS[idx], edgecolor=EGG_EDGES[idx], linewidth=2.8, zorder=1)

    if ca.shape[0] < 4:
        ax.text(
            0,
            0,
            "Not enough\nCα atoms",
            ha="center",
            va="center",
            fontsize=11,
            color="#a08090",
            zorder=5,
        )
    else:
        xy = pca_project_xy(ca) * 0.78
        n = len(xy)
        # Soft ribbon: line segments colored along sequence
        segs = np.stack([xy[:-1], xy[1:]], axis=1)
        t = np.linspace(0.15, 0.95, len(segs))
        cmap = colormaps[SCATTER_CMAPS[idx]]
        colors = cmap(t)
        lc = LineCollection(segs, colors=colors, linewidths=2.4, alpha=0.85, capstyle="round", zorder=3)
        ax.add_collection(lc)
        # Sparkle scatter on top
        sizes = 14 + 10 * np.sin(np.linspace(0, 3 * np.pi, n)) ** 2
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            c=np.arange(n),
            cmap=SCATTER_CMAPS[idx],
            s=sizes,
            alpha=0.9,
            edgecolors="white",
            linewidths=0.45,
            zorder=4,
        )

    # Tiny decorative dots (speckles)
    rng = np.random.default_rng(abs(hash(pdb_id)) % (2**32))
    for _ in range(28):
        te = rng.uniform(0, 2 * np.pi)
        rx = 0.95 * np.cos(te)
        ry = 1.18 * np.sin(te) * (1.0 + 0.22 * np.cos(te))
        if ry > -0.35 and abs(rx) > 0.12:
            continue
        ax.scatter(
            [rx * rng.uniform(0.82, 0.98)],
            [ry * rng.uniform(0.82, 0.98)],
            s=rng.uniform(4, 14),
            c=[RIBBON_COLORS[idx]],
            alpha=0.35,
            zorder=2,
        )

    ax.text(
        0,
        -1.38,
        pdb_id,
        ha="center",
        va="top",
        fontsize=13,
        fontweight="600",
        color=RIBBON_COLORS[idx],
        fontfamily="sans-serif",
        zorder=6,
    )
    ax.text(
        0,
        1.32,
        "PDB Pasqua",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#c9a8b8",
        fontfamily="sans-serif",
        style="italic",
        zorder=6,
    )

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.55, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout(pad=0.2)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=bg, edgecolor="none", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    buf.seek(0)
    return buf


def main() -> None:
    st.set_page_config(
        page_title="PDB Pasqua",
        page_icon="🥚",
        layout="centered",
    )

    st.markdown(
        """
        <style>
        .main-header { font-size: 2.1rem; font-weight: 700; color: #8b5a7c; margin-bottom: 0.2rem; }
        .subtle { color: #9a8a96; font-size: 1rem; }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.markdown('<p class="main-header">PDB Pasqua</p>', unsafe_allow_html=True)
    st.markdown(
        '<p class="subtle">Turn real protein structures into adorable Easter Egg art!!! Enter a four-character PDB ID (e.g. 1CRN, 2LYZ) to begin exploring!. These structures are derived from the secondary structure of the protein!</p>',
        unsafe_allow_html=True,
    )

    c1, c2, c3 = st.columns(3)
    with c1:
        id1 = st.text_input("Egg 1 — PDB ID", placeholder="1CRN", key="e1")
    with c2:
        id2 = st.text_input("Egg 2 — PDB ID", placeholder="2LYZ", key="e2")
    with c3:
        id3 = st.text_input("Egg 3 — PDB ID", placeholder="1PGB", key="e3")

    go = st.button("Hatch my eggs!", type="primary")

    if not go:
        st.info("Add one to three PDB IDs and click **Hatch my eggs!** to generate PNGs.")
        return

    raw_ids = [id1, id2, id3]
    normalized: list[tuple[str, str | None]] = []
    for r in raw_ids:
        n = normalize_pdb_id(r)
        normalized.append((r.strip(), n))

    errors: list[str] = []
    valid: list[str] = []
    for raw, nid in normalized:
        if not raw:
            continue
        if nid is None:
            errors.append(f"`{raw}` is not a valid PDB ID (use four letters/numbers).")
        elif nid in valid:
            errors.append(f"`{nid}` was entered more than once — skipping duplicate.")
        else:
            valid.append(nid)

    if errors:
        for e in errors:
            st.warning(e)

    if not valid:
        st.error("Please enter at least one valid four-character PDB ID.")
        return

    st.markdown("---")

    for i, pdb_id in enumerate(valid):
        with st.spinner(f"Fetching {pdb_id} from RCSB…"):
            try:
                text = fetch_pdb_text(pdb_id)
                parsed = parse_ca_trace(text, pdb_id)
            except requests.RequestException as e:
                st.error(f"Could not download **{pdb_id}**: {e}")
                continue

        if parsed.ca_coords.shape[0] == 0:
            st.error(f"**{pdb_id}**: no Cα atoms found in the file.")
            continue

        buf = render_egg_figure(parsed, style_index=i)
        data = buf.getvalue()
        st.subheader(f"Egg — {pdb_id}")
        st.caption(f"{parsed.ca_coords.shape[0]} Cα atoms · PCA projection")
        st.image(data, use_container_width=True)
        st.download_button(
            label=f"Save {pdb_id}_pasqua.png",
            data=data,
            file_name=f"{pdb_id}_pasqua.png",
            mime="image/png",
            key=f"dl_{pdb_id}_{i}",
        )

    st.markdown("")
    st.caption("Structures courtesy of the [Protein Data Bank](https://www.rcsb.org/).")


if __name__ == "__main__":
    main()
