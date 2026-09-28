import re
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
VAR = HERE / "data" / "GAMT_variants_export_safe.csv"
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)

PROTEIN_LENGTH = 236

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "savefig.dpi": 600,
})

def variant_class(consequence):
    c = str(consequence)
    if "splice_donor" in c or "splice_acceptor" in c:
        return "Canonical splice"
    if any(t in c for t in ["stop_gained", "frameshift", "start_lost"]):
        return "Truncating"
    if "missense" in c:
        return "Missense"
    return "Other (splice region)"

STYLE = {"Truncating": (ORANGE, "s"), "Missense": (BLUE, "o"),
         "Canonical splice": (AQUA, "D"), "Other (splice region)": (AQUA, "^")}

def protein_position(cdna):
    m = re.match(r"c\.(-?\d+)", str(cdna))
    return max(1, math.ceil(int(m.group(1)) / 3)) if m else np.nan

v = pd.read_csv(VAR, encoding="utf-8-sig")
v = v[~v["variant (GRCh38)"].str.startswith("All variants")].copy()
v["class"] = v["consequence"].map(variant_class)
v["aa"] = v["cDNA (MANE ENST00000252288)"].map(protein_position)
v["n"] = pd.to_numeric(v["carriers"].astype(str).str.replace(",", ""), errors="coerce")
common, rare = v[v["n"].notna()], v[v["n"].isna()]

fig, ax = plt.subplots(figsize=(6.8, 2.8))
ax.add_patch(plt.Rectangle((1, -0.18), PROTEIN_LENGTH - 1, 0.36, color="#dedcd6", zorder=1))
ax.text(PROTEIN_LENGTH / 2, 0, f"GAMT protein ({PROTEIN_LENGTH} aa)", ha="center", va="center",
        fontsize=7, color=INK2, zorder=2)

stack = {}
for _, r in rare.sort_values("aa").iterrows():
    col, mk = STYLE[r["class"]]
    b = int(r["aa"] // 4)
    level = stack.get(b, 0)
    stack[b] = level + 1
    yy = 0.42 + 0.2 * level
    ax.vlines(r["aa"], 0.18, yy, color=GRID, lw=0.6, zorder=1)
    ax.scatter(r["aa"], yy, marker=mk, s=16, zorder=3, facecolor="white", edgecolor=col, linewidth=1.0)

top = 0.42 + 0.2 * (max(stack.values()) if stack else 0) + 0.35
for k, (_, r) in enumerate(common.sort_values("aa").iterrows()):
    col, mk = STYLE[r["class"]]
    yy = top + (0.45 if k % 2 else 0)
    ax.vlines(r["aa"], 0.18, yy, color=MUTED, lw=0.7, zorder=2)
    ax.scatter(r["aa"], yy, marker=mk, s=42, zorder=4, facecolor=col, edgecolor="white", linewidth=1.0)
    ax.text(r["aa"], yy + 0.16, r["cDNA (MANE ENST00000252288)"], ha="center", va="bottom", fontsize=6.3,
            color=INK, zorder=5, bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))

ax.set_xlim(-4, PROTEIN_LENGTH + 6)
ax.set_ylim(-0.4, top + 0.95)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel("Protein position (splice variants placed at the nearest codon)")

handles = [Line2D([], [], marker=STYLE[k][1], ls="", color=STYLE[k][0], markersize=5.5, label=k)
           for k in ["Truncating", "Missense", "Canonical splice", "Other (splice region)"]]
handles += [Line2D([], [], marker="o", ls="", markerfacecolor=MUTED, markeredgecolor=MUTED,
                   markersize=5.5, label="≥20 carriers"),
            Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=MUTED,
                   markersize=4.5, label="<20 carriers")]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=6, frameon=False,
          fontsize=6.8, handletextpad=0.2, columnspacing=0.9)
fig.tight_layout()

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure3_variant_positions.{ext}", bbox_inches="tight")
plt.close(fig)
print(f"Saved Figure3_variant_positions.png and .pdf to {OUT}")
