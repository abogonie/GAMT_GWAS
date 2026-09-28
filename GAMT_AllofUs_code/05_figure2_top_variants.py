from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
VAR = HERE / "data" / "GAMT_variants_export_safe.csv"
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)

BLUE, ORANGE, AQUA, POOLED = "#2a78d6", "#eb6834", "#1baf7a", "#b9b8b2"
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

COLORS = {"Truncating": ORANGE, "Missense": BLUE, "Canonical splice": AQUA, "Other (splice region)": AQUA}

v = pd.read_csv(VAR, encoding="utf-8-sig")
is_pool = v["variant (GRCh38)"].str.startswith("All variants")
pool = v[is_pool].iloc[0]
v = v[~is_pool].copy()
v["class"] = v["consequence"].map(variant_class)
v["n"] = pd.to_numeric(v["carriers"].astype(str).str.replace(",", ""), errors="coerce")
common = v[v["n"].notna()].copy()
common["pct"] = common["% of pathogenic alleles"].str.rstrip("%").astype(float)

n_rare = int(pool["variant (GRCh38)"].split("(")[1].split()[0])
labels = [f"{r['cDNA (MANE ENST00000252288)']}  {r['protein']}" for _, r in common.iterrows()]
labels.append(f"{n_rare} rarer variants (<20 carriers each)")
values = list(common["pct"]) + [float(str(pool["% of pathogenic alleles"]).rstrip("%"))]
colors = [COLORS[c] for c in common["class"]] + [POOLED]

fig, ax = plt.subplots(figsize=(6.8, 3.0))
y = np.arange(len(values))[::-1]
ax.barh(y, values, color=colors, height=0.62, edgecolor="white", linewidth=2)
for yi, val in zip(y, values):
    ax.text(val + 0.5, yi, f"{val:.1f}%", va="center", fontsize=7.5, color=INK2)
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=7.5, color=INK)
ax.set_xlabel("Share of all pathogenic GAMT alleles (P/LP + LoF, n = 725)")
ax.set_xlim(0, 36)
ax.tick_params(axis="y", length=0)
ax.xaxis.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
shown = [k for k in ["Truncating", "Missense", "Other (splice region)", "Canonical splice"]
         if k in set(common["class"])]
ax.legend(handles=[Patch(color=COLORS[k], label=k) for k in shown] + [Patch(color=POOLED, label="Pooled rare variants")],
          loc="center right", frameon=False, fontsize=7.5)
fig.tight_layout()

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure2_top_variants.{ext}", bbox_inches="tight")
plt.close(fig)
print(f"Saved Figure2_top_variants.png and .pdf to {OUT}")
