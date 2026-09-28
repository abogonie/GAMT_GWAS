import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SUM = HERE / "data" / "GAMT_summary_export_safe.csv"
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "savefig.dpi": 600,
})

def one_in(text):
    m = re.search(r"1 in ([\d,]+)", str(text))
    return float(m.group(1).replace(",", "")) if m else np.nan

s = pd.read_csv(SUM, encoding="utf-8-sig")
s["anc"] = s["ancestry"].str.replace(r"Other combined.*", "Other combined*", regex=True)
s["est"] = 1e4 / s["carrier (1 in)"].map(one_in)
ci = s["carrier 95% CI"].str.replace(" – ", " to ", regex=False).str.split(" to ", expand=True)
s["hi"] = 1e4 / ci[0].map(one_in)
s["lo"] = 1e4 / ci[1].map(one_in)

order = ["ALL", "EUR", "AFR", "AMR", "Other combined*"]
labels = {"ALL": "All participants", "EUR": "European", "AFR": "African",
          "AMR": "Admixed American", "Other combined*": "Other groups combined*"}
npeople = s[s["set"] == "ClinVar P/LP"].set_index("anc")["n_people"]

fig, ax = plt.subplots(figsize=(6.8, 3.0))
series = [("ClinVar P/LP", BLUE, "o", -0.14, "ClinVar pathogenic / likely pathogenic"),
          ("P/LP + LoF", ORANGE, "s", 0.14, "P/LP + predicted loss-of-function")]
for name, col, marker, offset, legend_label in series:
    d = s[s["set"] == name].set_index("anc").loc[order]
    y = np.arange(len(order))[::-1] + offset
    ax.hlines(y, d["lo"], d["hi"], color=col, lw=2, zorder=2)
    ax.scatter(d["est"], y, s=34, marker=marker, color=col, edgecolor="white", linewidth=1.2,
               zorder=3, label=legend_label)
    for yi, (_, r) in zip(y, d.iterrows()):
        ax.text(r["hi"] + 0.25, yi, f"1 in {1e4 / r['est']:,.0f}", va="center", fontsize=7.5, color=INK2)

ax.set_yticks(np.arange(len(order))[::-1])
ax.set_yticklabels([f"{labels[a]}\n(n = {npeople[a]})" for a in order], color=INK)
ax.axhline(len(order) - 1.5, color=GRID, lw=0.8)
ax.set_xlabel("Carriers per 10,000 people (95% CI)")
ax.set_xlim(0, 22.5)
ax.xaxis.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
ax.tick_params(axis="y", length=0)
ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, frameon=False, fontsize=7.5,
          handletextpad=0.4)
fig.text(0.01, 0.01, "*East Asian, Middle Eastern, South Asian and other/unassigned genetic ancestry groups, "
         "combined to comply with the All of Us small-count policy.", fontsize=6.5, color=MUTED)
fig.tight_layout(rect=(0, 0.04, 1, 1))

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Figure1_carrier_frequency_by_ancestry.{ext}", bbox_inches="tight")
plt.close(fig)
print(f"Saved Figure1_carrier_frequency_by_ancestry.png and .pdf to {OUT}")
