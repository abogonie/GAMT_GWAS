from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, norm

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

s = pd.read_csv(OUT / "gnomad_vs_allofus_summary.csv", encoding="utf-8-sig")
s["n"] = s["n_people"].astype(str).str.replace(",", "").astype(int)


def wilson(x, n, z=norm.ppf(0.975)):
    p = x / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, c - h, c + h


def one_in(x):
    return f"1 in {1 / x:,.0f}"


PAIRS = [("Same variants: ClinVar P/LP", "ClinVar P/LP"),
         ("Same variants: P/LP + LoF", "P/LP + LoF"),
         ("Same criteria: ClinVar P/LP", "ClinVar P/LP"),
         ("Same criteria: P/LP + LoF", "P/LP + LoF")]

rows = []
for gset, aset in PAIRS:
    for anc in ["ALL", "AFR", "AMR", "EUR", "Other combined"]:
        a = s[(s.cohort == "All of Us v9") & (s.set == aset) & (s.ancestry == anc)].iloc[0]
        g = s[(s.cohort == "gnomAD v4") & (s.set == gset) & (s.ancestry == anc)].iloc[0]
        an_a, an_g = 2 * a.n, 2 * g.n
        table = [[a.allele_count, an_a - a.allele_count], [g.allele_count, an_g - g.allele_count]]
        p = chi2_contingency(table, correction=True)[1]
        q, lo, hi = wilson(a.allele_count + g.allele_count, an_a + an_g)
        rows.append({"comparison": "gnomAD " + gset[0].lower() + gset[1:], "ancestry": anc,
                     "people": f"{a.n + g.n:,}",
                     "AoU": a["carrier (1 in)"], "gnomAD": g["carrier (1 in)"], "het_p": f"{p:.2g}",
                     "pooled_carrier": one_in(1 - (1 - q) ** 2),
                     "pooled_CI": f"{one_in(1 - (1 - hi) ** 2)} to {one_in(1 - (1 - lo) ** 2)}",
                     "pooled_prev": f"1 in {1 / (q * q) / 1e6:.2f} M"})

res = pd.DataFrame(rows)
res.to_csv(OUT / "pooled_allofus_gnomad.csv", index=False, encoding="utf-8-sig")
print(res.to_string(index=False))
