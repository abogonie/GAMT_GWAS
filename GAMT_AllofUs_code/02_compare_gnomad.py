import glob
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
DATA, OUT = HERE / "data", HERE / "results"
OUT.mkdir(exist_ok=True)

LOF_TERMS = ("stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost")
LOF_MAX_AF = 0.005
GROUPS = {
    "AFR": ["African/African American"],
    "AMR": ["Admixed American"],
    "EUR": ["European (non-Finnish)", "European (Finnish)", "Ashkenazi Jewish"],
    "Other combined": ["East Asian", "Middle Eastern", "South Asian", "Amish", "Remaining"],
}

def wilson(x, n, z=norm.ppf(0.975)):
    if n <= 0:
        return np.nan, np.nan, np.nan
    p = x / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), c + h

def one_in(x):
    return "—" if not x or np.isnan(x) else f"1 in {1 / x:,.0f}"

def summarize(sub, label, cohort="gnomAD v4"):
    rows = []
    groups = {"ALL": None, **GROUPS}
    for g, members in groups.items():
        if members is None:
            ac, an = sub["Allele Count"], sub["Allele Number"]
        else:
            ac = sum(sub[f"Allele Count {m}"].fillna(0) for m in members)
            an = sum(sub[f"Allele Number {m}"].fillna(0) for m in members)
        an = an[an > 0]
        n_alleles = float(an.median()) if len(an) else 0.0
        x = float(ac.sum())
        q, lo, hi = wilson(x, n_alleles)
        rows.append({"cohort": cohort, "set": label, "ancestry": g, "n_people": int(n_alleles // 2),
                     "n_variants": int((ac > 0).sum()), "allele_count": int(x),
                     "carrier (1 in)": one_in(1 - (1 - q) ** 2),
                     "carrier 95% CI": f"{one_in(1 - (1 - hi) ** 2)} to {one_in(1 - (1 - lo) ** 2)}",
                     "prevalence (1 in)": one_in(q * q),
                     "prevalence 95% CI": f"{one_in(hi * hi)} to {one_in(lo * lo)}",
                     "_q": q})
    return rows

gpath = sorted(glob.glob(str(DATA / "gnomAD*.csv")))[-1]
g = pd.read_csv(gpath)
aou_v = pd.read_csv(DATA / "GAMT_variants_export_safe.csv", encoding="utf-8-sig")
aou_s = pd.read_csv(DATA / "GAMT_summary_export_safe.csv", encoding="utf-8-sig")
aou_v = aou_v[~aou_v["variant (GRCh38)"].str.startswith("All variants")].copy()
print(f"gnomAD file: {Path(gpath).name} ({len(g):,} variants)")

passes = (g["Filters - exomes"].fillna("").eq("PASS") | g["Filters - genomes"].fillna("").eq("PASS"))
g = g[passes].copy()
g["vid"] = g["gnomAD ID"]
print(f"PASS variants: {len(g):,}; homozygotes among them: {int(g['Homozygote Count'].fillna(0).sum())}")

aou_plp = set(aou_v.loc[aou_v["in set"] == "P/LP", "variant (GRCh38)"])
aou_all = set(aou_v["variant (GRCh38)"])
found = g[g["vid"].isin(aou_all)]
print(f"A) All of Us qualifying variants found in gnomAD: {found['vid'].nunique()} of {len(aou_all)}")

cl = g["ClinVar Germline Classification"].fillna("").str.lower()
g["clinvar_plp"] = cl.str.contains("pathogenic") & ~cl.str.contains("conflicting|benign|uncertain")
vep = g["VEP Annotation"].fillna("").str.lower()
on_mane = g["Transcript"].fillna("").str.startswith("ENST00000252288")
g["lof"] = on_mane & vep.apply(lambda s: any(t in s for t in LOF_TERMS)) & (g["Allele Frequency"].fillna(0) < LOF_MAX_AF)
qual = g[g["clinvar_plp"] | g["lof"]].copy()
qual["in_allofus"] = qual["vid"].isin(aou_all)
print(f"B) gnomAD variants meeting the criteria: {len(qual)} "
      f"({int(qual['clinvar_plp'].sum())} P/LP); {int(qual['in_allofus'].sum())} also qualifying in All of Us")

rows = []
rows += summarize(found[found["vid"].isin(aou_plp)], "Same variants: ClinVar P/LP")
rows += summarize(found, "Same variants: P/LP + LoF")
rows += summarize(qual[qual["clinvar_plp"]], "Same criteria: ClinVar P/LP")
rows += summarize(qual, "Same criteria: P/LP + LoF")
res = pd.DataFrame(rows)

a = aou_s.copy()
a["cohort"] = "All of Us v9"
a["ancestry"] = a["ancestry"].str.replace(r"Other combined.*", "Other combined", regex=True)
a["carrier 95% CI"] = a["carrier 95% CI"].astype(str).str.replace(" – ", " to ", regex=False)
a["prevalence 95% CI"] = a["prevalence 95% CI"].astype(str).str.replace(" – ", " to ", regex=False)
cols = ["cohort", "set", "ancestry", "n_people", "n_variants", "allele_count",
        "carrier (1 in)", "carrier 95% CI", "prevalence (1 in)", "prevalence 95% CI"]
both = pd.concat([a[cols], res[cols]], ignore_index=True)
both.to_csv(OUT / "gnomad_vs_allofus_summary.csv", index=False, encoding="utf-8-sig")

common = aou_v[aou_v["carriers"] != "<20"].copy()
m = common.merge(g[["vid", "Allele Count", "Allele Number", "Allele Frequency", "Homozygote Count"]],
                 left_on="variant (GRCh38)", right_on="vid", how="left")
m["gnomAD carrier freq (1 in)"] = m["Allele Frequency"].map(lambda f: one_in(2 * f) if f == f else "not in gnomAD")
m = m.rename(columns={"carriers": "All of Us carriers", "carrier freq (1 in)": "All of Us carrier freq (1 in)",
                      "Allele Count": "gnomAD AC", "Allele Number": "gnomAD AN"})
m[["variant (GRCh38)", "cDNA (MANE ENST00000252288)", "protein", "All of Us carriers",
   "All of Us carrier freq (1 in)", "gnomAD AC", "gnomAD AN", "gnomAD carrier freq (1 in)"]] \
    .to_csv(OUT / "gnomad_common_variant_comparison.csv", index=False, encoding="utf-8-sig")

keep = ["vid", "HGVS Consequence", "Protein Consequence", "VEP Annotation", "ClinVar Germline Classification",
        "Allele Count", "Allele Number", "Allele Frequency", "Homozygote Count", "clinvar_plp", "lof", "in_allofus"]
qual.sort_values("Allele Count", ascending=False)[keep].to_csv(OUT / "gnomad_qualifying_variants.csv",
                                                               index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(both[cols].to_string(index=False))
print(f"\nSaved 3 files to {OUT}")
