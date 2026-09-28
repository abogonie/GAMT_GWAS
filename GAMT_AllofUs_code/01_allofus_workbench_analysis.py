# GAMT carrier frequency analysis, All of Us CDR v9

# 1. Settings

import os, subprocess

PROJECT = None

VAT_FOLDER = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/vat/"

GENE = "GAMT"
CHROM, START, END = "chr19", 1_392_000, 1_406_600

LOF_MAX_AF = 0.005

MANE_TX = ("ENST00000252288", "NM_000156")

if PROJECT is None:
    PROJECT = os.environ.get("GOOGLE_CLOUD_PROJECT") or os.environ.get("GOOGLE_PROJECT")
if PROJECT is None:
    try:
        PROJECT = subprocess.run(["gcloud", "config", "get-value", "project"],
                                 capture_output=True, text=True).stdout.strip() or None
    except FileNotFoundError:
        pass
print("Billing project:", PROJECT)
assert PROJECT, "Could not detect the project. Set PROJECT above (Workspace > Overview shows the Google project ID)."

# 2. Helpers: ranged reads from Cloud Storage and a tabix-index reader

import io, gzip, struct, zlib

def _split(gs_path):
    b, _, name = gs_path[len("gs://"):].partition("/")
    return b, name

try:
    from google.cloud import storage
    _client = storage.Client(project=PROJECT)
    def read_range(gs_path, start, end):
        b, name = _split(gs_path)
        blob = _client.bucket(b, user_project=PROJECT).blob(name)
        return blob.download_as_bytes(start=start, end=end - 1)
    def list_folder(gs_folder):
        b, prefix = _split(gs_folder)
        return [f"gs://{b}/{x.name}" for x in
                _client.list_blobs(_client.bucket(b, user_project=PROJECT), prefix=prefix)]
    def object_size(gs_path):
        b, name = _split(gs_path)
        blob = _client.bucket(b, user_project=PROJECT).get_blob(name)
        return blob.size
    print("Using google-cloud-storage")
except ImportError:
    def read_range(gs_path, start, end):
        return subprocess.run(["gsutil", "-u", PROJECT, "cat", "-r", f"{start}-{end-1}", gs_path],
                              capture_output=True, check=True).stdout
    def list_folder(gs_folder):
        out = subprocess.run(["gsutil", "-u", PROJECT, "ls", gs_folder],
                             capture_output=True, text=True, check=True).stdout
        return [l.strip() for l in out.splitlines() if l.strip()]
    def object_size(gs_path):
        out = subprocess.run(["gsutil", "-u", PROJECT, "du", gs_path],
                             capture_output=True, text=True, check=True).stdout
        return int(out.split()[0])
    print("Using gsutil")

def bgzf_blocks(buf):
    pos = 0
    while pos + 18 <= len(buf):
        if buf[pos:pos + 2] != b"\x1f\x8b":
            raise ValueError(f"Not a BGZF block at offset {pos}")
        xlen = struct.unpack_from("<H", buf, pos + 10)[0]
        extra, bsize, p = buf[pos + 12: pos + 12 + xlen], None, 0
        while p < len(extra):
            si1, si2, slen = extra[p], extra[p + 1], struct.unpack_from("<H", extra, p + 2)[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack_from("<H", extra, p + 4)[0]
            p += 4 + slen
        total = bsize + 1
        if pos + total > len(buf):
            return
        yield pos, zlib.decompress(buf[pos:pos + total], 31)
        pos += total

def parse_tbi(raw):
    data = gzip.decompress(raw)
    assert data[:4] == b"TBI\x01", "Not a tabix index"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<8i", data, 4)
    names = data[36:36 + l_nm].split(b"\x00")[:n_ref]
    p, refs = 36 + l_nm, {}
    for i in range(n_ref):
        n_bin = struct.unpack_from("<i", data, p)[0]; p += 4
        for _ in range(n_bin):
            _, n_chunk = struct.unpack_from("<Ii", data, p); p += 8 + 16 * n_chunk
        n_intv = struct.unpack_from("<i", data, p)[0]; p += 4
        ioff = struct.unpack_from(f"<{n_intv}Q", data, p); p += 8 * n_intv
        refs[names[i].decode()] = ioff
    return dict(col_seq=col_seq, col_beg=col_beg, meta=chr(meta), skip=skip, refs=refs)

def fetch_region(vat, idx, chrom, start, end, chunk=4 << 20):
    ioff = idx["refs"][chrom]
    w = min(start >> 14, len(ioff) - 1)
    voff = next((v for v in ioff[w:] if v), 0)
    coff, uoff = voff >> 16, voff & 0xFFFF
    cs, cb = idx["col_seq"] - 1, idx["col_beg"] - 1

    out, tail, first, done = [], b"", True, False
    buf, base = b"", coff
    while not done:
        new = read_range(vat, base + len(buf), base + len(buf) + chunk)
        if not new:
            break
        buf += new
        used = 0
        for off, block in bgzf_blocks(buf):
            if first:
                block, first = block[uoff:], False
            text = tail + block
            lines = text.split(b"\n")
            tail = lines.pop()
            for ln in lines:
                f = ln.split(b"\t", max(cs, cb) + 1)
                c, ps = f[cs].decode(), int(f[cb])
                if c != chrom or ps > end:
                    done = True; break
                if ps >= start:
                    out.append(ln.decode())
            used = off + (struct.unpack_from("<H", buf, off + 16)[0] + 1)
            if done or len(block) == 0:
                done = True; break
        base, buf = base + used, buf[used:]
    return out

def read_header(vat, meta="#"):
    buf = read_range(vat, 0, 1 << 20)
    text = b"".join(b for _, b in bgzf_blocks(buf))
    first = text.split(b"\n", 1)[0].decode()
    return first.lstrip(meta).split("\t")

# 3. Locate the VAT and its index, read the GAMT region

files = list_folder(VAT_FOLDER)
for f in files: print(f)
tbi = next(f for f in files if f.endswith(".tbi"))
vat = next(f for f in files if not f.endswith(".tbi") and f + ".tbi" == tbi) if any(f + ".tbi" == tbi for f in files) \
      else next(f for f in files if not f.endswith(".tbi") and not f.endswith("/"))
print("\nVAT:  ", vat, "\nIndex:", tbi)

idx = parse_tbi(read_range(tbi, 0, object_size(tbi)))
if CHROM not in idx["refs"]:
    CHROM = CHROM.replace("chr", "") if CHROM.startswith("chr") else "chr" + CHROM
print("Contigs in index:", list(idx["refs"])[:5], "...  using", CHROM)

header = read_header(vat, idx["meta"])
lines = fetch_region(vat, idx, CHROM, START, END)
print(f"{len(header)} columns; {len(lines):,} VAT rows in {CHROM}:{START:,}-{END:,}")

# 4. Build the variant table

import pandas as pd, numpy as np

raw = pd.DataFrame([l.split("\t") for l in lines], columns=header)
print("Columns:", ", ".join(raw.columns))

def col(*cands, required=True):
    for c in cands:
        if c in raw.columns: return c
    for c in cands:
        hit = [x for x in raw.columns if c in x]
        if hit: return hit[0]
    if required: raise KeyError(f"None of {cands} found in VAT columns")
    return None

C_VID   = col("vid")
C_GENE  = col("gene_symbol")
C_CSQ   = col("consequence")
C_TX    = col("transcript", required=False)
C_CANON = col("is_canonical_transcript", required=False)
C_CLNV  = col("clinvar_classification", "clinvar_significance", "clinvar")
C_PCHG  = col("aa_change", required=False)
C_DNA   = col("dna_change_in_transcript", required=False)

gene = raw[raw[C_GENE] == GENE].copy()
print(f"{len(gene):,} rows annotated to {GENE}, {gene[C_VID].nunique():,} distinct variants")

print("GAMT transcripts in VAT:", gene[C_TX].str.split(".").str[0].value_counts().to_dict() if C_TX else "n/a")
gene["_rank"] = 2
if C_TX:    gene.loc[gene[C_TX].str.startswith(MANE_TX, na=False), "_rank"] = 0
if C_CANON: gene.loc[(gene["_rank"] > 0) & gene[C_CANON].str.lower().isin(["true", "1", "yes"]), "_rank"] = 1
n_mane = int((gene["_rank"] == 0).sum())
print(f"Rows on MANE transcript {MANE_TX[0]}: {n_mane:,}")
if n_mane == 0:
    print("WARNING: MANE transcript not found; LoF calls will be empty. Check the transcript list above.")
var = gene.sort_values("_rank").drop_duplicates(C_VID).copy()
var["on_mane"] = var["_rank"] == 0

POPS = [p for p in ["afr", "amr", "eas", "eur", "mid", "sas", "oth"] if f"gvs_{p}_ac" in raw.columns]
for p in ["all"] + POPS:
    for s in ["ac", "an", "af", "sc"]:
        c = f"gvs_{p}_{s}"
        if c in var.columns:
            var[c] = pd.to_numeric(var[c].replace("", np.nan), errors="coerce")
print("Ancestry groups in VAT:", POPS)

# 5. Classify variants

LOF_TERMS = ["stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost"]

cl = var[C_CLNV].fillna("").str.lower()
var["clinvar_plp"] = cl.str.contains("pathogenic") & ~cl.str.contains("conflicting|benign|uncertain")
csq = var[C_CSQ].fillna("").str.lower()
var["lof"] = var["on_mane"] & csq.apply(lambda s: any(t in s for t in LOF_TERMS)) \
             & (var["gvs_all_af"].fillna(0) < LOF_MAX_AF)

show = [c for c in [C_VID, C_DNA, C_PCHG, C_CSQ, C_CLNV, "gvs_all_ac", "gvs_all_an",
                    "gvs_all_af", "gvs_all_sc", "on_mane", "clinvar_plp", "lof"] if c]
qual = var[var.clinvar_plp | var.lof]
print(f"ClinVar P/LP: {var.clinvar_plp.sum()}   LoF: {var.lof.sum()}   Either: {len(qual)}")
qual[show].sort_values("gvs_all_af", ascending=False)

# 6. Carrier frequency and predicted prevalence

from scipy.stats import norm

def wilson(x, n, z=norm.ppf(0.975)):
    if n == 0: return (np.nan, np.nan, np.nan)
    p = x / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d; h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return p, max(0, c - h), c + h

def summarize(sub, label, groups=None):
    rows = []
    for p in list(groups) if groups else ["all"] + POPS:
        members = groups[p] if groups else [p]
        ac = sum(sub[f"gvs_{m}_ac"].fillna(0) for m in members)
        an = sum(sub[f"gvs_{m}_an"].fillna(0) for m in members)
        an = an[an > 0]
        sc = sum(sub[f"gvs_{m}_sc"].fillna(0) for m in members) if all(f"gvs_{m}_sc" in sub for m in members) \
             else pd.Series(dtype=float)
        N = float(an.median()) if len(an) else 0.0
        X = float(ac.sum())
        q, lo, hi = wilson(X, N)
        rows.append(dict(set=label, ancestry=p if groups else p.upper(), n_variants=int((ac > 0).sum()),
                         allele_count=int(X), median_AN=int(N), n_people=int(N // 2),
                         observed_carriers=int(sc.sum()) if len(sc) else np.nan,
                         q=q, carrier_freq=1-(1-q)**2, carrier_lo=1-(1-lo)**2, carrier_hi=1-(1-hi)**2,
                         prevalence=q*q, prev_lo=lo*lo, prev_hi=hi*hi))
    return pd.DataFrame(rows)

res = pd.concat([summarize(var[var.clinvar_plp], "ClinVar P/LP"),
                 summarize(var[var.clinvar_plp | var.lof], "P/LP + LoF")], ignore_index=True)

def one_in(x): return "—" if not x or np.isnan(x) else f"1 in {1/x:,.0f}"
pretty = res.assign(**{"carrier (1 in)": res.carrier_freq.map(one_in),
                       "carrier 95% CI": [f"{one_in(h)} – {one_in(l)}" for l, h in zip(res.carrier_lo, res.carrier_hi)],
                       "prevalence (1 in)": res.prevalence.map(one_in),
                       "prevalence 95% CI": [f"{one_in(h)} – {one_in(l)}" for l, h in zip(res.prev_lo, res.prev_hi)]})
pretty[["set", "ancestry", "n_people", "n_variants", "allele_count", "observed_carriers",
        "carrier (1 in)", "carrier 95% CI", "prevalence (1 in)", "prevalence 95% CI"]]

# 7. Export-safe summary

BIG = [p for p in ["afr", "amr", "eur"] if p in POPS]
SMALL = [p for p in POPS if p not in BIG]
groups = {"ALL": ["all"], **{p.upper(): [p] for p in BIG},
          "Other combined (" + "/".join(p.upper() for p in SMALL) + ")": SMALL}

res_safe = pd.concat([summarize(var[var.clinvar_plp], "ClinVar P/LP", groups),
                      summarize(var[var.clinvar_plp | var.lof], "P/LP + LoF", groups)], ignore_index=True)

def cnt(v): return "—" if pd.isna(v) else ("<20" if 0 < v < 20 else f"{int(v):,}")

def fmt(r):
    small = 0 < r.allele_count < 20 or (pd.notna(r.observed_carriers) and 0 < r.observed_carriers < 20)
    return pd.Series({
        "set": r.set, "ancestry": r.ancestry, "n_people": cnt(r.n_people), "n_variants": r.n_variants,
        "allele_count": cnt(r.allele_count), "observed_carriers": cnt(r.observed_carriers),
        "carrier (1 in)": "suppressed" if small else one_in(r.carrier_freq),
        "carrier 95% CI": "suppressed" if small else f"{one_in(r.carrier_hi)} – {one_in(r.carrier_lo)}",
        "prevalence (1 in)": "suppressed" if small else one_in(r.prevalence),
        "prevalence 95% CI": "suppressed" if small else f"{one_in(r.prev_hi)} – {one_in(r.prev_lo)}"})

safe = res_safe.apply(fmt, axis=1)

all_ok = True
for label, g in res_safe.groupby("set"):
    total = g.loc[g.ancestry == "ALL", "allele_count"].iat[0]
    parts = g.loc[g.ancestry != "ALL", "allele_count"]
    hidden = total - parts[(parts == 0) | (parts >= 20)].sum()
    ok = hidden == 0 or hidden >= 20
    all_ok &= ok
    print(f"{label}: total {total:,}, remainder not shown in any displayed row: "
          f"{'none' if hidden == 0 else ('>=20' if hidden >= 20 else '1-19')} -> {'OK' if ok else 'PROBLEM'}")
print("\nExport-safe:", "YES" if all_ok else "NO - do not share this table yet")

safe.to_csv("GAMT_summary_export_safe.csv", index=False)
safe

# 8. Export-safe variant table

def hgvs_c(x):  return x.split(":", 1)[-1] if isinstance(x, str) else ""
def hgvs_p(x):
    if not isinstance(x, str): return ""
    p = x.split(":", 1)[-1]
    i = p.find("p.")
    if i < 0: return ""
    q = p[i:].rstrip(")")
    return q + (")" if q.count("(") > q.count(")") else "")

GROUPS_V = {"AFR": ["afr"], "AMR": ["amr"], "EUR": ["eur"],
            "Other": [p for p in POPS if p not in ("afr", "amr", "eur")]}

qv = var[var.clinvar_plp | var.lof].copy()
qv["carriers"] = qv["gvs_all_sc"].fillna(qv["gvs_all_ac"]).fillna(0).astype(int)
qv = qv.sort_values("carriers", ascending=False)
total_alleles = int(qv["gvs_all_ac"].fillna(0).sum())

rows = []
for _, r in qv.iterrows():
    n = r.carriers
    common = n >= 20
    row = {
        "variant (GRCh38)": r[C_VID],
        "cDNA (MANE ENST00000252288)": hgvs_c(r[C_DNA]) if C_DNA else "",
        "protein": hgvs_p(r[C_PCHG]) if C_PCHG else "",
        "consequence": r[C_CSQ],
        "ClinVar": r[C_CLNV] if r[C_CLNV] else "not in ClinVar",
        "in set": "P/LP" if r.clinvar_plp else "LoF only",
        "carriers": f"{n:,}" if common else "<20",
        "allele freq": f"{r.gvs_all_af:.2e}" if common else "suppressed",
        "carrier freq (1 in)": one_in(2 * r.gvs_all_af) if common else "suppressed",
        "% of pathogenic alleles": f"{100 * r.gvs_all_ac / total_alleles:.1f}%" if common else "suppressed",
    }
    g = {k: int(sum(r.get(f"gvs_{m}_sc", r.get(f"gvs_{m}_ac", 0)) or 0 for m in ms)) for k, ms in GROUPS_V.items()}
    if common:
        masked = {k: (0 < v < 20) for k, v in g.items()}
        unassigned = n - sum(g.values())
        def hidden_total():
            return sum(v for k, v in g.items() if masked[k]) + max(unassigned, 0)
        while 0 < hidden_total() < 20:
            shown = [(v, k) for k, v in g.items() if not masked[k] and v > 0]
            if not shown:
                break
            masked[min(shown)[1]] = True
        if 0 < hidden_total() < 20:
            masked = {k: True for k in g}
        for k, v in g.items():
            row[f"carriers {k}"] = f"{v:,}" if not masked[k] else ("<20" if 0 < v < 20 else "masked")
    else:
        for k in g: row[f"carriers {k}"] = "—"
    rows.append(row)

vt = pd.DataFrame(rows)

rare = qv[qv.carriers < 20]
pool = int(rare.carriers.sum())
pool_row = {c: "" for c in vt.columns}
pool_row.update({"variant (GRCh38)": f"All variants with <20 carriers ({len(rare)} variants)",
                 "carriers": f"{pool:,}" if (pool == 0 or pool >= 20) else "<20",
                 "% of pathogenic alleles": f"{100 * rare.gvs_all_ac.sum() / total_alleles:.1f}%"
                                            if pool >= 20 else "suppressed"})
vt = pd.concat([vt, pd.DataFrame([pool_row])], ignore_index=True)

print(f"{len(qv)} qualifying variants: {int((qv.carriers >= 20).sum())} with >=20 carriers, {len(rare)} with <20.")
print("Pooled rare-variant carriers:", "OK to show" if (pool == 0 or pool >= 20) else "masked (<20)")
vt.to_csv("GAMT_variants_export_safe.csv", index=False)
vt

# 9. Individuals with two qualifying alleles

import re, numpy as np

SNPINDEL = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/"

def gs_list_recursive(prefix):
    try:
        return list_folder(prefix)
    except Exception:
        out = subprocess.run(["gsutil", "-u", PROJECT, "ls", prefix + "**"],
                             capture_output=True, text=True).stdout
        return [l.strip() for l in out.splitlines() if l.strip()]

def download(gs_path, local):
    try:
        b, name = _split(gs_path)
        _client.bucket(b, user_project=PROJECT).blob(name).download_to_filename(local)
    except NameError:
        subprocess.run(["gsutil", "-u", PROJECT, "cp", gs_path, local], check=True)
    return local

print("Searching the exome callset for chr19 genotype files (can take a minute)...")
files = gs_list_recursive(SNPINDEL + "exome")
if not files:
    files = gs_list_recursive(SNPINDEL)
bed_hits = [f for f in files if re.search(r"(^|[/._-])chr19\.bed$", f)]
mt_hits = sorted({f.split(".mt/")[0] + ".mt" for f in files if ".mt/" in f})
print("chr19 PLINK .bed:", bed_hits or "none")
print("Hail MatrixTables:", mt_hits or "none")

def key_of(vid):
    c, p, r, a = vid.split("-", 3)
    return int(p), r, a
QKEYS = {key_of(v): v for v in qv[C_VID]}
print(f"{len(QKEYS)} qualifying variants to genotype")

def mask_n(n): return "0" if n == 0 else ("<20" if n < 20 else f"{n:,}")

carriers = {}

if bed_hits:
    BED = bed_hits[0]; BIM = BED[:-4] + ".bim"; FAM = BED[:-4] + ".fam"
    print("Using", BED)
    fam = pd.read_csv(download(FAM, "chr19.fam"), sep=r"\s+", header=None, usecols=[1], dtype=str)
    N = len(fam); B = (N + 3) // 4
    bim_local = download(BIM, "chr19.bim")
    bim_rows = []
    for ch in pd.read_csv(bim_local, sep="\t", header=None, usecols=[1, 3, 4, 5],
                          names=["id", "pos", "a1", "a2"], dtype={"pos": int, "a1": str, "a2": str},
                          chunksize=1_000_000):
        ch = ch.reset_index().rename(columns={"index": "idx"})
        bim_rows.append(ch[(ch.pos >= START) & (ch.pos <= END)])
    bim = pd.concat(bim_rows)
    assert read_range(BED, 0, 3) == b"\x6c\x1b\x01", "Not a SNP-major PLINK .bed"

    EXPECT = dict(zip(var[C_VID], var["gvs_all_sc"].fillna(var["gvs_all_ac"])))
    best = {}
    rejected = []
    for _, b in bim.iterrows():
        cands = [(b.a1, b.a2), (b.a2, b.a1)]
        cands = [(r, a) for r, a in cands if (b.pos, r, a) in QKEYS]
        if not cands:
            continue
        raw = np.frombuffer(read_range(BED, 3 + b.idx * B, 3 + (b.idx + 1) * B), dtype=np.uint8)
        codes = np.stack([(raw >> s) & 3 for s in (0, 2, 4, 6)], axis=1).reshape(-1)[:N]
        for ref, alt in cands:
            vid = QKEYS[(b.pos, ref, alt)]
            lut = np.array([2, -1, 1, 0]) if b.a1 == alt else np.array([0, -1, 1, 2])
            dos = lut[codes]
            idx = np.nonzero(dos > 0)[0]
            exp = float(EXPECT.get(vid, np.nan))
            diff = abs(len(idx) - exp) if exp == exp else 0
            if exp == exp and len(idx) > 2 * exp + 10:
                rejected.append(vid)
                continue
            if vid not in best or diff < best[vid][0]:
                best[vid] = (diff, idx, dos[idx])
    carriers = {v: (i, d) for v, (_, i, d) in best.items()}
    matched = len(carriers)
    rej = sorted(set(rejected) - set(carriers))
    if rej:
        print("Could not match these variants to a plausible genotype row (excluded):", ", ".join(rej))
    print(f"Genotyped {matched} of {len(QKEYS)} qualifying variants in {N:,} samples")

elif mt_hits:
    import hail as hl
    hl.init(default_reference="GRCh38", gcs_requester_pays_configuration=PROJECT, idempotent=True)
    mt = hl.read_matrix_table(mt_hits[0])
    mt = hl.filter_intervals(mt, [hl.locus_interval(CHROM, START, END, reference_genome="GRCh38")])
    if "was_split" not in mt.row:
        mt = hl.split_multi_hts(mt)
    keys = hl.literal(set(QKEYS))
    mt = mt.filter_rows(keys.contains((mt.locus.position, mt.alleles[0], mt.alleles[1])))
    mt = mt.annotate_rows(vid=hl.str("19-") + hl.str(mt.locus.position) + "-" + mt.alleles[0] + "-" + mt.alleles[1])
    ent = mt.filter_entries(mt.GT.is_non_ref()).entries()
    ent = ent.select(vid=ent.vid, s=ent.s, d=ent.GT.n_alt_alleles()).collect()
    N = mt.count_cols()
    df_e = pd.DataFrame([(e.vid, e.s, e.d) for e in ent], columns=["vid", "s", "d"])
    s_codes = {s: i for i, s in enumerate(df_e.s.unique())}
    for v, g in df_e.groupby("vid"):
        carriers[v] = (g.s.map(s_codes).to_numpy(), g.d.to_numpy())
    print(f"Genotyped {len(carriers)} of {len(QKEYS)} qualifying variants in {N:,} samples (Hail)")
else:
    raise SystemExit("No chr19 PLINK file or Hail MatrixTable found. Paste the 'chr19 PLINK' and "
                     "'Hail MatrixTables' lines above (paths only) so the code can be pointed at the right files.")

agree = sum(1 for v, (i, d) in carriers.items()
            if abs(len(i) - int(var.loc[var[C_VID] == v, "gvs_all_sc"].fillna(-1).iat[0])) <= 2)
print(f"Carrier counts agree with the VAT (within 2) for {agree} of {len(carriers)} variants")

from collections import defaultdict
per_sample = defaultdict(list)
for v, (idx, dos) in carriers.items():
    for i, d in zip(idx, dos):
        per_sample[int(i)].append((v, int(d)))

hom = {s for s, L in per_sample.items() if any(d == 2 for _, d in L)}
multi = {s for s, L in per_sample.items() if len({v for v, _ in L}) >= 2}
two_plus = hom | multi
plp_ids = set(qv.loc[qv.clinvar_plp, C_VID])
multi_plp = {s for s in multi if len({v for v, _ in per_sample[s]} & plp_ids) >= 2}

print("People carrying >=1 qualifying allele: ", mask_n(len(per_sample)))
print("Homozygous for a qualifying variant:   ", mask_n(len(hom)))
print("Two different qualifying variants:     ", mask_n(len(multi)))
print("   ...both ClinVar P/LP:               ", mask_n(len(multi_plp)))
print("Any person with >=2 qualifying alleles:", mask_n(len(two_plus)))

if multi:
    pairs = defaultdict(int)
    for s in multi:
        vs = sorted({v for v, _ in per_sample[s]})
        pairs[tuple(vs)] += 1
    print("\nVariant combinations seen (counts masked):")
    for vs, n in pairs.items():
        print("  ", " + ".join(vs), "->", mask_n(n))
    print("\nNext: check phase (cis vs trans) in the v9 phasing data before calling anyone affected.")
