#!/usr/bin/env python3
"""
Composition-preserving shuffled-proteome control for the motif screen.

This is the control listed as deferred work in Section 3.13. It requires no
structural modelling -- only sequence scanning -- so it can be run now.

Method
------
1. Download the reviewed (Swiss-Prot) sequences for the taxa in Table 10 from
   UniProt. Herpesviridae (TaxID 10292) no longer resolves; the manuscript
   already identifies Orthoherpesviridae (3044472) as the current family, which
   is used here.
2. Count matches of the three structural-screen motifs (Table 5) in the real
   sequences, counting overlapping matches at every eligible start position.
3. Build a null by shuffling each sequence independently. Per-sequence shuffling
   preserves each protein's own amino-acid composition and length exactly, and
   destroys only the ordering -- which is the property the motif is supposed to
   detect.
4. Repeat, and compare the observed count with the null distribution.
5. Additionally run a scrambled-motif control: permute the order of the
   constrained classes within each motif, keeping the composition of constraints
   identical, and compare yields.

Important caveat, stated in the output: this is a CONTEMPORARY reconstruction of
the search space. The original ScanProsite run used the database as it stood
then, and taxonomy has since drifted. The observed count here is not expected to
equal the published 1,092, and the comparison that matters is
observed-vs-shuffled within this reconstruction, not either against 1,092.
"""
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
import re, sys, json, random, subprocess, statistics
from collections import defaultdict

RES = _os.path.join(_REPO, "results")
CACHE = _os.path.join(_REPO, "data", "proteome_cache")
_os.makedirs(CACHE, exist_ok=True)

TAXA = [("Mumps virus", 2560602), ("Rubella virus", 11041), ("SARS-CoV-2", 2697049),
        ("Enterovirus", 12059), ("Rotavirus", 10912), ("Orthoherpesviridae", 3044472),
        ("Bacteroides dorei", 357276), ("Clostridium sporogenes", 1505),
        ("Klebsiella oxytoca", 571), ("Parabacteroides distasonis", 823)]

# Table 5, structural screen. Written literally from the printed ScanProsite patterns.
MOTIFS = {
    "1": ".{5}[AGST][LIV][ED][AGST][ED].{4}",
    "2": ".{5}[AGST][LIV].[AGST][ED].[LIV].{2}",
    "5": ".{3}[ED][LIV].{3}.[AGST][LIV][ED].{3}",
}
# the constrained classes of each motif, in order, for the scrambled-motif control
CLASSES = {
    "1": [("x",5),("AGST",1),("LIV",1),("ED",1),("AGST",1),("ED",1),("x",4)],
    "2": [("x",5),("AGST",1),("LIV",1),("x",1),("AGST",1),("ED",1),("x",1),("LIV",1),("x",2)],
    "5": [("x",3),("ED",1),("LIV",1),("x",3),("x",1),("AGST",1),("LIV",1),("ED",1),("x",3)],
}

def fetch(tid, name):
    path = _os.path.join(CACHE, f"{tid}.fasta")
    if _os.path.exists(path) and _os.path.getsize(path) > 0:
        return path
    url = (f"https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+"
           f"taxonomy_id:{tid}&format=fasta")
    subprocess.run(["curl", "-s", "--max-time", "180", url, "-o", path], check=False)
    return path

def read_fasta(path):
    seqs, cur = [], []
    if not _os.path.exists(path): return seqs
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            else:
                cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return [s for s in seqs if s]

def count_overlapping(seqs, pattern):
    rx = re.compile(f"(?=({pattern}))")
    return sum(len(rx.findall(s)) for s in seqs)

def eligible_starts(seqs, mlen):
    return sum(max(0, len(s) - mlen + 1) for s in seqs)

# ------------------------------------------------------------------ gather
print("Downloading reviewed sequences (cached after first run)...", flush=True)
all_seqs, per_taxon = [], []
for name, tid in TAXA:
    p = fetch(tid, name)
    s = read_fasta(p)
    res = sum(len(x) for x in s)
    per_taxon.append((name, tid, len(s), res))
    all_seqs.extend(s)
    print(f"  {name:26s} tax={tid:<8} seqs={len(s):<6} residues={res:,}", flush=True)

total_res = sum(len(s) for s in all_seqs)
print(f"\nReconstructed search space: {len(all_seqs):,} sequences, {total_res:,} residues")
if total_res == 0:
    print("No sequences retrieved; aborting."); sys.exit(1)

# ------------------------------------------------------------------ observed
print("\n" + "="*74); print("OBSERVED MATCHES (contemporary reconstruction)"); print("="*74)
obs = {}
for mid, pat in MOTIFS.items():
    mlen = len(re.sub(r"\[[^\]]*\]", "X", pat).replace(".", "X").replace("{","").replace("}",""))
    n = count_overlapping(all_seqs, pat)
    obs[mid] = n
    print(f"  motif {mid}: {n:,} matches")
obs_total = sum(obs.values())
print(f"  TOTAL across the three structural motifs: {obs_total:,}")
print(f"  (published original-run Tier 1 yield, different database snapshot: 1,092)")

# ------------------------------------------------------------------ null
N_SHUFFLE = int(_os.environ.get("N_SHUFFLE", "30"))
print("\n" + "="*74)
print(f"COMPOSITION-PRESERVING SHUFFLED-PROTEOME NULL ({N_SHUFFLE} replicates)")
print("="*74, flush=True)
rng = random.Random(0)
null_tot, null_by = [], defaultdict(list)
for i in range(N_SHUFFLE):
    sh = []
    for s in all_seqs:
        l = list(s); rng.shuffle(l); sh.append("".join(l))
    t = 0
    for mid, pat in MOTIFS.items():
        c = count_overlapping(sh, pat); null_by[mid].append(c); t += c
    null_tot.append(t)
    if (i+1) % 10 == 0: print(f"  ...{i+1}/{N_SHUFFLE}", flush=True)

mu = statistics.mean(null_tot); sd = statistics.stdev(null_tot) if len(null_tot) > 1 else 0.0
print(f"\n  null total: mean {mu:,.0f}   sd {sd:,.0f}   range {min(null_tot):,}-{max(null_tot):,}")
print(f"  observed  : {obs_total:,}")
if sd > 0:
    print(f"  z = (observed - null_mean)/null_sd = {(obs_total-mu)/sd:+.2f}")
ge = sum(1 for v in null_tot if v >= obs_total)
print(f"  empirical one-sided p (null >= observed): {(ge+1)/(N_SHUFFLE+1):.4f}")
print(f"  ratio observed/null_mean: {obs_total/mu:.3f}" if mu else "")
print("\n  per motif:")
for mid in MOTIFS:
    m2 = statistics.mean(null_by[mid])
    print(f"    motif {mid}: observed {obs[mid]:,}  null mean {m2:,.0f}  ratio {obs[mid]/m2:.3f}")

# ------------------------------------------------------------------ scrambled motif
print("\n" + "="*74); print("SCRAMBLED-MOTIF CONTROL (real sequences, permuted constraint order)")
print("="*74, flush=True)
def build(parts):
    out = []
    for cls, n in parts:
        tok = "." if cls == "x" else f"[{cls}]"
        out.append(tok * n if n <= 2 else (tok + "{" + str(n) + "}" if cls == "x" else tok * n))
    return "".join(out)
for mid, parts in CLASSES.items():
    yields = []
    for _ in range(12):
        p = parts[:]; rng.shuffle(p)
        yields.append(count_overlapping(all_seqs, build(p)))
    print(f"  motif {mid}: real {obs[mid]:,}   scrambled-order median {statistics.median(yields):,.0f} "
          f"(range {min(yields):,}-{max(yields):,})")

out = dict(reconstruction=dict(sequences=len(all_seqs), residues=total_res,
                               per_taxon=[dict(name=a,taxid=b,seqs=c,residues=d) for a,b,c,d in per_taxon]),
           observed=obs, observed_total=obs_total,
           null=dict(replicates=N_SHUFFLE, mean=mu, sd=sd, min=min(null_tot), max=max(null_tot),
                     empirical_p=(ge+1)/(N_SHUFFLE+1)))
with open(_os.path.join(RES, "motif_null_control.json"), "w") as f:
    json.dump(out, f, indent=2)
print(f"\nWrote {RES}/motif_null_control.json")
print("\nCAVEAT: contemporary reconstruction, not the original search snapshot.")
