#!/usr/bin/env python3
"""
Phase 1.3b -- rebuild Figure 3 from the verified per_seed_metrics.csv, replacing
the retracted screenshot-transcribed version.

Recomputes, from raw per-seed data (not from any prior hardcoded figure script):
  - pooled candidate ipTM values per system (for the boxplot + strip)
  - 5-number boxplot summary per system (candidates, and negative control)
  - native-HIP mean and 95% CI band per system (Eq. 2, t_crit=2.776, n=5)
  - jittered scatter coordinates (seeded RNG, same visual approach as before)

Outputs a fresh fig_distribution.tex (TikZ/pgfplots) and compiles it with
tectonic. Does not touch the manuscript.
"""
import csv, math, subprocess
from collections import defaultdict
import numpy as np
import os as _os
_REPO = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

DATA = _os.path.join(_REPO, "data")
OUT_TEX_DIR = _os.path.join(_REPO, "results", "figure3_rebuild")
import os
os.makedirs(OUT_TEX_DIR, exist_ok=True)

T_CRIT = 2.776

rows = list(csv.DictReader(open(f"{DATA}/per_seed_metrics.csv")))
iptm_rows = [r for r in rows if r["metric"] == "TCR_pMHC_ipTM" and r["provenance"] == "WORKBOOK"]

grouped = defaultdict(dict)
for r in iptm_rows:
    key = (r["system"], r["role"], r["peptide_label"], r["sequence"])
    grouped[key][int(r["seed"])] = float(r["value"])

SYSTEMS = ["HIP6/A2.11/DQ8", "HIP11/8.E3/DQ8", "HIP11/8.E3/DQ8-trans"]

def vals_for(system, role):
    out = []
    for (s, r, label, seq), seeds in grouped.items():
        if s == system and r == role:
            out.append([seeds[i] for i in range(1, 6)])
    return out

def box_stats(pooled_vals):
    a = np.array(pooled_vals)
    q1, med, q3 = np.percentile(a, [25, 50, 75])
    iqr = q3 - q1
    lo = a[a >= q1 - 1.5 * iqr].min()
    hi = a[a <= q3 + 1.5 * iqr].max()
    return dict(lower=round(float(lo), 3), q1=round(float(q1), 3), median=round(float(med), 3),
                q3=round(float(q3), 3), upper=round(float(hi), 3))

def native_ci(native_vals_5):
    a = np.array(native_vals_5)
    m = a.mean()
    sd = a.std(ddof=1)
    hw = T_CRIT * sd / math.sqrt(5)
    return m, m - hw, m + hw

rng = np.random.default_rng(11)  # same seed as the original figure, for visual continuity

def jitter_coords(pooled, x0, width=0.16):
    j = rng.uniform(-width, width, size=len(pooled))
    return [(round(x0 + jj, 4), v) for jj, v in zip(j, pooled)]

panel_data = {}
for system in SYSTEMS:
    cand_vals = vals_for(system, "candidate")
    pooled = [v for c in cand_vals for v in c]
    ctrl_vals = vals_for(system, "negative_control")[0]
    native_vals = vals_for(system, "native_HIP")[0]

    bstats_cand = box_stats(pooled)
    bstats_ctrl = box_stats(ctrl_vals)
    n_mean, n_lo, n_hi = native_ci(native_vals)
    cand_pts = jitter_coords(pooled, 0.35)
    ctrl_pts = jitter_coords(ctrl_vals, 1.35, width=0.08)
    panel_data[system] = dict(
        n=len(cand_vals), bstats_cand=bstats_cand, bstats_ctrl=bstats_ctrl,
        native_mean=n_mean, native_lo=n_lo, native_hi=n_hi,
        cand_pts=cand_pts, ctrl_pts=ctrl_pts, ctrl_mean=round(float(np.mean(ctrl_vals)), 3),
    )
    print(f"{system}: n_candidates={len(cand_vals)} native=[{n_lo:.3f},{n_mean:.3f},{n_hi:.3f}] "
          f"cand_box={bstats_cand} ctrl_box={bstats_ctrl} ctrl_mean={panel_data[system]['ctrl_mean']}")

def coords_str(pts):
    return " ".join(f"({x},{y})" for x, y in pts)

TEMPLATE = r"""\documentclass[tikz,border=2pt]{standalone}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\usepgfplotslibrary{groupplots}
\usepgfplotslibrary{statistics}
\definecolor{seriesblue}{HTML}{2a78d6}
\definecolor{seriesred}{HTML}{c0392b}
\definecolor{bandblue}{HTML}{cfe0f5}

\begin{document}
\begin{tikzpicture}[font=\sffamily]
\begin{groupplot}[
  group style={group size=3 by 1, horizontal sep=14pt, group name=distgrid},
  width=5.6cm, height=8cm,
  ymin=0.10, ymax=0.98,
  xmin=-0.15, xmax=1.75,
  axis lines=left,
  tick align=outside,
  ytick={0.2,0.3,...,0.9},
  xtick={0.35,1.35},
  xticklabels={Candidates, Negative\\control},
  x tick label style={align=center, font=\scriptsize},
  y tick label style={font=\scriptsize},
  tick style={line width=0.5pt, black},
  axis line style={line width=0.6pt},
  ytick pos=left,
]

%%PANEL1%%
%%PANEL2%%
%%PANEL3%%

\end{groupplot}

\begin{scope}[shift={(distgrid c3r1.north east)}, xshift=18pt]
\node[circle, fill=seriesblue, fill opacity=0.55, minimum size=3pt, inner sep=0, anchor=center] (l1) at (0.35cm,-0.3cm) {};
\node[font=\scriptsize, anchor=west] at (0.5cm,-0.3cm) {Individual seed value};

\draw[fill=white, draw=black, line width=0.6pt] (0.28cm,-0.85cm) rectangle (0.42cm,-1.15cm);
\draw[black, line width=0.6pt] (0.28cm,-1.0cm) -- (0.42cm,-1.0cm);
\node[font=\scriptsize, anchor=west] at (0.5cm,-1.0cm) {Candidate seeds (boxplot)};

\draw[fill=seriesred, fill opacity=0.25, draw=seriesred, line width=0.6pt] (0.28cm,-1.55cm) rectangle (0.42cm,-1.85cm);
\draw[seriesred, line width=0.6pt] (0.28cm,-1.7cm) -- (0.42cm,-1.7cm);
\node[font=\scriptsize, anchor=west] at (0.5cm,-1.7cm) {Negative control (boxplot)};

\node[rectangle, rotate=45, fill=seriesred, draw=seriesred, minimum size=3.2pt, inner sep=0, anchor=center] at (0.35cm,-2.35cm) {};
\node[font=\scriptsize, anchor=west] at (0.5cm,-2.35cm) {Control seed value};

\draw[black, line width=1pt] (0.28cm,-2.9cm) -- (0.42cm,-2.9cm);
\node[font=\scriptsize, anchor=west] at (0.5cm,-2.9cm) {Native HIP mean};

\draw[fill=bandblue, draw=none, opacity=0.6] (0.28cm,-3.25cm) rectangle (0.42cm,-3.55cm);
\node[font=\scriptsize, anchor=west, align=left] at (0.5cm,-3.4cm) {Native HIP 95\% CI};
\end{scope}

\end{tikzpicture}
\end{document}
"""

panel_titles = {"HIP6/A2.11/DQ8": "HIP6/A2.11/DQ8", "HIP11/8.E3/DQ8": "HIP11/8.E3/DQ8",
                 "HIP11/8.E3/DQ8-trans": "HIP11/8.E3/DQ8-trans"}

panels_tex = []
for idx, system in enumerate(SYSTEMS, start=1):
    d = panel_data[system]
    bc, bo = d["bstats_cand"], d["bstats_ctrl"]
    opts = "title={\\small %s}" % panel_titles[system]
    if idx == 1:
        opts += r", ylabel={TCR--pMHC ipTM (individual seed values)}, ylabel style={font=\scriptsize}"
    else:
        opts += ", yticklabels={}"
    tex = f"\\nextgroupplot[{opts}]\n"
    tex += f"\\draw[fill=bandblue, draw=none, opacity=0.6] (axis cs:-0.15,{d['native_lo']:.4f}) rectangle (axis cs:1.75,{d['native_hi']:.4f});\n"
    tex += f"\\draw[black, line width=1pt] (axis cs:-0.15,{d['native_mean']:.4f}) -- (axis cs:1.75,{d['native_mean']:.4f});\n"
    tex += (f"\\addplot+[boxplot prepared={{lower whisker={bc['lower']}, lower quartile={bc['q1']}, "
            f"median={bc['median']}, upper quartile={bc['q3']}, upper whisker={bc['upper']}, draw position=0.35}}, "
            "boxplot/box extend=0.42, fill=white, draw=black, line width=0.6pt, mark=none] coordinates {};\n")
    tex += (f"\\addplot+[boxplot prepared={{lower whisker={bo['lower']}, lower quartile={bo['q1']}, "
            f"median={bo['median']}, upper quartile={bo['q3']}, upper whisker={bo['upper']}, draw position=1.35}}, "
            "boxplot/box extend=0.26, fill=seriesred, fill opacity=0.25, draw=seriesred, line width=0.6pt, mark=none] coordinates {};\n")
    tex += ("\\addplot[only marks, mark=*, mark size=0.9pt, draw=none, fill=seriesblue, fill opacity=0.55] coordinates {\n"
            f"{coords_str(d['cand_pts'])}\n}};\n")
    tex += ("\\addplot[only marks, mark=diamond*, mark size=1.6pt, draw=seriesred, line width=0.4pt, fill=seriesred] coordinates {\n"
            f"{coords_str(d['ctrl_pts'])}\n}};\n")
    tex += f"\\node[font=\\scriptsize, anchor=south, fill=white, inner sep=1pt] at (axis cs:1.35,{d['ctrl_mean']+0.06:.3f}) {{${d['ctrl_mean']}$}};\n"
    panels_tex.append(tex)

final_tex = TEMPLATE.replace("%%PANEL1%%", panels_tex[0]).replace("%%PANEL2%%", panels_tex[1]).replace("%%PANEL3%%", panels_tex[2])

tex_path = f"{OUT_TEX_DIR}/fig_distribution_rebuilt.tex"
with open(tex_path, "w") as f:
    f.write(final_tex)
print(f"\nWrote {tex_path}")

r = subprocess.run(["tectonic", "fig_distribution_rebuilt.tex"], cwd=OUT_TEX_DIR, capture_output=True, text=True)
print(r.stdout[-1500:])
print(r.stderr[-1500:])
print("Return code:", r.returncode)
