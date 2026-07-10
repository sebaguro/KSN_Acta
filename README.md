# Kardashev's Conundrum / KSN — Acta Astronautica build

Complete, compile-ready LaTeX project for the reframed paper, in the Elsevier
`elsarticle` class. This revision fixes several issues found by diffing the
compiled PDF against your original full-article PDF (see "What changed" below).

## Upload to Overleaf

1. Keep this folder zipped (or re-zip it).
2. Overleaf: **New Project → Upload Project** → choose the zip.
3. It auto-detects `main.tex`. Compiler: **pdfLaTeX** (default).
4. **Recompile.** Builds out of the box — `elsarticle` + `elsarticle-harv` ship
   with Overleaf.

## File tree

```
KSN_Acta/
├── main.tex                  <- manuscript (elsarticle). Title/affiliation here.
├── body.tex                  <- all sections (the content).
├── references.bib            <- 37 references (author–year).
├── figures/                  <- your four figures (extracted from your PDF).
│   ├── fig1_energy_loglog.png   (Figure 1, log axis)
│   ├── fig2_energy_linear.png   (Figure 2, linear zoom)
│   ├── fig3_deltaP.png          (Figure 3, ΔP scatter)
│   └── fig4_ksn.png             (Figure 4, KSN B(t))
└── make_errorbar_figures.py  <- regenerates Fig 1 & Fig 2 WITH error bars
```

## The error bars — already done

The two energy figures in `figures/` are **already generated with error bars on the real
OWID data** (I pulled the public `owid-energy-data.csv` and ran the script for you). They
match your manuscript exactly: P(1965)=4.95 TW, P(2024)=20.16 TW, OLS slope
b=2.44×10¹¹ W/yr, R²=0.987. So the project compiles with error bars out of the box — no
action needed.

Each energy figure now shows **both** linear fits, with their slope values in the legend:
the unweighted OLS fit (b=2.44×10¹¹, your §2.3 headline) and the 1/σᵢ²-weighted
fit (b=2.56×10¹¹, §3.3) — so the figure agrees with every number in the paper. Figures 3
(ΔP) and 4 (KSN) deliberately carry **no** error bars; the reasons are now stated in their
captions (systematic errors cancel in the year-to-year difference for ΔP; sub-symbol on the
14-decade log axis for KSN).

If you want to **re-run or tune** it (e.g. change the adopted uncertainty model in
`fractional_sigma()`):

```bash
pip install numpy scipy matplotlib
python make_errorbar_figures.py --csv /path/to/owid-energy-data.csv
```

It writes both PNG and PDF versions of the two energy figures; the project uses the
**PNG** files (`fig1_energy_loglog.png`, `fig2_energy_linear.png`) so it renders out of the
box in Overleaf — drop them into `figures/` and recompile. It also prints the robustness
numbers. `--demo` runs on synthetic data for testing.

### Important framing for the error bars
OWID / the Energy Institute / the EIA publish **no** formal per-year uncertainties. The
bars are an **adopted** heteroscedastic model (≈10% mid-1960s → ≈1.5% post-2000), used for
a transparent sensitivity check — not measured errors. The captions and §3.3 say this; keep
it explicit. Tune it in `fractional_sigma()`.

## What changed in this revision (diff vs your original full article)

1. **Restored §3.6 "The Doomsday Clock".** It was dropped during reconstruction
   (along with the Sagan 1985 citation it carries). Now back in, after §3.5
   (Drake). ⚠ Its final sentence — that decentralised proof-of-work "may offer a
   partial technological foundation for the amity Sagan envisioned" — is the most
   advocacy-adjacent line in the paper and sits in mild tension with the
   Bitcoin-as-proxy framing of the reframe.
2. **Fixed figure order.** original has Fig 2 = linear-zoom, Fig 3 = ΔP. The
   first compile had them swapped (and referenced out of order). Now matches
   original: Fig 1 log, Fig 2 linear, Fig 3 ΔP, Fig 4 KSN.
3. **Fixed citation rendering.** Added the `authoryear` class option, so citations
   render as "Kardashev (1964)" etc. The first compile showed broken numeric
   forms — "Ćirković's [2018]", "Type II (16 and 27)", "attributed to 24".
   (If Acta's editor insists on numbered references, remove `authoryear` from the
   class options in `main.tex` and switch to `\bibliographystyle{elsarticle-num}`.)

All key numbers were checked against original and match (P(1965)=4.95 TW,
P(2024)=20.16 TW, b=2.44×10¹¹, r=2.01%, ΔWAIC=5.5, W=0.925, p=0.0014,
substitution 7.66 TW / 3.12×10¹¹ / 1.85% / 8.9×10⁴ H₀⁻¹, etc.).


