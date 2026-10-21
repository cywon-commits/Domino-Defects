# Electrostatics of forced defects in dimer coverings — code

Code and data behind the paper

> C. Won, *Electrostatics of forced defects in dimer coverings: corner charges, image charges, and boundary chemical potentials* (manuscript in `paper/`).

All numbers in the paper come from **exact** computations (transfer matrices, Kasteleyn determinants, LGV determinants, certified ball arithmetic), not Monte Carlo. Every script is small and self-contained; results are printed to stdout and/or written as JSON next to the script.

## Requirements

```
python >= 3.10
pip install -r requirements.txt      # numpy scipy sympy mpmath matplotlib python-flint
gcc                                  # for the C transfer-matrix codes in dp/
```

`python-flint` provides exact integer/rational linear algebra (FLINT) and certified ball arithmetic (Arb). It is essential for the Aztec-diamond, half-plane and quadrant computations, where double precision fails.

## Layout

| Folder | Contents | Paper |
|---|---|---|
| `dp/` | Broken-profile transfer-matrix counts in C with `__int128` (`cnt.c`, `cntf.c`, `full.c`, `one.c`, `mut.c`); Pfaffian control for the triangular lattice (`tri.c`); Python reference DP (`mut.py`) | Sec. II, VII |
| `lattice/` | Kasteleyn machinery: sparse log-counts on boards (`kast.py`) and general regions (`kreg.py`, incl. `aztec(n)`), rectangles with strings (`rect.py`), dense-inverse board with bordered determinants (`mutil.py`), Aztec regions (`azr.py`), Coulomb-gas maps and densities (`theory.py`) | Sec. II–III |
| `coulomb/` | Runs and analyses for single corner, mutilated board, odd boards, rectangles, holes, boundary/bulk pairs; `export.py` builds `data/figdata.json`, `figfix_data.py` adds the full-resolution maps of Figs. 3, 4 and 6(a) | Sec. III, Figs. 2–6 |
| `exact/` | Exact plane Green function `G = a + b/π` (`exactG.py`); certified determinants for bulk, half plane and quadrant (`exact_bulk.py`, `exact_half.py`, `exact_quad.py`); quadrant representation and Theorem A (`thmB*.py`, `halfplane.py`, `bulk.py`, `amp.py`) | Sec. IV |
| `symbol/` | Scalar-symbol reductions R_bulk = ½\|det c_{j−k−1}\|, R_half = ½\|det(c_{j−k−1} − c_{j+k+3})\|; Toeplitz+Hankel minors and D_r asymptotics (`th1–th3.py`, `rh1–rh2.py`) | Sec. IV B–C |
| `aztec/` | Exact Aztec-diamond boundary/interior weights (`az_*.py`, `run_az.py`); Kuo-condensation and zipper checks for the elementary proofs (`kuo*.py`); the Chhita–Young pair formula (`r12a.py`, `r12b.py`, `r12m.py`) | Sec. V, App. A |
| `limitshape/` | Cohn–Kenyon–Propp limit shapes by FEM (`varsolve.py`, `var2.py`); boundary shifts (`var_arc.py`); tangent continuation of the dual variables (`tangent_lib.py`, `tan_az.py`, `tan_house.py`); house-domain exact rates (`house_exact_k.py`) and diagnostics (`mu_*.py`, `az_diag.py`, `house_refine.py`, `orient.py`) | Sec. VI A–D |
| `lozenge/` | Hexagon Kasteleyn counts (`hexlib.py`); dent-pair tables (`pair.py`); corner law (`row1.py`); LGV formula and its derivation/verification (`lgv.py`, `lgv2.py`); rates and the layer rule (`rate.py`, `lozfig.py`) | Sec. VI E, Fig. 9, Table I |
| `data/` | Small result files used for the figures (JSON, one FEM limit shape as pickle) | — |
| `paper/` | `main.tex` (REVTeX 4.2), `mkfigs.py`, `figs/` | — |

Scripts import helpers from sibling folders through a one-line path header, so they can be run from anywhere. A few scripts `exec` a sibling file by relative name (`symbol/th2.py`, `symbol/rh1.py`, `lozenge/pair.py`, `lozenge/rate.py`, `lozenge/lozfig.py`); **run those from inside their own folder**.

## Reproducing key results

```bash
# 8x8 board: 12988816 tilings; mutilated board (a1, h8 removed): 0 tilings
cd dp && gcc -O2 -o cnt cnt.c && ./cnt 8 8 && ./cnt 8 8 0 0 7 7

# Chhita–Young boundary pair formula (Eq. 13), exact for n <= 20
cd aztec && python3 r12a.py && python3 r12m.py

# Elementary proof ingredients (zippers, Kuo condensation) for Eqs. (11)-(12)
cd aztec && python3 kuo3.py && python3 kuo5.py

# Certified half-plane / bulk ratios (Arb), e.g. r = 15..1023
cd exact && python3 exact_half.py 15 31 63 127 && python3 exact_bulk.py 15 31 63

# Scalar-symbol reduction check (half plane and bulk)
cd symbol && python3 th2.py && python3 th3.py

# Aztec limit shape + tangent continuation (Fig. 7)
cd limitshape && python3 fem_save.py az 48 && python3 tan_az.py

# House domain: exact lattice rates and parameter-free prediction (Fig. 8)
cd limitshape && python3 house_exact_k.py 0.5 32 48 64 96 128
python3 fem_save.py house 64 0.5 && python3 tan_house.py fem_house64_0.5.pkl

# Lozenge hexagon: LGV vs Kasteleyn, then the layer-rule test (Fig. 9, Table I)
cd lozenge && python3 lgv2.py && python3 rate.py && python3 lozfig.py

# Figures and paper (a compiled main.pdf is included)
cd paper && python3 mkfigs.py && latexmk -pdf main.tex
```

Typical run times are seconds to minutes on a laptop. The largest runs, such as N = 96 FEM limit shapes, n ≥ 200 LGV determinants and r = 1024 Arb determinants, take up to tens of minutes.

## Conventions

- Cell `(r, c)` is black if `r + c` is even. White defects have charge `e = −1`, black defects `e = +1`.
- Kasteleyn weights are 1 on horizontal edges and ±1 (alternating by column) or `i` on vertical edges. Interior defects are handled by sign-flip strings.
- The Aztec diamond is `aztec(n)` in `lattice/kreg.py`: cells with |c − n + ½| + |r − n + ½| ≤ n.
- The hexagon is `hexagon(a, b, c)` in `lozenge/hexlib.py`, with sides a (bottom), b (right) and c. Dents are counted from the bottom-right corner.

## Status of results

See Table II of the paper. It distinguishes rigorous proofs, exact arithmetic, numerical evidence and conjectures (the layer rule, Eq. 18).

## Changes

See `CHANGELOG.md`.

## License

MIT (see `LICENSE`).
