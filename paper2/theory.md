# Paper II, Task 0 — Coulomb-gas predictions as a function of stiffness g

Conventions of Paper I: normalized height h = h_T/2; a defect (monomer) is a vortex of h
(Thurston winding ±4); in the dual field a bulk monomer creates (2/π) ln|z−z0|; a boundary
height jump J at c creates (J/π) ln|z−c| in the conjugate field. Free-fermion point: g = 1/2.

## 1. Stiffness enters only as an overall factor

For a fixed configuration of vortex-type charges (the defect e and the corner charges q_c),
the Gaussian weight is exp(−g·E), with E depending only on the geometry (positions,
conformal map, and the cutoff self-energy that produces the |f'| factor). Corner charges are
fixed by topology, q_c = −J_c/4 (boundary height jump / monomer winding), independent of g.
Hence, within the Gaussian theory,

    ρ_g(z) ∝ [ρ_{1/2}(z)]^{2g}.

With Paper I's Eq. (3), ρ_{1/2} ∝ |f'|^{1/4} (Im f)^{1/4} Π_c |f − f_c|^{e q_c}:

    ρ_g(z) ∝ |f'|^{g/2} (Im f)^{g/2} Π_c |f − f_c|^{2g e q_c}.            (T1)

The monomer scaling dimension is x_m = g/2 (x_m = 1/4 at g = 1/2), and a pair of defects
carries |f1 − f2|^{g e1 e2} |f1 − conj f2|^{g e1 e2}.

## 2. Corner exponent

Near a corner of interior angle θ, f ~ z^{π/θ}:

    γ(g) = 2g [ (π/θ)(1/2 + e q_c) − 1/4 ]  =  (π/θ)(g + 2g e q_c) − g/2.   (T2)

g = 1/2 gives Paper I Eq. (4), γ = (π/θ)(1/2 + e q_c) − 1/4.

## 3. Checks at g = 1/2 against Paper I (all reproduced)

| quantity | general g | g = 1/2 | Paper I |
|---|---|---|---|
| white defect at removed corner (θ=π/2, q_c=3/4, e=−1) | γ = −3g/2 | −3/4 | −3/4 (slopes 0.754, 0.750) |
| odd×odd board, corner/centre ratio | N^{−g/2} | N^{−1/4} | N^{−1/4} |
| bulk pair correlation | d^{−g} | d^{−1/2} | d^{−1/2} |
| pair along a free edge | d^{−2g} | d^{−1} | d^{−1} |

Check of the removed-corner case: (π/θ)=2, e q_c = −3/4 → 2g[2(1/2 − 3/4) − 1/4] = 2g(−3/4) = −3g/2.

## 4. What is and is not universal (predictions for Task 4)

- Exponents: all scale by the single factor 2g — the ratio of any two exponents is g-independent
  (a sharp, parameter-free test for the interacting-dimer Monte Carlo).
- q_c, the sign of the image (Neumann for the dual field), and the conformal map f: g-independent.
- Amplitudes: not predicted.

## 5. Assumptions and caveats (to be checked, not assumed)

1. The dual field keeps Neumann boundary conditions at a free edge for all g (same-sign image).
2. No electric (height-locking) operators are relevant: fails near the columnar transition,
   where cos(4πh)-type terms become relevant.
3. Mapping to Alet et al. (PRL 94, 235702, 2005; arXiv cond-mat/0501241) — CONFIRMED from the
   original: alpha(e,m) = e^2/g + g m^2, alpha_m = g (monomer M(x) ~ x^{-alpha_m}), alpha_d = 1/g,
   g = 1/2 at T = infinity. Aligning v = -1: g grows to g(Tc) = 4 at Tc = 0.65(1) (|v|/Tc ~ 1.54),
   KT via the p = 4 locking term. Repulsive v = +1: g decreases to 0 at T_t = 0.449(1).
   So our g is exactly theirs; Task 4 can measure g from winding fluctuations, their Eq. (5).
4. Frozen regions / μ (Aztec) are outside this Gaussian argument.

## 6. Honeycomb (Task 0, item 2)

- Winding: honeycomb Thurston units (empty edge +1, dimer −2 around a vertex of degree 3):
  a monomer winds by ±3, so the corner rule generalizes to q_c = −J_c / w with w = 3 when J_c is
  measured in the same units. The monomer is still a unit magnetic charge, so x_m = g/2 = 1/4 at
  the free-fermion point (consistent with Ciucu's q1 q2/2 hole interactions; to be checked
  numerically with a bulk pair on a large flake).
- Which boundaries are flat? In brick-wall form the honeycomb is the square lattice with one of the
  four dimer types missing, so the slope polygon shrinks from a square to a triangle. The square
  lattice's zero tilt (forced by an axis-parallel rectangle) lies ON an edge of that triangle:
  brick-wall rectangles are frozen/degenerate. A flat liquid needs boundaries whose tilt is interior
  (all three orientations 1/3).
- Numerical test (paper2/honeyflake.py, honey_diag.py; exact inverse of the biadjacency matrix):
  armchair D6 flakes are flat liquids: bulk mean|P−1/3| = 0.035, 0.017, 0.012 (core) and 0.052,
  0.028, 0.021 (mid) for R = 12, 20, 28 (366–2106 sites), no near-frozen edge, deviation confined
  to a lattice-scale boundary layer. Zigzag D6 flakes (= lozenge hexagons) are not: mid-annulus
  deviation stays ~0.06 and near-frozen edges grow 18 → 108 (arctic behaviour).
- Pitfall found: averaging P by orientation over annuli hides the anisotropy (D6 symmetry makes each
  annulus average 1/3). Use local |P − 1/3|.
- Next: boundary height jumps J_c at the 120° corners of armchair flakes → q_c; then Task 2 proper
  (one forced defect: remove a site, density vs Eq. T1 with f for the hexagon).

## 7. Corner charges of armchair flakes (paper2/honey_corner.py)

Boundary height from the dual (triangular-lattice) boundary cycle, step +1 (−1) when the interior
vertex next to the boundary edge is A (B); configuration-independent. Averaging over each straight
armchair edge (12° away from the corners) and differencing gives integer jumps. Scan R = 8..30:

    J_c = 0, ±1 or ±2  (alternating in sign around the flake; total 0)   →   q_c = −J_c/3 ∈ {0, ∓1/3, ∓2/3}

The value is set by the microscopic corner termination, not by the edge length (R = 9 and 10 both
have 42 boundary steps but J = 0 and ±2). Three corner types, versus two (±1/4) on the square lattice.
Adjacent corners have opposite sign because a 60° rotation about the face centre swaps A and B.
Prediction for Task 2: a defect of charge e near a 120° corner (π/θ = 3/2) has
γ = 2g[(3/2)(1/2 + e q_c) − 1/4] = (3/2)(1/2 + e q_c) − 1/4 at g = 1/2, i.e. γ ∈ {1/2, 0, 1} for
e q_c = 0, −1/3, +1/3 and γ ∈ {−1/2, 3/2} for e q_c = −2/3, +2/3.

## 8. Task 2 first results — the naive J rule FAILS on armchair corners; a half-unit correction works

Setup (paper2/honey_defect.py): A-centred armchair flakes have n_A = n_B + 1, so one A defect
(e = +1) is forced without removing any site. Z(flake − v) by sparse LU with a downward string;
checked against brute-force enumeration for every defect site (R = 3.5, 4.5). Two families by
corner termination (J measured as in §7, mod 3 because of the defect's winding):
  type I  (R = 10, 15, 18, 20, 23, 28, 31, 33, 36): J = 0 (30°,150°,270°), +1 (90°,210°,330°)
  type II (R = 9, 12, 14, 17, 22, 25, 27, 30, 35, 38): J = −1 and +2.

Observations:
- type I: all boundary A sites have exactly the same defect probability, and the profiles from a
  J = 0 corner and a J = +1 corner into the centre coincide to 3 digits, although the corners differ
  geometrically (B count 3 vs 5 within distance 3). The A set is C6-symmetric; the defect density
  inherits that symmetry → type I is Temperleyan (density governed by the A lattice), ρ ∝ R^{1/4}.
- Corner exponent from corner/centre ~ R^{−γ} (paper2/honey_scaling.py, R ≤ 38):
    type I   30°, 90°:  γ = 0.29 (fit), local 0.24–0.36
    type II  30°:       γ = 0.76 (fit), local 0.65–0.79
    type II  90°:       γ = −0.25 (fit), local −0.22 to −0.26
  Naive q_c = −J/3 predicts 1/2, 0 (type I) and 1, −1/2 (type II): wrong.
- Effective rule that fits all four: J_eff = J + 1/2 at the 30°-type corners and J − 1/2 at the
  90°-type corners, q_eff = −J_eff/3:
    type I  q_eff = −1/6 at all six corners → γ = 1/4 (Temperley)
    type II q_eff = +1/6, −1/2              → γ = 3/4, −1/4
  Totals: Σ q_eff = −1 in both families (neutralizes the defect).
- Interpretation (hypothesis): on armchair edges the boundary height oscillates on the lattice
  scale, and averaging the dual-edge midpoint heights sets the reference differently by 1/2 unit
  depending on the corner termination. Needs a principled coarse-graining (e.g. average over the
  faces adjacent to the edge) before it can be called a rule. Local exponents still oscillate;
  larger R and the full density formula (Schwarz–Christoffel map of the hexagon) are next.
- Same-sign image test at edge midpoints: not done yet (in type I edge sites equal corner sites).

## 9. Principled corner jump from the interior mean height (paper2/honey_meanh.py, honey_fitJ.py)

Balanced face-centred armchair flakes (no defect, so <h> is single-valued). Crossing a bond with A on
the left, dh = 1 − 3x, so <dh> = 1 − 3P with exact P from the inverse; integrated on the dual graph
(curl error ≤ 4e−15). The boundary faces themselves give levels that differ from the dual-edge-midpoint
reading by half units — the source of the ambiguity. The continuum boundary value of an edge is the
value the interior mean field approaches, so fit interior faces (hexagonal radius < 0.7) to
<h> = c + a u(z), u = harmonic function on the regular hexagon with boundary values ±1 on alternating
edges (finite differences, 401²). J_eff at the 30° corner = −2a.

  family (naive J at 30°)   predicted J_eff   measured (R = ...)
  J = 0   (R 12..30)        +1/2              +0.32, (+0.69 at R=17), +0.37, +0.42, +0.42, +0.46
  J = −1  (R 11..29)        −1/2              −0.24, −0.35, −0.61, −0.61, −0.42, −0.46
  J = −2  (R 15..28)        −3/2              −1.32, −1.33, −1.09, −1.39, −1.24

All three are clearly away from the naive integers and on the predicted half-integer side; J = 0
approaches +1/2 with R. Fit-window systematics ±0.1–0.15; J = −2 not yet at −1.5.
Conclusion: the half-unit rule is supported by an independent, fit-free-in-J definition, but the
convergence is not yet precise. Next: exact harmonic measure (Schwarz–Christoffel) instead of the
grid u, fit that excludes a lattice-scale boundary layer of fixed width, and R up to ~60.

## 10. Same-sign image on armchair edges (paper2/honey_image.py)

Edge-midpoint A site vs centre, A-centred flakes (one forced A defect):
  type I   edge/centre 0.914 → 0.628 (R 10 → 36): exponent −0.29 (local −0.24…−0.36)
  type II  edge/centre 0.825 → 0.633 (R 9 → 38):  exponent −0.24 (local values oscillate with R)
Inward profile from the edge midpoint (R = 36, 38): d ln ρ / d ln d = +0.20, +0.21 over 2 < d < R/4
(ρ increases away from the edge).
Same-sign image (dual Neumann) predicts R^{−1/4} and d^{+1/4}; opposite sign predicts R^{+1/4},
d^{−1/4}. Verdict: SAME SIGN on honeycomb armchair edges; the slightly small inward slope is the
expected saturation of the (Im f) factor at d ~ R/4.

## 11. Novelty check against Ciucu's programme (web, 2026-09-29)

Known (cite, do not claim):
- Bulk Coulomb law for gaps/holes on the hexagonal lattice, charge = (#down − #up triangles)
  (Ciucu, Mem. AMS 2005; PNAS 2008; Trans. AMS 2010).
- Image rule by boundary type: constrained zig-zag boundary → image with the SAME orientation
  (Neumann / insulating); free boundary (lozenges may protrude) → OPPOSITE orientation, gap attracted
  like a charge to a straight conductor, correlation ~ 1/(4πr) (Ciucu–Krattenthaler, CMP 302 (2011);
  Ciucu, 90° wedge with mixed boundary conditions, arXiv 1403.0925).
- 60° angle with zig-zag sides: correlation = product of distances to five images
  (Ciucu–Fischer, arXiv 2014) — i.e. images only, q_c = 0 in our language.
Our §10 (same sign on armchair edges) is therefore an extension of the constrained-boundary rule,
not a new rule.

Likely new (to be stated carefully in Paper II):
- Corner charges from the boundary-height jump, and on the honeycomb the half-unit correction
  (J_eff = J ± 1/2, §8–9); corner types {−1/6} (Temperley) and {+1/6, −1/2}.
- Forced single defects in odd (unbalanced) regions and their full density, corner exponents γ.
- Temperley / non-Temperley classification of armchair flakes by corner termination.
- The single g-covariant formula γ(g) tying square, honeycomb and (Task 4) interacting dimers.
Check to do: does our framework reproduce Ciucu–Fischer's 60° zig-zag result with q_c = 0?
(a consistency test of the corner-charge rule on a boundary type Ciucu solved exactly).

## 12. 60° armchair corners (Ciucu–Fischer boundary type) — rule-independent test (honey_triangle.py)

Note: Ciucu–Fischer's "zig-zag lines" are zig-zag paths on the TRIANGULAR lattice = armchair edges of the
honeycomb, i.e. our boundary type (flat liquid, same-sign images). Their 60° wedge result (gap + five
images) is the q_c = 0 case.
C3-symmetric armchair triangles centred on an A site: n_A − n_B = 1 for every R = 6..59 (one forced A
defect). By symmetry all three corners are equal, and neutrality forces q_c = −1/3 each — independent of
any rule for J. Naive boundary jump: J = +1 at all three corners (the −2 at 180° includes the cut: −2+3),
so q_c = −J/3 = −1/3: at 60° corners the naive rule is exact (no half-unit shift).
Prediction (g = 1/2, π/θ = 3): γ = 3(1/2 − 1/3) − 1/4 = 1/4 (vs 5/4 without corner charge).
Measured corner/centre, R = 8..58 (up to 13 262 sites): γ = 0.219 (fit); same-phase subsequence
R = 8, 16, 24, 32, 40: local 0.19, 0.22, 0.23, 0.25; R = 12, 20, 28, 36: 0.11, 0.27, 0.27.
Verdict: corner charge + same-sign images + γ formula confirmed at 60° corners; 5/4 excluded.
Consequence: the half-unit shift found at 120° armchair corners (§8–9) is specific to those corners
(where the two edges' lattice oscillations meet with a phase offset), not a failure of the framework.

## 13. Edge-class boundary-layer shift: explanation and a blind prediction test (honey_para.py, honey_para2.py)

Boundary heights are NOT ambiguous: along an armchair edge the steps alternate +1, −1, so the midpoint
height is exactly constant (half-integer); a corner shows a jump only where two equal steps meet.
(The earlier "phase offset" guess is rejected.)
Explanation: armchair edges come in two classes by outward normal (0/120/240° = class +, 60/180/300° =
class −). A 60° rotation maps one class to the other and swaps A↔B, which flips h → −h. Hence the
lattice boundary layer shifts the continuum boundary value by +δ on class + and −δ on class − edges.
  - 120° corners join different classes → J_eff = J ± 2δ;  60° corners join the same class → J_eff = J.
  - Calibration (§8–9): 2δ = 1/2; going ccw from a class + edge to a class − edge adds +1/2, the reverse −1/2.
  This explains §8–9 (120°: ±1/2) and §12 (60°: 0) at once.
Blind test: armchair parallelograms (A-centred, one forced A defect), two 60° and two 120° corners; exact
naive J with the branch cut fixed at the 120°-direction corner (ΣJ = 3). Rule fixed beforehand.
  F1 (R = 9..38), J(30,120,210,300) = (−1, 1, 2, 1):
     predicted γ = +0.75, +0.25, −0.25, +0.25 (naive: +1.00, +0.25, −0.50, +0.25)
     measured  γ = +0.71, +0.24, −0.27, +0.24 (last 5 R)
  F2 (R = 10..36), J = (0, 1, 1, 1):
     predicted +0.25 at all four (naive: +0.50, +0.25, 0, +0.25);  measured +0.28 at all four
Verdict: the edge-class shift rule predicts all eight corners within 0.02–0.05; the naive rule fails
at every 120° corner by 0.25–0.3. Honeycomb corner charge: q_c = −(J + σ/2)/3, σ = ±1 at 120° corners
by the class order, σ = 0 at 60° corners.
Open: derive δ = 1/4 from the lattice (e.g. mean height next to a single armchair edge in a half plane).
