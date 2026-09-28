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
3. Mapping to Alet et al. (PRL 94, 235702, 2005): the bulk monomer correlation decays as
   r^{−g}. I recall that their normalization has dimer ~ r^{−1/g}, monomer ~ r^{−g} with g = 1/2
   at v = 0, which matches ours; this must be confirmed from the original before Task 4.
4. Frozen regions / μ (Aztec) are outside this Gaussian argument.
