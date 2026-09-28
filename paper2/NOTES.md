# Paper II Task 1 notes (branch paper2)
- pftri.py: Pfaffian orientation, square + NE diagonal weight t. Validated vs dp/tri_t.c and brute force.
- pairfast.py: Schur-complement pair correlation. Accurate for t <= 0.3 on 96x96;
  for t >= 0.5 at large r it loses digits (g=G(u,v) ~ e^{-kappa r} tiny) -> data/paper2_pair96.json
  entries t>=0.5 are NOT reliable. Use data/paper2_pair_bulk.json (64x64 direct slogdet) there.
- deriv.py uses one-sided differences (bias ~0.2% at N=45). ln R is exactly even in t
  (Z(-t)=±Z(t)); use central differences in future.
- Corner ratio: ln R(N,t)-ln R(N,0) = a2(N) t^2 + ..., a2 = c N^2 + e N lnN + f, c~0.015 same for both
  corners; fit on N<=81 predicts N=101 to 0.1% (A) / 0.01% (B). Consistent with nu=1.
- t=1/2: C = [1-(sqrt3-2)^(|dr|+|dc|)]/6 for dr*dc>=0 (checked to ~1e-10), not for dr*dc<0.
- Monomer decay (commit after b78fa42): pairdirect.py = sparse LU of the correct string matrix (no
  cancellation; agrees with dense slogdet to 1e-11 for all t). Data: data/paper2_pair_big.json
  (L=256; t=0.1 on L=384; r<=80; values reproducible across L=128/160/256 to 1e-12).
  For t<=0.3: C - Cinf = (-1)^(r+1) A r^-1 e^{-2 kappa r} (1 + O(1/r)), kappa = artanh t.
  Residual rate delta0 = 0 within 5e-4..2e-3 (mono_delta0.py); candidate 2artanh(t/sqrt(1-t^2)) excluded.
  A(t) ~ 0.2 for small t, A -> 0 at t=1/2 (pure e^{-ln(2+sqrt3) r} there), sign change suggested at t=0.6.
  t>~0.6: gap minimum leaves phi=0 -> incommensurate oscillation (bulk, L-independent).
- pairsparse.py / pairfast.py keep the Schur-complement route; unreliable when g=G(u,v) is tiny.
