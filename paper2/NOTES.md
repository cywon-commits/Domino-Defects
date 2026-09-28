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
