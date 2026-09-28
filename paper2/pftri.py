"""Pfaffian (Kasteleyn) orientation for the square lattice + NE diagonals (weight t),
interpolating square (t=0) and triangular (t=1) dimers on an Lx x Ly board.
Site (r,c), index r*Lx+c, y=r upward.  Orientation (derived from 'odd clockwise per face'):
  horizontal (r,c)->(r,c+1) always east;  vertical in column c: up if c even, down if c odd;
  diagonal (r,c)-(r+1,c+1): NE iff vertical of column c+1 is up.
Monomer at interior v: flip sign of edges crossed by the half-line x=c+eps, y<r
(horizontal (r',c)-(r',c+1) and diagonal (r',c)-(r'+1,c+1), r'<r), which restores the
odd-clockwise condition on the face created by deleting v.
Z(board minus v)^2 = det K_v  (K_v = string-modified K with row/col v deleted)."""
import numpy as np

def kmat(Lx, Ly, t, v=None):
    n = Lx*Ly; K = np.zeros((n, n))
    vr, vc = (v if v is not None else (-10, -10))
    def add(a, b, w):   # oriented a -> b
        K[a, b] += w; K[b, a] -= w
    up = lambda c: 1 if c % 2 == 0 else -1
    for r in range(Ly):
        for c in range(Lx):
            i = r*Lx + c
            flip = -1 if (c == vc and r < vr) else 1
            if c+1 < Lx:
                add(i, i+1, flip)
            if r+1 < Ly:
                j = i + Lx
                if up(c) == 1: add(i, j, 1)
                else: add(j, i, 1)
            if t != 0 and c+1 < Lx and r+1 < Ly:
                j = i + Lx + 1
                if up(c+1) == 1: add(i, j, t*flip)
                else: add(j, i, t*flip)
    return K

def logZ(Lx, Ly, t, v=None):
    """ln Z for the board (v=None, needs even size) or board minus site v=(r,c)."""
    K = kmat(Lx, Ly, t, v)
    if v is not None:
        i = v[0]*Lx + v[1]
        keep = np.r_[0:i, i+1:Lx*Ly]
        K = K[np.ix_(keep, keep)]
    s, ld = np.linalg.slogdet(K)
    return 0.5*ld if s > 0 else -np.inf

def kmat_multi(Lx, Ly, t, mons):
    """K with a downward string from every monomer in `mons` (list of (r,c))."""
    n = Lx*Ly; K = np.zeros((n, n))
    up = lambda c: 1 if c % 2 == 0 else -1
    def flip(r, c):
        f = 1
        for (vr, vc) in mons:
            if c == vc and r < vr: f = -f
        return f
    for r in range(Ly):
        for c in range(Lx):
            i = r*Lx + c; f = flip(r, c)
            if c+1 < Lx: K[i, i+1] += f; K[i+1, i] -= f
            if r+1 < Ly:
                j = i+Lx
                if up(c) == 1: K[i, j] += 1; K[j, i] -= 1
                else: K[j, i] += 1; K[i, j] -= 1
            if t != 0 and c+1 < Lx and r+1 < Ly:
                j = i+Lx+1
                if up(c+1) == 1: K[i, j] += t*f; K[j, i] -= t*f
                else: K[j, i] += t*f; K[i, j] -= t*f
    return K

def logZ_multi(Lx, Ly, t, mons):
    K = kmat_multi(Lx, Ly, t, mons)
    rm = sorted(r*Lx+c for r, c in mons)
    keep = np.setdiff1d(np.arange(Lx*Ly), rm)
    s, ld = np.linalg.slogdet(K[np.ix_(keep, keep)])
    return 0.5*ld if s > 0 else -np.inf
