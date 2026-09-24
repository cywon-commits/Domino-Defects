import sys as _s, os as _o; _s.path[:0]=[_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..',_d) for _d in ('lattice','exact','aztec','limitshape','lozenge')]  # repo path setup
from exactG import *
from flint import arb, ctx
ctx.prec=200
cols=plane_table(6,10); pi=arb.pi()
import numpy as np
from bulk import Plane
