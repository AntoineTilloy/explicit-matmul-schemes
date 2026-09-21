"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from flint import fmpq as Q, fmpq_mat as QM
import numpy as np
ZERO = Q(0)

def qm(a):
    return QM(a.tolist())

def from_qm(a):
    return np.array([[a[i, j] for j in range(a.ncols())] for i in range(a.nrows())], dtype=object)
