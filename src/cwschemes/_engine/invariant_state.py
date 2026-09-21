"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
import numpy as np

class State:

    def reindex(self):
        self.patterns = [[H for H in self.kept if self.ordinary[H] & I == I] for I in range(1 << self.k)]
        self.lookups = [{H: j for j, H in enumerate(ps)} for ps in self.patterns]

    def zeros(self):
        nonzero = [set() for _ in range(3)]
        for mode in range(3):
            for I, ps in enumerate(self.patterns):
                nonzero[mode].update((H for j, H in enumerate(ps) if np.any(self.coefs[mode][I][:, j])))
        return [H for H in self.kept if any((H not in present for present in nonzero))]

    def orbit_size(self, I):
        if not hasattr(self, 'qs'):
            return self.q ** I.bit_count()
        value = 1
        for i in range(self.k):
            if I >> i & 1:
                value *= self.qs[i]
        return value

    def rank(self):
        return sum((self.orbit_size(self.ordinary[H]) for H in self.kept))
