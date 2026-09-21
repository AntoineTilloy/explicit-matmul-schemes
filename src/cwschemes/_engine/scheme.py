"""Exact coefficient evaluation / certificate replay. No discovery entry points."""
from .._resources import root as resource_root
from bisect import bisect_right
from fractions import Fraction as F
from .restrictions import completion, SEED, audit, rank_bound
from .patches import packet_plan, group_cost, patch_cost

class CompactPacket:

    def __init__(self, q, rows=SEED, m=2, patched=True):
        if q < 1:
            raise ValueError('q must be positive')
        self.q = q
        self.rows = rows
        self.m = m
        self.patched = patched
        if patched:
            self.plan = packet_plan(rows, m, q)
        else:
            witness = audit(rows, m)
            self.plan = dict(witness=witness, q=q, n=q ** m, original_rank=rank_bound(q, witness))
            self.plan['groups'] = [dict(pure=s, zero=False, patches=[], remaining=[i for i in range(3 * m) if not s >> i & 1], families=[sorted({r[j] for r in rows if r[j] & s == s}) for j in range(3)]) for s in self.plan['witness']['surviving_pure_subsets']]
            self.plan['rank'] = self.plan['original_rank']
        self.forms, self.weights = completion(q)
        self.rank = self.plan['rank']
        self._ends = []
        offset = 0
        for group in self.plan['groups']:
            offset += group_cost(q, group)
            self._ends.append(offset)
        assert offset == self.rank
        self._families = [[frozenset(s) for s in g['families']] for g in self.plan['groups']]

    def term(self, index):
        if not 0 <= index < self.rank:
            raise IndexError(index)
        q = self.q
        gi = bisect_right(self._ends, index)
        g = self.plan['groups'][gi]
        if gi:
            index -= self._ends[gi - 1]
        blocks = [dict(kind='local', positions=[i]) for i in g['remaining']] + g['patches']
        radices = [q + 1 if b['kind'] == 'local' else patch_cost(q, b) for b in blocks]
        choices = [0] * len(blocks)
        for i in reversed(range(len(blocks))):
            index, choices[i] = divmod(index, radices[i])
        return (gi, list(zip(blocks, choices)))

    def matrix_variable(self, product_index, mode, row, col):
        """Map A[row,col], B[row,col], or output C[row,col] to local indices."""
        if not 0 <= product_index < len(self.rows) or mode not in (0, 1, 2):
            raise IndexError((product_index, mode))
        n = self.q ** self.m
        if not (0 <= row < n and 0 <= col < n):
            raise IndexError((row, col))
        masks = self.rows[product_index]
        rowmask, colmask = [(masks[1], masks[2]), (masks[2], masks[0]), (masks[1], masks[0])][mode]
        values = [0] * (3 * self.m)
        for mask, index in [(rowmask, row), (colmask, col)]:
            positions = [i for i in range(3 * self.m) if mask >> i & 1]
            for i in reversed(positions):
                index, digit = divmod(index, self.q)
                values[i] = digit + 1
        return tuple(values)

    def factor_entry(self, mode, index, variable):
        if mode not in (0, 1, 2):
            raise ValueError('mode must be 0,1,2')
        if len(variable) != 3 * self.m or any((not 0 <= x <= self.q for x in variable)):
            raise ValueError('invalid variable')
        gi, blocks = self.term(index)
        g = self.plan['groups'][gi]
        q = self.q
        mask = sum((1 << i for i, x in enumerate(variable) if x == 0))
        if mask not in self._families[gi][mode]:
            return F(0)
        value = self.weights[-1] ** g['pure'].bit_count() if mode == 2 else F(1)
        for patch, choice in blocks:
            pos = patch['positions']
            kind = patch['kind']
            if kind == 'local':
                a = self.forms[choice][variable[pos[0]]]
                if mode == 2:
                    a *= self.weights[choice]
            elif kind == 'identity':
                x = variable[pos[0]]
                a = F(x == 0 if mode == patch['zero_mode'] else x == choice + 1)
            elif kind == 'two_matvec':
                a = self._two_matvec_entry(mode, patch, choice, variable)
            else:
                t, z = divmod(choice, q + 1)
                pair = patch['pair']
                other = 3 - sum(pair)
                i, j = pos
                a = self.forms[t][variable[i]] * self.forms[z][variable[j]]
                if mode == 2:
                    a *= self.weights[t] * self.weights[z]
                if mode == other:
                    alpha = F((3 if t < q else -4 * q) * (3 if z < q else -4 * q), 16 * q * q)
                    if 2 in pair:
                        alpha *= self.weights[q] ** 2 / (self.weights[t] * self.weights[z])
                    last = self.forms[q][variable[i]] * self.forms[q][variable[j]]
                    if mode == 2:
                        last *= self.weights[q] ** 2
                    a -= alpha * last
            value *= a
            if not value:
                break
        return value

    def _two_matvec_entry(self, mode, patch, term, x):
        q = self.q
        i, j = patch['positions']
        h = patch['matrix_mode']
        a, b = (x[i] - 1, x[j] - 1)
        if term < q * q:
            r, c = divmod(term, q)
            if mode == h:
                return F(a == r and b == c)
            return F(a == r and b == -1 or (a == -1 and b == c))
        if term < q * q + q:
            r = term - q * q
            if mode == h:
                return F(-int(a == r and b >= 0))
            return F(a == r and b == -1)
        c = term - q * q - q
        if mode == h:
            return F(-int(b == c and a >= 0))
        return F(b == c and a == -1)

    def certificate(self):
        return dict(format='cw-square-packet-v1', coefficient_field='Q', patched=self.patched, definition='scheme.CompactPacket.factor_entry(mode, term, variable)', local_completion='sum Li^3 - (4q/3)Lstar^3 + (q/3)x0*y0*z0', characteristic_exclusions='denominators require characteristic not dividing 6q', **self.plan)
