# Small rational rectangular compositions

The initial size-106 and size-146 schemes used **424,690** and **1,044,715** products.
They improve the strengthened rational composition reference counts 431,011
and 1,063,741 by 1.47% and 1.79%. Size 432 now uses **19,665,128** products,
through the routine product of the current published LITA size-18 and size-24
schemes, with ranks 3,031 and 6,488. This refresh is distinct from the new
uneven rectangular compositions.

The subsequent [extended rectangular search](extended_rectangular.md) lowers
size 146 to **1,040,656** products and updates several other entries below.
The counts in this page document the earlier search stages.

For a split `n = p + q`, `p >= q`, Strassen's identity can use different
rectangular inner algorithms for its seven calls. Its cost is bounded by

```text
R(p,p,p) + 3 R(p,p,q) + 3 R(p,q,q).
```

The concrete recipes retain each call's orientation and exact dimensions.
At size 106 the split is 54+52; at size 146 it is 74+72. Rectangular products
use padding, tensor products and axis splits, with cheap published narrow
rectangular schemes for the remainders. The square output has exact rational
coefficient access through the usual `cwschemes.load(n)` API.

The format `exact-rectangular-recipe-v1` is a directed acyclic graph. Nodes
are frozen sparse rational coefficient objects, classical products, tensor
products, splits along one matrix-index axis, padding, index permutations,
or uneven Strassen compositions. Every node checks its dimensions and rank.
The immutable primitive files are individually hash checked and materialized
on demand from the package's content-addressed store.

The primitive LITA revision is
[472a79f](https://github.com/khoruzhii/lita/tree/472a79fbc8248423a97bebb3aacbc89e769b9b3e).
All saved LITA coefficients were regenerated from the pinned generator,
with exact equality of every stored array. Their identities rely on the
published LITA construction theorem. This regeneration is not described as
an exhaustive check of the large LITA tensors. The source snapshot is in
`sources/lita_current.py`, with its original attribution and license scope
described in [third_party.md](third_party.md).

The imported rectangular primitive tensors were checked completely over Q,
using exact integer arithmetic on these integer coefficient schemes.
Small uneven compositions at sizes 3, 4, 5 and 6 also passed complete tensor
checks, independently of the dimension/count argument. The package test
checks all 15,625 size-five coefficients using unequal 3+2 blocks.

Evidence is in `checks/small1024_rectangular*.json`. Some historical research
records contain the original staging paths; the current immutable hashes
and catalogue entries identify the published objects. Huge target tensors
are established compositionally and have not been expanded in full.

## Independent rectangular splits and recursive feedback

For a rectangle with dimensions `(p0+q0, p1+q1, p2+q2)`, with each
`pi >= qi > 0`, the seven Strassen calls have dimensions

```text
(p0,p1,p2), (q0,p1,p2), (p0,p1,q2), (q0,q1,p2),
(p0,q1,q2), (q0,p1,q2), (p0,q1,p2).
```

These are all corners except `(q0,q1,q2)`. Zero padding and cropping use
the independent split on each matrix-index axis. Every multiplication has
an A-linear form on the left and a B-linear form on the right; the identity
therefore permits recursive substitution by matrix blocks. The recipe's
`p` and `q` fields accept three-element vectors. Existing scalar fields keep
their original meaning of equal splits on all three axes.

The October 13 search iterates this rule within a rectangular table through
dimension 128, together with padding, axis splits, exact-divisor tensor
products and feedback from square improvements. Each changed cost receives
a new immutable recipe version, so later padding cannot create a cycle.
Square results are compared with the minimum of the previous full rectangular
closure, the public dense closure and the catalogue at commit `e313024`.

| Square size | Previous envelope | New products |
|---:|---:|---:|
| 125 | 681,408 | 680,814 |
| 154 | 1,214,333 | 1,211,453 |
| 162 | 1,401,348 | 1,394,733 |
| 171 | 1,643,187 | 1,641,699 |
| 178 | 1,819,017 | 1,814,857 |
| 186 | 2,055,249 | 2,052,507 |
| 187 | 2,111,184 | 2,109,033 |
| 194 | 2,322,311 | 2,309,010 |
| 195 | 2,372,355 | 2,370,114 |
| 218 | 3,206,454 | 3,205,147 |
| 234 | 3,912,786 | 3,911,301 |
| 242 | 4,261,876 | 4,255,382 |

A matched control with feedback and tensor products but without the new
rectangular Strassen rule improves sizes 194 and 242. The new rule improves
the control at eleven sizes and ties it at 242; the latter gain comes from
the shared closure changes. There is no improvement at size 1024.

All saved search nodes pass an independent integer recipe and dependency
audit. The coefficient engine is tested by reconstructing every coefficient
of MM(3,4,5), with splits `(2+1,3+1,3+2)`, under all six axis permutations.
Primitive coefficient files and proof scope are recorded in
`checks/asymmetric_strassen.json`. As above, the large final tensors are
certified by composition and the stated primitive identities.

Larger current-LITA leaves use indexed access to the pinned generator, avoiding
expansion of hundreds of megabytes of primitive arrays. Each leaf binds the
source hash, revision, dimensions and exact rank. Its accessor is compared
coefficient-for-coefficient with the original serializer at small sizes and
against frozen larger primitive data. These even-family identities retain
the published construction theorem as their proof basis. The odd size-31
leaf is frozen as sparse rational coefficients and separately passes a full
exact tensor check of all 887,503,681 coefficients, using signed 128-bit
accumulation with a verified overflow bound.

To replay the complete size-31 primitive check from the packaged coefficients
(requires a C++17 compiler):

```bash
PYTHONPATH=src python checks/replay_odd31.py
```
