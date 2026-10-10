# Small rational rectangular compositions

The size-106 and size-146 schemes use **424,690** and **1,044,715** products.
They improve the strengthened rational composition reference counts 431,011
and 1,063,741 by 1.47% and 1.79%. Size 432 now uses **19,665,128** products,
through the routine product of the current published LITA size-18 and size-24
schemes, with ranks 3,031 and 6,488. This refresh is distinct from the new
uneven rectangular compositions.

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
