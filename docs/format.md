# Scheme and coefficient format

`catalogue.json` and `schemes/<n>.json` identify each finite construction. Dimensions and ranks are decimal strings to preserve arbitrary-size integers in JavaScript and other readers. The field is Q. `root_object` is the SHA-256 of the uncompressed JSON object.

`src/cwschemes/data/index.json` maps object hashes to byte counts and stable logical resource names. The corresponding file is `data/objects/<hash>.gz`. Gzip is deterministic (`mtime=0`). The digest refers to uncompressed bytes. Identical objects are shared across scheme entries.

The root object contains one of the explicit construction formats interpreted by `_engine`: ordinary products, CW batches, mixtures, fixed-sector quotients, or compressed fixed-sector quotients, or explicit sparse schemes (`explicit-qcsr-square-v1`: a `source` path `seeds/explicit<n>_<rank>.npz` holding one CSR matrix per factor with int64 `indptr`, `indices` (row-major entry index), `numerators` and positive `denominators`, plus `metadata_json` with format `fmmp.qcsr.v1`; U and V read A and B, W writes C directly). Its recursive leaves are embedded constructions; references to external proof objects use `proofs/<sha256>.json`. A `source` reference without an extension denotes the same JSON resource. All such resources must occur in the index. No Python code is embedded in the objects.

**Z4-quotient schemes.** A `cw-z4-quotient-square-v1` root object has `n`, `full_n` = 2q³ (`n < full_n` means zero padding), `rank`, `outer_recipe` (always Strassen, `{"kind": "lille", "n": 2, "rank": 7}`) and `unit`. The unit is `{q, m: 3, rank, support, relations}`. Composite entries are ordinary `cw-mixed-square-v1` objects whose `packets` carry the same unit fields plus `"compression": "z4-quotient-v1"`; each such packet covers 7 consecutive outer products. `support` and `relations` are logical paths `z4q/support-<sha16>.json` and `z4q/relations-<sha16>.json`, shared by every size.

A `cw-z4-real-character-square-v1` object has the same fields and requires `unit.input_character = "odd"`. Both input row embeddings carry the character χ(g) = −1, while the output embedding is even. The reader excludes generator-fixed column orbits and retains the character in mixed certificates. Existing units without this field use the even character and their historical term ordering.

A `cw-z4-gauged-square-v1` object wraps an unchanged odd-character unit with an explicit permutation of the seven Strassen calls and two lists of nonzero rational input gauges. The output gauge is their reciprocal product. `deleted_terms` lists sorted original unit indices whose factors vanish after this square embedding, and `deletion_proof` names a hash-bound resource. Retained square term indices skip these entries. `audit_deletions()` independently rebuilds every declared zero factor from the public coefficient tables; the direct-sum unit's rank is unchanged.

- The support object stores `N = 9`, the Z4 generator `gens` and the 28 `rows` of zero masks (A, B, C). Rows are grouped into Z4-orbits in stored order; row orbit h is outer product h.
- The relation object has format `cw-z4-quotient-relations-v1`. It stores the mode permutation `perm` (packet mode j reads original mode `perm[j]`), the permuted `rows` and two `steps`. Each step lists its `pair`, `kept`, `deleted` and `zero_patterns` pattern codes (base-3 digits, 0 = pure, 1 = ordinary, 2 = special), plus per-block rational `basis`/`deleted`/`coefficients`. The historical object records its replayed q values in a `certified_q` list; new finite-q objects have a single positive integer `certified_q`. The requested q must be certified. Optional `rank_polynomials` stores the ten descending coefficients of `4 R_even(q)` and `4 R_odd(q)`; the audit reconstructs both from the retained patterns. These count polynomials do not extend coefficient validity beyond `certified_q`. `z4_quotient.replay_document(path, q)` instantiates the object as a `cw-compact-orbit-relations-v1` proof at an allowed q.
- Unit terms are Z4-orbits of columns (pattern, labels) of the compressed packet, ordered by pattern-orbit representative and then by an explicit unranking of label orbits (mixed radix; unordered pairs for stabiliser {e, g²}; necklaces for Z4). An input coefficient is the representative column's factor summed over the 4 rows of the product's row orbit. An output coefficient sums the factors of all members of the column orbit at the representative row.

These internal proof formats retain exact JSON integer literals and rational strings from the mathematical certificates. A non-Python implementation must parse those integer literals without floating-point rounding. They are not all bounded by 2^53. Rational strings are integers or `numerator/denominator`, never approximate decimals.

The reader verifies object hashes and allowed resource paths before unpacking dependencies into its cache. Coefficient tables generated later are disposable caches, not proof inputs. Delete the cache to reproduce a read from the bundled objects.

The `exact-rectangular-recipe-v1` square format stores a graph with a `root`
node, `nodes`, `tensor: [n,n,n]`, `rank` and `coefficient_field: "Q"`.
Node operations are `qcsr`, `classical`, `product`, `split`, `pad`, `permute`
and `uneven_strassen`. Frozen `qcsr` primitives use logical paths
`rectbases/<sha256>.npz` in the same object index, with one sparse CSR matrix
per factor. Rectangular dimensions refer to row, contraction and column
indices; factor edges are `(0,1)`, `(1,2)` and `(0,2)`. See
[rectangular.md](rectangular.md) for the construction and verification scope.

For a catalogue entry `s`, `s.coefficient("U", t, i, j)`, and similarly V or W, returns a `fractions.Fraction`. Indices are zero-based; `0 <= t < s.rank` and `0 <= i,j < s.n`. W describes the output directly, without transposing C. A requested batch is bounded by `max_queries`; full-array export is intentionally not automatic.

The inner data formats are documented by the reference coefficient implementations and the construction note. They are versioned by their `format` tags. Unknown versions must fail rather than silently change the decomposition.
