# Exact finite graded CW constructions

These entries cross effective exponent 2.65 with rational bilinear schemes.
They use the explicitly saved puzzle rows, integer potentials, and finite
composition recipes in their proof objects. They are valid recursively when
entries are themselves matrices. Their very large dimensions do not improve
the existing size-1,024 entry, which still costs 198,679,192 multiplications.

The puzzle examples were transcribed from the work of
[Anderson and Le on simplifiable strong uniquely solvable puzzles](https://arxiv.org/abs/2307.06463).
We verify the concrete rows and the complete tensor support directly; the
finite counts below do not follow merely by extrapolating an asymptotic exponent.
This is a new construction in this catalogue, without a claim of a globally
smallest possible crossing or a practical numerical speedup.

## Rational local tensor and exact grading

For an integer q at least 2, put s = e₁ + … + e_q and

```
L_i = e_0 + e_i + s/q,                 i = 1,...,q,
L_* = e_0 + 3s/(2q),
G_q = sum_i L_i ⊗ L_i ⊗ L_i - (4q/3) L_* ⊗ L_* ⊗ L_*.
```

This is an explicit rational decomposition with q+1 products. A coefficient
with exactly one zero index is the equality delta on the other two indices;
with two zero indices it vanishes; at (0,0,0) it equals −q/3. We retain all
sectors with no zero indices as potentially nonzero.

Each puzzle row is a word over {1,2,3} specifying which tensor mode has a
zero in each coordinate. Each mode's zero masks must be distinct. A triple
of masks A,B,C can survive only when A∩B = A∩C = B∩C. The checker enumerates
**every** such triple, including unwanted triples with three simultaneous
zeros. It checks integer potentials a_i,b_i,c_i satisfying

```
a_i + b_i + c_i = 0                       for each desired row,
1 <= a_i + b_j + c_k <= d                 for every unwanted surviving triple.
```

Scaling mask variables by powers of a nonzero rational parameter z gives a
polynomial tensor whose constant coefficient is exactly the wanted direct
sum. Individual factor entries can have negative powers, but every complete
tensor coefficient has degree between zero and d.

## Keeping all balanced products

Take t copies of each of the three cyclic orientations, for 3t copies in
total. If a row's zero counts are (a,b,c), let
H(x,y,z) = Σ_rows xᵃ yᵇ zᶜ. The number of retained square products is exactly

```
M = [x^(wt) y^(wt) z^(wt)] (H(x,y,z) H(y,z,x) H(z,x,y))^t,
```

where w is the puzzle width. Each product has side q^(wt). Rectangular
tuples that do not satisfy the balance condition are discarded. Distinct
selected tuples have disjoint variable sets in every mode by mask
injectivity.

The coefficient reader counts suffix completions by an integer dynamic
program and un-ranks a selected tuple lexicographically. Discovery used
carry-free integer encoding; an independent verifier recomputed capacity
by ordinary profile counting and a direct sum over profile pairs. Thus the
saved capacity has independent counting implementations.

The total degree is at most D = 3td. Constant-coefficient interpolation at
the nonzero nodes 1,…,D+1 uses weights
λ_j = (−1)^(j−1) binomial(D+1,j). This gives an exact finite packet cost

```
P = (3td + 1) (q + 1)^(3wt).
```

There is no limiting step and no uncharged interpolation overhead.

## Square composition and the strict threshold

Let an existing certified square scheme have side b and cost r. Assign
its ordered block products to the M packet slots, using as many packets as
needed. Unused input slots are zero. The resulting finite scheme has

```
n = b q^(wt),             R = ceil(r/M) P.
```

The outer scheme is an explicitly recorded tensor product of a saved
catalogue base, the rank-14,196 size-32 scheme, and Strassen's rank-seven
size-two scheme. References are immutable, hash-checked proof objects.
Every packet, including a partially filled final packet, is charged in full.

The claim logₙ R < 2.65 is verified by the strict integer comparison
**R²⁰ < n⁵³**. A companion entry pads inputs from the smallest valid size for
that same rank ceiling, namely floor((R²⁰)^(1/53))+1. Its printed floating
exponent can round to 2.65; the integer comparison remains strict. This
cropping minimum concerns that particular rank ceiling, not all schemes.

## Ordered blocks and coefficient access

In each constituent, the first mode is indexed (I,J), the second (J,K),
and the third tensor mode (K,I) supplies the coefficient of output (I,K).
The third index convention never transposes matrix-block values. Each
product has a linear form in A on the left and a linear form in B on the
right, with rational scalar coefficients. The identities therefore hold
over any associative algebra over Q, including square matrix blocks.

The source defines every coefficient as an explicit finite sum over outer
calls, without allocating a full factor matrix. **Evaluating such a sum at
the new enormous dimensions is generally computationally prohibitive.**
These entries do not provide fast random access to their coefficients.
The usual `coefficient` method refuses sums with more than 10,000 active
slots. Passing `max_slot_terms=None` explicitly requests the full sum;
`coefficient_summands` exposes its iterator. Even finding or evaluating
individual streamed summands can be expensive through the outer scheme.

`cwschemes verify N` replays the grading, exact balanced capacity, rational
construction counts, referenced base metadata, and strict threshold. It
does not attempt the immense coefficient sums. `--exact` uses that same
finite proof replay for these entries. The verification basis is the
documented tensor identity and previously certified outer schemes, rather
than expansion of the complete large tensor.

Release checks include complete rational tensor identities on small
instances, exhaustive selected-tuple unranking, all three coordinate
orientations, multibatch composition, and padding. They also check the actual
large recipe's counts and selected coordinate mappings. The finite search
uses bounded candidate families; its results establish constructions, not
a global lower bound on the size of a sub-2.65 algorithm.

## Saved crossing

Solving exact integer linear programs over the Aronov–Landsberg puzzles lowers
the certified positive grading degrees on width-6 `(14, 6)` (`204 -> 203`),
width-7 `(23, 7)` (`1,056 -> 436`, and `411` on a 22-row subset), width-8
`(36, 8)` (`950`), and width-9 `(55, 9)` (`1,874`). Pairing the 14-row,
width-six puzzle (`d = 203`, `t = 11`, `q = 6`, interpolation factor
`3 * 11 * 203 + 1 = 6,700`) with the new single-packet Strassen-tensor outer
at `n = 27,113,235,502,176 = 96 * 9^12` (exponent `2.679931`) requires three
packets and gives

```
n = 96 * 9^12 * 6^66
R = 3 * 6700 * 7^198
log_n(R) = 2.6490125772706076.
```

The natural side is approximately `6.18 × 10⁶⁴` (65 decimal digits); cropping
gives approximately `5.85 × 10⁶⁴` (`13,446.8×` smaller than the preceding
`7.86 × 10⁶⁸` crossing) under the same strict integer check `R^20 < n^53` and
minimality check `(n - 1)^53 <= R^20`. Complete decimal integers and SHA-256
proof hashes are in `catalogue.json` and
[checks/nov16_above1024.json](../checks/nov16_above1024.json).
