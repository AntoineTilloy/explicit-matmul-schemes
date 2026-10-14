# Explicit matrix multiplication schemes

A catalogue of **99 explicit schemes for square matrix multiplication over R and C**, with rational coefficients, exact multiplication counts, compact coefficient descriptions, and finite proof witnesses.

A new finite graded CW construction reaches exponent **2.649062317** at
**n = 2⁴⁸ · 6⁷⁰ ≈ 8.32 × 10⁶⁸**, using **4 · 31,681 · 7²¹⁰** products.
Its companion cropped entry crosses 2.65 at approximately **7.86 × 10⁶⁸**.
Both satisfy the strict integer check **R²⁰ < n⁵³** and work recursively on
matrix blocks. These are enormous algebraic constructions; the size-1,024
count remains **198,679,192**. See the [finite proof, coefficient definition,
and verification limits](docs/graded_cw.md). This is the smallest crossing
found in the bounded search, without a global minimality claim.

Most schemes use Coppersmith–Winograd (CW) tensor restrictions, rational compression and composition with published algorithms. Three small explicit schemes, at sizes 14, 16 and 32, use one product fewer than the pinned LITA comparison: **1,593**, **2,236** and **14,196**. Current upstream LITA also attains these counts. Every coefficient is accessible as an exact rational number.

Between sizes 1,024 and 65,536, the best entries come from a Z4-symmetric CW packet quotiented to 7 products and fed to Strassen's algorithm. Size **9,826** needs **74,352,484,826** multiplications (effective exponent **2.723013**), and size **1,024** now needs **198,679,192**. A new rational compression saves **4,096** products against the preceding entry; together with the odd input character and square cancellations, the saving against the original even-character entry is **4,744**. This is 1.43% fewer than 14,197² = 201,554,809 obtained by squaring the pinned LITA size-32 scheme. The best exponent in the earlier size range is **2.692242** at size **12,310,160,190**, following the power-21 support and recursive composition improvements of October 2026.

![Effective exponents of our saved schemes: detail and wider size range](figures/exponents.svg)

The three panels preserve detail through size 16,384, show the earlier range
through 30 billion, and show the full range including the new graded schemes.
All use logarithmic size axes. The main figure contains only rational schemes
valid over R and C; characteristic-two results remain separate.

The exponent of a finite scheme is **logₙ R**, where R counts bilinear scalar multiplications. Additions and multiplications by fixed constants are excluded. These are algebraic reference constructions; the counts do not establish a practical speedup over numerical matrix multiplication libraries.

The dashed reference is the current size-32 LITA recursion exponent,
**log₃₂(14,196) = 2.758639**. Pure recursion attains this value at powers of
32, not at every finite size: blocking and padding must be charged explicitly.
The new rectangular points improve counts at their own sizes while remaining
above this exponent. At size 1,024, current LITA squared costs **201,526,416**;
our existing **198,679,192** construction saves **1.41%**. The rectangular
continuation does not improve that result. See the
[same-size blocking comparison](docs/extended_rectangular.md#comparison-with-size-32-lita-blocking).

## Quick start

Requires Python 3.11 or newer. From a checkout:

```bash
python -m pip install .
cwschemes info 16384
```

```python
import cwschemes

s = cwschemes.load(16384)
assert s.n == 16384
assert s.rank == 317774199680
u = s.coefficient("U", term=0, row=0, col=0)  # fractions.Fraction
print(u)
```

No search needs to be rerun. All required scheme data is included, compressed and hash-checked. Resources are unpacked on demand into `~/.cache/cwschemes`; set `CWSCHEMES_CACHE` to choose another directory. The full coefficient arrays are never allocated. For the enormous graded entries, coefficients have an explicit finite slot-sum definition, whose evaluation is generally computationally prohibitive. Their default coefficient query guard prevents accidentally starting such a sum; see [coefficient access](docs/graded_cw.md#ordered-blocks-and-coefficient-access).

First use of a construction can take minutes to build exact intermediate tables. Subsequent queries reuse those tables. See the [release checks](checks/summary.json) for measured verification results.

## Catalogue

| Square size | Multiplications | Effective exponent |
|---:|---:|---:|
| [14](schemes/14.json) | 1,593 | 2.793942453 |
| [16](schemes/16.json) | 2,236 | 2.781676118 |
| [32](schemes/32.json) | 14,196 | 2.758639372 |
| [686](schemes/686.json) | 70,153,048 | 2.766272903 |
| [1,024](schemes/1024.json) | 198,679,192 | 2.756586554 |
| [2,025](schemes/2025.json) | 1,303,338,368 | 2.756771167 |
| [2,048](schemes/2048.json) | 1,303,338,368 | 2.752687685 |
| [8,192](schemes/8192.json) | 45,396,314,240 | 2.723219701 |
| [9,826](schemes/9826.json) | 74,352,484,826 | 2.723013418 |
| [16,384](schemes/16384.json) | 317,774,199,680 | 2.729229360 |
| [65,536](schemes/65536.json) | 14,360,405,808,994 | 2.731697610 |
| [528,108,169](schemes/528108169.json) | 355,919,505,257,948,041,899,273 | < 2.7, exactly checked |
| [536,870,912](schemes/536870912.json) | 355,919,505,257,948,041,899,273 | 2.697789558 |
| [12,310,160,190](schemes/12310160190.json) | 1,463,616,394,485,527,559,895,695,360 | 2.692241557 |
| [13,468,840,704](schemes/13468840704.json) | 2,014,523,380,268,092,193,411,235,280 | 2.695555553 |
| [≈ 8.32 × 10⁶⁸](docs/graded_cw.md) | 4 · 31,681 · 7²¹⁰ | 2.649062317 |
| [≈ 7.86 × 10⁶⁸ (cropped)](docs/graded_cw.md) | 4 · 31,681 · 7²¹⁰ | < 2.65, exactly checked |

See [catalogue.json](catalogue.json) for all sizes, or run `cwschemes list`.

Rectangular compositions at sizes **106** and **146** use **424,690**
and **1,040,656** products, saving **1.47%** and **2.17%** against the strengthened
rational composition reference. Size **432** now uses **19,665,128**, through
a routine refresh to the current published LITA primitives. All three have
exact rational coefficients in the installed package; see the
[rectangular construction and check scope](docs/rectangular.md).

The first search applying unequal Strassen splits inside intermediate rectangles,
then feeding the improved costs back into composition, added **12 square schemes
at sizes 125–242**. At that stage the largest relative gain was at **194**, from **2,322,311** to
**2,309,010** products (about **0.57%**). These comparisons use the full previous
strengthened rational envelope. The size-1,024 count remains **198,679,192**.
All new schemes are bilinear and can be substituted recursively on matrix blocks.

Extending the intermediate rectangular table to dimension 256 improves **48
square counts at sizes 137–370**, including ten previous catalogue entries.
Size **291** falls from **7,273,834** to **7,144,555** products (**1.78%**).
Some other gains are only one product; the full counts, comparison baselines,
and exact verification scope are in [the extended search record](docs/extended_rectangular.md).
The coefficient engine is unchanged, and size 1,024 remains at **198,679,192**.

![Verified improvements over R and C through size 1,024](figures/small1024.svg)

The first sub-2.7 size **in this catalogue** is now **528,108,169**, padded from a natural size-536,870,912 construction. Its rank R = 355,919,505,257,948,041,899,273 satisfies

```text
528,108,168^27 <= R^10 < 528,108,169^27.
```

This is the strict threshold for this particular scheme. It is **91.01% smaller** than the previous catalogue threshold of 5,873,299,070. The natural construction uses an outer size-32 scheme with 14,196 calls, 9 power-18 packets of 1,569 independent products each, and 75 explicitly specified recursive remainders. The best exponent and the threshold are different catalogue points. A sub-2.7 construction at size 1,024 remains open.

## Z4-quotient schemes (sizes 1,024 to 65,536)

These entries (formats `cw-z4-quotient-square-v1`, `cw-z4-real-character-square-v1` and `cw-z4-gauged-square-v1`, or `z4-quotient-v1` units inside a mixed outer scheme) share one CW support. The larger entries retain the original relation certificate; size 1,024 uses a new certificate bound to q=8.

- **Packet.** A power-9 CW support with 28 rows is invariant under a cyclic group Z4 permuting the 9 coordinates. Its rows form 7 free Z4-orbits. The packet tensor is the local completion G_q^{⊗9} restricted to the support, plus Möbius corrections on the 36 sunflower kernels; it equals 28 independent ⟨q³⟩ products.
- **Compression.** In the joint triangular representation, patterns whose factor vanishes on the support are removed. A Z4-invariant pair Gram elimination then deletes 300 more patterns, and its relations are absorbed into the output factor. Both steps delete whole Z4-orbits and leave the two input factors equivariant.
- **Quotient.** On Z4-invariant inputs, all columns of one orbit give the same product. One multiplication per column orbit, whose output factor sums over the orbit, computes the 7 row-orbit products independently. The unit rank is

  ```text
  R7(q) = (q^9 + 17q^8 + 122q^7 + 392q^6 + 493q^5 + 131q^4 + 2q^3 + 2q^2) / 4.
  ```

- **Sizes.** Strassen's 7 products give ⟨2q³⟩ with R7(q) multiplications for q = 8, …, 20; sizes 4,096, 5,324 and 13,500 are padded from these. Sizes 10,976 and 16,384 place 7 units inside Strassen squared. Size 32,768 places 48 units inside ⟨8;336⟩. Size 65,536 is padded from 24·14³ = 65,856 and uses 929 units inside the original rank-6,503 size-24 scheme, with no remainder products.

The original relation data does not depend on q. It was replayed exactly over Q for every q from 2 to 25 before inclusion; `cwschemes verify <n> --exact` replays it again at the q in use. The new size-1,024 certificate has separate exact q=2 and q=8 witnesses and rejects other q values.

The size-1,024 entry uses the character χ(g) = −1 for a generator of Z4 on both inputs and the even output character. Columns fixed by g have zero input factors over Q, removing q³ + q² products. The new compression deletes **304** patterns rather than 300, changing the q⁴ coefficient in the displayed numerator from **131 to 127**. At q=8 the odd unit has rank **198,679,264**. A permuted and signed Strassen embedding makes another **72** input factors vanish, lowering the square count to **198,679,192**. See [the new compression and its finite-q proof scope](docs/c4_recompression.md) and [the real-character argument](docs/real_character.md). These arguments give rational schemes; characteristic-2 results are reported separately.

The October 2026 update rebuilds the older large recipes with the current certified leaves and adds two independently audited power-21 supports. Recursive recipes retain the exact primitive versions used by their coefficient generators. The update and its validation scope are recorded in [checks/oct10.json](checks/oct10.json) and [verification.md](docs/verification.md).

## Explicit schemes at sizes 14, 16 and 32

These three entries store every coefficient explicitly (sparse rational factor matrices, format `explicit-qcsr-square-v1`). They are points of the exact (γ, r₀) parameter family of the [LITA](https://github.com/khoruzhii/lita) row–column aggregation scheme. Two of LITA's three parameter views are unchanged; the third is re-anchored so that one more raw term vanishes (at size 32: 14,314 raw terms − 104 vanishing − 14 merged = 14,196). The same choice gives one product fewer than the pinned LITA counts 1,594, 2,237 and 14,197. The newer upstream revision [472a79f](https://github.com/khoruzhii/lita/tree/472a79fbc8248423a97bebb3aacbc89e769b9b3e) also attains the saved ranks; these three entries are not claimed as improvements over that revision.

`cwschemes verify 32` evaluates the algorithm on random matrices modulo 2³¹−1; `cwschemes verify 32 --exact` checks every tensor coefficient modulo enough primes to certify equality over Q (a few minutes of CPU at size 32).

## Verify

```bash
# Integrity, exact support/completion identities, counts and coefficient queries:
cwschemes verify 16384
cwschemes verify --all

# Reconstruct and check the rational compression relations as well:
cwschemes verify 432 --exact --threads 1
cwschemes verify 9826 --exact

# Verify the new threshold and its compressed recursive leaf:
cwschemes verify 528108169 --exact
```

The default checks are not advertised as a full reduction proof. `--exact` additionally rebuilds the relevant Gram matrices and checks the supplied rational relations, or checks exact zero-factor norms. It can take substantially longer and use more memory, especially for refined orbit proofs. Verification must run without Python's `-O` flag.

The mathematical basis consists of the documented CW identities, group projections, composition rules and published primitive schemes. This is an arithmetic certificate checker, not a formal proof-assistant development. Read [verification.md](docs/verification.md) for the precise scope and the release checks actually performed.

## Literature comparison

The main figure displays our saved R/C schemes; the three explicit small schemes are drawn as diamonds, the Z4-quotient schemes as squares, and the rectangular compositions as triangles. All 97 saved points beat their stated pinned comparison in [counts.csv](comparisons/counts.csv). The catalogue does not claim globally optimal ranks or independently established records.

- Up to 16,384, the comparison data uses the same product, padding and scalar-peeling closure for published constructions alone and for published constructions augmented with our saved schemes.
- At larger sizes, each comparison is one explicit published construction, not an exhaustive literature optimum.

The comparison pins [LITA](https://github.com/khoruzhii/lita/tree/c1dd9225df98676e385b53ae7517ff2ea0ec5779) at rank **14,197** for size 32 and includes its even and applicable odd families; the explicit size-14, 16 and 32 entries are one below these pinned counts. Other sources are the [Lille catalogue](https://fmm.univ-lille.fr/) and the documented [Schwartz–Zwecher recompression](https://arxiv.org/abs/2508.01748). Existing scheme certificates retain the older primitive versions with which they were constructed.

The saved points have [exact comparison recipes](comparisons/recipes.json). [Source revisions](comparisons/sources.json), [all point counts](comparisons/counts.csv), and the [dense finite closure](comparisons/closure.csv) are included.

To rebuild the figure from these pinned inputs:

```bash
python -m pip install '.[plot]'
python comparisons/build.py
python comparisons/small1024.py
```

Characteristic-two experiments are secondary leads toward faster rational schemes.
Their counts, certificates and separate figure are in
[characteristic2/](characteristic2/README.md); the best saving here is **2.43%**
at size 243. No corresponding R/C lift is claimed.

## Format and reuse

The public API uses the direct-output convention

```text
C[i,j] = Σₜ W[t,i,j] (Σₐᵦ U[t,a,b] A[a,b]) (Σ꜀𝒹 V[t,c,d] B[c,d]).
```

Scheme descriptors point to immutable compressed objects. Shared proof objects and recursive leaves are stored once; no files from the discovery repository are required. See [format.md](docs/format.md), [construction.md](docs/construction.md) and the [small runnable example](examples/strassen.py).

The original search programs, cluster launchers, logs and experimental dependencies are deliberately absent. Third-party source and seed provenance is documented in [third_party.md](docs/third_party.md). See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff) for reuse and citation.
