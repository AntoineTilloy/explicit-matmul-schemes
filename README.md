# Explicit matrix multiplication schemes

A catalogue of **33 explicit rational schemes for square matrix multiplication**, with exact multiplication counts, compact coefficient descriptions, and finite proof witnesses.

The schemes use Coppersmith–Winograd (CW) tensor restrictions, rational compression and composition with published algorithms. Every coefficient is accessible as an exact rational number. The catalogue includes a size-16,384 scheme with effective exponent **2.731191**, and a size-13,468,840,704 scheme with exponent **2.698590**.

![Effective exponents and pinned literature comparisons](figures/exponents.svg)

The exponent of a finite scheme is **logₙ R**, where R counts bilinear scalar multiplications. Additions and multiplications by fixed constants are excluded. These are algebraic reference constructions; the counts do not establish a practical speedup over numerical matrix multiplication libraries.

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
assert s.rank == 323882771968
u = s.coefficient("U", term=0, row=0, col=0)  # fractions.Fraction
print(u)
```

No search needs to be rerun. All required scheme data is included, compressed and hash-checked. Resources are unpacked on demand into `~/.cache/cwschemes`; set `CWSCHEMES_CACHE` to choose another directory. The full coefficient arrays are never allocated.

First use of a construction can take minutes to build exact intermediate tables. Subsequent queries reuse those tables. See the [release checks](checks/summary.json) for measured verification results.

## Catalogue

| Square size | Multiplications | Effective exponent |
|---:|---:|---:|
| [686](schemes/686.json) | 70,153,048 | 2.766272903 |
| [2,048](schemes/2048.json) | 1,303,338,368 | 2.752687685 |
| [16,384](schemes/16384.json) | 323,882,771,968 | 2.731191484 |
| [536,870,912](schemes/536870912.json) | 390,821,774,252,911,273,874,842 | 2.702443346 |
| [5,873,299,070](schemes/5873299070.json) | 237,674,020,365,443,824,643,488,557 | < 2.7, exactly checked |
| [13,468,840,704](schemes/13468840704.json) | 2,162,272,334,177,007,487,103,296,378 | 2.698590115 |

See [catalogue.json](catalogue.json) for all sizes, or run `cwschemes list`.

The first sub-2.7 size **in this catalogue** is obtained by padding to a natural size-5,905,580,032 construction. Its unchanged rank R satisfies

```text
5,873,299,069^27 <= R^10 < 5,873,299,070^27.
```

This establishes the strict threshold for that scheme, not the smallest possible matrix size. Its exponent at the natural size is 2.699342235. The new power-21 support computes 5,032 independent products; one packet and 123 recursive remainders supply the 5,155 calls of the original size-22 outer scheme.

## Verify

```bash
# Integrity, exact support/completion identities, counts and coefficient queries:
cwschemes verify 16384
cwschemes verify --all

# Reconstruct and check the rational compression relations as well:
cwschemes verify 432 --exact --threads 1

# An uncompressed large construction has no compression relations to replay:
cwschemes verify 5873299070 --exact
```

The default checks are not advertised as a full reduction proof. `--exact` additionally rebuilds the relevant Gram matrices and checks the supplied rational relations, or checks exact zero-factor norms. It can take substantially longer and use more memory, especially for refined orbit proofs. Verification must run without Python's `-O` flag.

The mathematical basis consists of the documented CW identities, group projections, composition rules and published primitive schemes. This is an arithmetic certificate checker, not a formal proof-assistant development. Read [verification.md](docs/verification.md) for the precise scope and the release checks actually performed.

## Literature comparison

The figure shows **best available bounds within the stated catalogue and constructions**, not globally optimal ranks or independently established record claims. The size-432 and size-1,024 candidates lose their finite-size comparisons and remain explicitly marked.

- Up to 16,384, both curves use the same product, padding and scalar-peeling closure. One starts with published constructions; the other also includes our saved schemes.
- At larger sizes, each diamond is one explicit published construction, not an exhaustive literature optimum. Sparse points are not interpolated into unproved bounds.
- The horizontal LITA32 line is an exponent reference, not a finite-size rank bound at arbitrary n.

The comparison pins [LITA](https://github.com/khoruzhii/lita/tree/c1dd9225df98676e385b53ae7517ff2ea0ec5779) at rank **14,197** for size 32 and includes its even and applicable odd families. Other sources are the [Lille catalogue](https://fmm.univ-lille.fr/) and the documented [Schwartz–Zwecher recompression](https://arxiv.org/abs/2508.01748). Existing scheme certificates retain the older primitive versions with which they were constructed.

The saved points have [exact comparison recipes](comparisons/recipes.json); the dense curves are reproduced by the comparison builder. [Source revisions](comparisons/sources.json), [all point counts](comparisons/counts.csv), and the [dense finite closure](comparisons/closure.csv) are included.

To rebuild the figure from these pinned inputs:

```bash
python -m pip install '.[plot]'
python comparisons/build.py
```

## Format and reuse

The public API uses the direct-output convention

```text
C[i,j] = Σₜ W[t,i,j] (Σₐᵦ U[t,a,b] A[a,b]) (Σ꜀𝒹 V[t,c,d] B[c,d]).
```

Scheme descriptors point to immutable compressed objects. Shared proof objects and recursive leaves are stored once; no files from the discovery repository are required. See [format.md](docs/format.md), [construction.md](docs/construction.md) and the [small runnable example](examples/strassen.py).

The original search programs, cluster launchers, logs and experimental dependencies are deliberately absent. Third-party source and seed provenance is documented in [third_party.md](docs/third_party.md). See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff) for reuse and citation.
