# What verification establishes

There are three distinct kinds of evidence in this release.

1. **Integrity and structure.** Object hashes, dependency resolution, catalogue dimensions and multiplication counts must agree. A changed object fails before use.
2. **Mathematical certificate checks.** The retained CW sectors must form the required induced matching; local rational completions are checked exactly. With `--exact`, each reduction relation is reconstructed and checked against its pair Gram matrix. These checks use the documented base identities and the imported published primitive schemes.
3. **Regression tests.** Exact coefficient queries, independently evaluated small multiplication examples and malformed-input tests check the software interface. They do not by themselves prove an enormous target tensor identity.

The two graded CW entries additionally use complete integer grading checks,
exact balanced-tuple counts, finite rational interpolation, and charged
outer-call batching. Their public verifier replays that finite construction
and the hash-bound outer metadata, without evaluating the enormous composed
coefficient sums. `--exact` has the same scope for these entries. Independent
capacity counting and small full tensor controls are recorded under
[checks/graded](../checks/graded/); see [the construction and its limits](graded_cw.md).

No target matrix with billions of rows is expanded. The full algorithm follows by induction through the explicit composition tree, using tensor products, disjoint batches and restriction by zero padding. Fixed-sector schemes additionally use the finite involution projection and the specified eigenspace decomposition. Their formulas and wiring are part of the reference implementation and construction argument; this is not a formal proof-assistant certificate.

## Reduction replay

For compact and refined orbit proofs, the checker reconstructs rational factor intertwiners, computes pair Gram matrices, and verifies every supplied relation `G[:, B] C = G[:, J]`. It checks the kept/deleted/zero partition and updates the free factor with the supplied exact C. A stored `complete` flag is not used as a substitute for those equalities.

For quotient proofs, the checker reconstructs the representation geometry. Original pair Grams are recovered with integer CRT and a proved absolute entry bound. Successive pair Grams incorporate previously verified factor changes; their clearing denominator and magnitude bounds are recomputed. The modulus exceeds twice the absolute bound. Final zero-factor proofs rebuild the small exact norm calculations. This is exact reconstruction, not a few random-prime tests.

The source uses assertions for mathematical invariants and refuses verification with optimized Python (`python -O`). `--threads` bounds the internal parallel arithmetic; the default is one. Large refined states and successive Grams can require substantial time and RAM. The quick default and exact mode therefore have different labels.

## Published primitives

The catalogue explicitly depends on the cited primitive identities, including Strassen, the rank-48 size-four scheme, small Lille schemes and the original pinned LITA coefficient data/formula. The checker does not re-prove the published general LITA family theorem symbolically. The included primitive coefficients are the exact data used by the constructions; refreshed comparison-only LITA ranks are not silently substituted into a scheme.

## Release evidence

The [version 0.1.0 results](../checks/summary.json) cover the 33 CW-based entries of that version; the explicit size-14, 16 and 32 entries and the Z4-quotient entries were added later (see below), and the Z4-quotient entries replace 15 of the version-0.1.0 entries. Selected exact-mode checks passed for sizes 432, 1,370 and 5,873,299,070, taking approximately 250, 376 and 280 seconds respectively on one Atlas core each. The first two exercise quotient and refined-orbit compression; the last exercises the uncompressed power-21 construction and its strict sub-2.7 threshold. These timings include loading and coefficient evaluation.

The checker rejected a deliberately altered rational relation. Six regression tests passed against the installed wheel from outside the repository, as did the small exact multiplication example and a fresh-cache size-1,024 smoke check. Forty-five nonzero coefficient comparisons checked the formula-based size-86 LITA reader against the original coefficient data. This last check is a regression comparison, not an exhaustive tensor-identity proof.

Machine-readable results under `checks/` record the checks actually run on the standalone package. The final release summary distinguishes all-scheme integrity/support/coefficient checks, selected complete reduction replays, tests of rejection, and the wheel installation test. A full `--exact --all` run is not implied unless explicitly recorded there. Prior research verification does not replace the release's fresh-install checks.

To reproduce without the original project, install the release, choose an empty `CWSCHEMES_CACHE`, change to another directory and run `cwschemes verify --all`. Exact replays use the same offline object bundle. No search executable, researcher home directory, cluster credential or hidden coefficient cache is needed.

The small rational rectangular update checks every graph node's dimensions
and count and every frozen coefficient object's hash. LITA primitive identities
use the pinned published construction theorem and complete coefficient
regeneration; the seven imported rectangular tensors and small uneven wiring
examples were checked completely over Q. `--exact` does not expand the large
LITA primitive tensors. See [rectangular.md](rectangular.md) and
`checks/small1024_rectangular*.json` for this scope.

## Explicit sparse schemes (sizes 14, 16, 32)

Entries in format `explicit-qcsr-square-v1` bundle every coefficient. The default check evaluates the algorithm on two random
input pairs modulo 2³¹−1 (Schwartz–Zippel). `--exact` checks every one of the n⁶ tensor coefficients modulo primes whose product
exceeds 2Δ(S+1), where Δ is the product of the per-factor denominator lcms and S = Σₜ‖uₜ‖∞‖vₜ‖∞‖wₜ‖∞; this certifies equality over Q.
Results are recorded in [checks/explicit.jsonl](../checks/explicit.jsonl). A single perturbed coefficient makes both checks fail.

## Z4-quotient schemes (sizes 1,024 to 65,536)

`cwschemes verify <n>` walks the entry's composition tree and checks each distinct unit (including its q and input character) in three ways.

- **Audit.** It checks the support structure: 28 CW rows in 7 free Z4-orbits, and the sunflower kernels with their Möbius coefficients. It checks the relation object: the kept, zero and deleted pattern sets partition each step and are Z4-stable, and relation steps modify only the output mode. The rank is computed three independent ways: orbit unranking, Burnside's lemma over the kept patterns, and the closed-form polynomial.
- **Exact packet identity.** It evaluates the packet tensor Σ_columns u(x) v(y) w(z) exactly over Q, with Fractions, at four sampled entries: a valid triple, a sunflower "garbage" block triple where the corrections must cancel, a label-perturbed diagonal entry and a uniform entry. Each value is compared with the 28 copies of ⟨q³⟩. The label sums factor coordinate-wise; packet modes 1 and 2 use the generator's own coefficient tables.
- **Orbit lemma.** Every member of a term's column orbit has the input factor prescribed by its character: equal for the even character and differing by the exact character sign for the odd one.

`--exact` additionally instantiates the relation object at each unit's q. A finite-q object rejects any q other than its saved `certified_q`; legacy objects retain their documented replay scope. It replays both steps exactly over Q with the same checker as other compact-orbit proofs (`replay_orbit`): zero patterns must be zero, and every relation G[:, B] C = G[:, J] must hold for the reconstructed pair Gram. One historical replay takes 8 to 9 minutes and about 3.2 GB of memory on one core.

The quotient step itself is covered by the orbit lemma above and the construction argument in construction.md, not by the replay. The same holds for the Strassen or outer wiring.

Release evidence is recorded in [checks/z4_quotient.json](../checks/z4_quotient.json) and [checks/z4_verify.jsonl](../checks/z4_verify.jsonl):

- **Faithful port.** At q = 2 (n = 16, rank 15,866), every one of the 3 × 15,866 × 256 coefficients equals the research generator's value. At q = 3, 8 and 17, 3,000 random coefficients each were compared the same way.
- **All coefficients at q = 2.** All n⁶ coefficients of the q = 2 instance (same support and relation objects) were checked modulo primes whose product exceeds 2Δ(S+1), which makes the check exact over Q. The vectorised factor matrices used for this check agree modulo 2³¹−1 with all published Fraction coefficients.
- **Random tensor entries.** Random full-tensor entries T[a,b,c] (uniform, valid, near-valid and garbage-targeted) were checked modulo 2³¹−1 and 2³¹−19 at q = 3, 8 and 17.
- **Exact packet entries and orbit lemma** were checked at q = 2, 3, 8 and 17.
- **Negative controls.** A +1 perturbation of the packet value, of the tensor value or of one factor coefficient makes each check fail.
- **Replay.** `cwschemes verify <n> --exact` passed for all 20 Z4 entries.

## October 2026 update

The new power-21 witnesses and recursive composition trees were audited independently of the discovery code. The independent pure-subset computation uses a Boolean superset transform; the generator uses downsets. Both mode-completion roles are checked. Exact local rational identities were checked at several q, and exhaustive small square examples cover padding, partially filled packets, leftover ordinary products and fully filled batches.

For every new threshold the integer inequality (n−1)^27 ≤ R^10 < n^27 is checked. Recursive audits verify the certified primitive rank, every inner dimension, each packet's exact rank, contiguous call allocation and the full square rank. Rational coefficient samples exercise the actual published loader. Relevant Z4 relation replays use the existing exact rational checker.

[checks/oct10.json](../checks/oct10.json) records changed objects and headline counts. Companion October check files record independent audits, package tests and fresh wheel loading. Earlier check files describe their historical object versions. No exhaustive expansion of a large target tensor is claimed.

## Small-square real-character update

At size 1,024 the odd-character Z4 unit removes 576 products. The generator-fixed input factors vanish over Q; all other orbit products use the same exact relation object and the documented signed row embedding. The uniform argument is in [real_character.md](real_character.md).

The q = 2 instance has rank 15,854 at size 16. Every one of its 16⁶ tensor coefficients passed an exact rational CRT check with a rigorously computed clearing denominator and entry bound. The port comparison checks every q = 2 input coefficient, every underlying rational coefficient table, every orbit unranking and 576 exact Fraction samples. All twelve excluded fixed columns and all thirty-four length-two stabilizer identities are checked. Larger dimensions follow from the uniform construction, rather than exhaustive expansion at size 1,024. The complete regression suite passed after the character and mixed-certificate changes.

Evidence is recorded in the `checks/small1024*` files. Earlier Z4 results refer to the even-character object versions and retain their historical counts.

The subsequent gauged square embedding removes 72 further products at q=8 without changing the odd direct-sum unit. The default public verifier replays every declared zero factor exactly over Q from the public packet tables, even when relation replay is not requested. The six q=2 deletions also pass a full rational CRT tensor check at rank 15,848; scalar queries agree with the independently assembled vectorized factors. The hash-bound proof and independent results are in `checks/small1024_gauge.json`. See [the square cancellation argument](real_character.md#further-square-only-cancellation).

## Further C4 rational compression

The current size-1,024 entry uses a new retained set and separate q=2/q=8
coefficient witnesses. Source rational replays check every original and
averaged relation; the public exact-mode check replays the transported q=8
object. The audit reconstructs the entire advertised count polynomial from
fixed patterns, and the reader rejects uses outside `certified_q`.

The independent q=2 full tensor gate covers the odd square with ordinary
Strassen wiring at rank 15,838. Its eight-prime CRT modulus exceeds the
rigorous denominator and entry bound, and every one of 16⁶ coefficients
agrees. This gate is not a full contraction of the gauged rank-15,832 variant.
Exact coordinate transport, complete rational MM2 wiring and separately
replayed physical zero factors certify the gauged construction. All 72 q=8
deletions bind to the new unit rank 198,679,264, giving 198,679,192 at size
1,024. See [c4_recompression.md](c4_recompression.md) and
`checks/small1024_recompression*.json` for the current object versions and
distinct check scopes. Historical evidence above retains its historical ranks.

## Strassen-tensor (S_q) packets

Packets with `"local_tensor": "strassen-completion-v1"` (sizes 8,388,608, 13,468,840,704 and 13,947,137,604) are checked by `cwschemes verify <n>` as follows; `--exact` adds nothing, since these packets have no relation proofs.

- **Support audit**, at load and again from the root object. Every row partitions the N = 3m coordinates into masks of size m, uses only the letters its coordinate kinds allow, and has distinct blocks in each leg. A cross triple must have B_v on kind 0 equal to the complement of A_u there, so grouping rows by that part enumerates **every** triple (u,v,w) with A_u, B_v, C_w partitioning the coordinates. The audit requires exactly the M diagonal triples.
- **Exact local identity.** From the generator's own local factors, all three kinds are rebuilt exactly at the packet's q and compared with Str_q plus Σₜ aₜbₜcₜ. This also checks that no local block has two zeros.
- **Counts.** The packet rank (q+1)^N, the batch allocation and the square rank are recomputed at load; unused rows of the final packet receive zero input.
- **Exact packet entries.** At sampled entries (a valid entry of a row, a changed contracted label and random block triples), Σ over all (q+1)^N terms of U·V·W is evaluated exactly. The term set is the product of local choices, so the sum factors coordinate-wise.

The full composition then follows from the documented batching, as for the CW packets. The default coefficient queries at the last matrix entry are zero for these entries; nonzero rational probes are recorded separately.

Release evidence for these entries is in `checks/sq/`:

- `export.json`: input hashes of the study certificates, rows rebuilt from the RS witnesses, the package audit, a separate O(M²) cross-triple enumeration, the local identity and every count. `rs_sets.json` holds the RS witnesses.
- `independent.json`: a separate support checker on both supports, an end-to-end exact evaluation of all 2²⁷ terms of the power-27 packet at q = 1 (every output block of the 9,842 rows), and rejection of four corrupted supports. It also records the earlier independent study checks, which were not repeated here: the full power-18 packet at q = 2 (3¹⁸ terms) and an exact CRT replay of the size-36 current-LITA tensor over Q.
- `lita36_index.json`: the indexed current-LITA accessor used as the outer scheme agrees with the pinned serializer on all 19,666 terms.
- `probes.json`: 24 nonzero coefficients through the published loader, each equal to a separate evaluation of the row wiring, local factors and outer coefficients.
- `verify.jsonl` and `tests.txt`: the default verifier output and the new regression tests (complete tensors at N = 6 and N = 9, a complete size-8 composed square, negative controls).

Exactness for every q follows from the uniform local identity together with the support audit, which does not depend on q. No target tensor of these sizes is expanded.

## November 2026 updates

- **Below 1,024 (`checks/nov16_below1024.json`).** All 158 reachable non-LITA `qcsr` leaves passed multi-prime C++ CRT verification over Q exceeding `2ΔS`; all `641,243,628` search DAG nodes passed independent dimension/rank/acyclicity audits; every `serendipitous` bud decomposition is verified over Q against the exact Brent equations at load time; and all 217 packaged square certificates (`81 .. 1,010`) passed dense end-to-end random matrix multiplication modulo `2^31 - 1` with exact scalar multiplication count checks.
- **Coordinate-8 refined Z4 units and Pair-Schur towers (`checks/nov16_above1024.json`).** Relation objects with `"coordinates": [8]` are bound to `certified_q` in `{3, 8, 9, ..., 20}` and replayed over Q by `--exact`; the `q^2 + q` gauged square deletions are rebuilt and verified to vanish over Q by default. At `q = 3` (`n = 54`, rank 203,307), an independent audit verified 3-prime CRT random matrix products and all `54^6 = 24,794,911,296` tensor entries modulo `1,048,573`. Pure-integer checks certify `(n - 1)^27 <= R^10 < n^27` at `n = 240,619,650` and `(n - 1)^53 <= R^20 < n^53` at the `5.85 × 10^64` graded entry.
- **Zero-trimmed S_q packets (`zero-trimmed-v1`).** `cwschemes verify <n>` audits the rebuilt RS supports (`N = 18, 27, 36`), verifies the retained mask count against the exact 2D downsets `(D_A, D_B, D_C)`, checks that discarded special-choice masks vanish, and evaluates sampled exact packet entries over Q. Complete small tensor identities (`m = 2, 3`) and partial-batch `MixedSquare` expansions are covered in `tests/test_nov16_improvements.py`.
