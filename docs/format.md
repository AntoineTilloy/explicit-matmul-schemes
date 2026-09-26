# Scheme and coefficient format

`catalogue.json` and `schemes/<n>.json` identify each finite construction. Dimensions and ranks are decimal strings to preserve arbitrary-size integers in JavaScript and other readers. The field is Q. `root_object` is the SHA-256 of the uncompressed JSON object.

`src/cwschemes/data/index.json` maps object hashes to byte counts and stable logical resource names. The corresponding file is `data/objects/<hash>.gz`. Gzip is deterministic (`mtime=0`). The digest refers to uncompressed bytes. Identical objects are shared across scheme entries.

The root object contains one of the explicit construction formats interpreted by `_engine`: ordinary products, CW batches, mixtures, fixed-sector quotients, or compressed fixed-sector quotients, or explicit sparse schemes (`explicit-qcsr-square-v1`: a `source` path `seeds/explicit<n>_<rank>.npz` holding one CSR matrix per factor with int64 `indptr`, `indices` (row-major entry index), `numerators` and positive `denominators`, plus `metadata_json` with format `fmmp.qcsr.v1`; U and V read A and B, W writes C directly). Its recursive leaves are embedded constructions; references to external proof objects use `proofs/<sha256>.json`. A `source` reference without an extension denotes the same JSON resource. All such resources must occur in the index. No Python code is embedded in the objects.

These internal proof formats retain exact JSON integer literals and rational strings from the mathematical certificates. A non-Python implementation must parse those integer literals without floating-point rounding. They are not all bounded by 2^53. Rational strings are integers or `numerator/denominator`, never approximate decimals.

The reader verifies object hashes and allowed resource paths before unpacking dependencies into its cache. Coefficient tables generated later are disposable caches, not proof inputs. Delete the cache to reproduce a read from the bundled objects.

For a catalogue entry `s`, `s.coefficient("U", t, i, j)`, and similarly V or W, returns a `fractions.Fraction`. Indices are zero-based; `0 <= t < s.rank` and `0 <= i,j < s.n`. W describes the output directly, without transposing C. A requested batch is bounded by `max_queries`; full-array export is intentionally not automatic.

The inner data formats are documented by the reference coefficient implementations and the construction note. They are versioned by their `format` tags. Unknown versions must fail rather than silently change the decomposition.
