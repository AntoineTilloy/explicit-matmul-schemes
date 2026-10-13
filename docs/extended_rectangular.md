# Extended rectangular compositions through size 370

The 256-grid continuation improves 48 square counts at sizes 137–370, adding
38 catalogue entries and replacing ten previous entries. All coefficients
are rational, every product has an A-linear left factor and a B-linear right
factor, and the schemes permit recursive substitution by matrix blocks.
Size 1024 remains at **198,679,192** products; sub-2.7 at that size is open.

The comparison is the complete previous rational square envelope, including
all twelve improvements published in commit `b647d95`. It is stronger than
comparing only with the previously saved catalogue points. Several gains
are just one product and are reported as such.

## How the larger table helps

The search retains the earlier rectangular table through dimension 128 and
extends each axis to 256. Newly introduced rectangles start from classical
multiplication. It iterates padding, axis splits, exact-divisor tensor
products, independent unequal Strassen splits, and square feedback. The
fourth pass makes no further change. Final square sizes stay at most 1024.

At size 291, the split is **147+144**. Its seven calls cost

```text
1,063,741 + 3 × 1,036,956 + 3 × 989,982 = 7,144,555.
```

The mixed rectangles exploit inexpensive thin borders around a known
144-square algorithm. This saves **129,279 products (1.78%)** against the
previous count 7,273,834. The gain comes from the improved rectangular
composition; it does not require a new primitive multiplication identity.

The one-product savings have a simpler explanation. At size 137, unequal
Strassen with split 136+1 costs `R(136) + 3 × 136² + 3 × 136`, saving the
single corner product included in ordinary classical border peeling.

## Comparison with size-32 LITA blocking

The current size-32 scheme uses 14,196 products, with recursion exponent
`log_32(14196) = 2.7586393718`. These new points have larger effective
exponents, so they do not improve that asymptotic exponent. Nevertheless,
they improve the counts at their particular finite sizes. The expression
`n^log_32(14196)` is not itself a construction for arbitrary n: pure
size-32 recursion realizes it at powers of 32.

For a concrete uniform blocking comparison, pad to `32 * ceil(n/32)` and
use our current certified square count at `k = ceil(n/32)` for the other
factor. The product count is `14196 * R(k)`:

| Size n | k | R(k) | Uniform size-32 blocking | Current saved count |
|---:|---:|---:|---:|---:|
| 125 | 4 | 48 | 681,408 | 680,814 |
| 291 | 10 | 651 | 9,241,596 | 7,144,555 |
| 1,024 | 32 | 14,196 | 201,526,416 | 198,679,192 |

All 48 new counts beat this uniform comparator. This is not an exhaustive
comparison against every nonuniform blocking strategy. The reported 1.78%
gain at 291 uses the stronger previous envelope (7,273,834), not the weaker
uniform size-32 count. Size 125 belongs to the preceding search; size 1,024
is unchanged and beats current size-32 recursion by about 1.41%.

The full comparison receipt is [saved here](../checks/lita32_block_comparison.json).

## Exactness and reproducibility

The continuation preserves all original leaf meanings and every earlier
positive recipe ID. New IDs are appended, so padding cannot introduce
cycles. The independent audit checks all **5,988,636 nodes**, **11,460,687
child references**, every saved rectangular permutation, and every square
snapshot reference. Hashes bind the previous leaves, remapped history,
comparison envelope and completed snapshot.

The exported certificates retain only their reachable construction graphs.
Every non-LITA primitive is checked by complete tensor expansion over exact
integers after clearing rational denominators. Current LITA leaves use the
pinned published construction theorem, with regenerated small coefficient
objects or the unchanged hash-bound indexed accessor. Large final tensors
are proved compositionally; they have not been expanded in full.

The package uses the same coefficient engine as commit `b647d95`.
[The complete check record](../checks/extended_rectangular.json) includes
primitive proof scope, audit hashes, and each new root object. Counts below
refer to this continuation's baseline; the comparison figure retains the
earliest stated reference when successive searches improve the same size.

## All improved counts

| Square size | Previous envelope | New products | Saving |
|---:|---:|---:|---:|
| 137 | 861,481 | 861,480 | 1 |
| 145 | 1,006,609 | 1,006,608 | 1 |
| 146 | 1,044,715 | 1,040,656 | 4,059 |
| 153 | 1,167,097 | 1,167,096 | 1 |
| 154 | 1,211,453 | 1,204,996 | 6,457 |
| 161 | 1,343,713 | 1,343,712 | 1 |
| 162 | 1,394,733 | 1,385,624 | 9,109 |
| 169 | 1,537,225 | 1,537,224 | 1 |
| 177 | 1,748,401 | 1,748,400 | 1 |
| 178 | 1,814,857 | 1,798,990 | 15,867 |
| 185 | 1,978,009 | 1,978,008 | 1 |
| 186 | 2,052,507 | 2,033,284 | 19,223 |
| 187 | 2,109,033 | 2,104,959 | 4,074 |
| 193 | 2,226,817 | 2,226,816 | 1 |
| 194 | 2,309,010 | 2,286,898 | 22,112 |
| 195 | 2,370,114 | 2,364,873 | 5,241 |
| 201 | 2,495,593 | 2,495,592 | 1 |
| 202 | 2,567,799 | 2,560,780 | 7,019 |
| 209 | 2,785,105 | 2,785,104 | 1 |
| 211 | 2,946,274 | 2,946,273 | 1 |
| 217 | 3,096,121 | 3,096,120 | 1 |
| 218 | 3,205,147 | 3,172,036 | 33,111 |
| 225 | 3,429,409 | 3,429,408 | 1 |
| 233 | 3,785,737 | 3,785,736 | 1 |
| 234 | 3,911,301 | 3,873,196 | 38,105 |
| 241 | 4,165,873 | 4,165,872 | 1 |
| 249 | 4,570,585 | 4,570,584 | 1 |
| 257 | 4,967,233 | 4,967,232 | 1 |
| 258 | 5,116,608 | 5,064,264 | 52,344 |
| 266 | 5,639,088 | 5,598,054 | 41,034 |
| 267 | 5,639,088 | 5,633,862 | 5,226 |
| 273 | 5,861,857 | 5,861,856 | 1 |
| 282 | 6,607,776 | 6,579,849 | 27,927 |
| 289 | 6,857,473 | 6,857,472 | 1 |
| 290 | 7,020,663 | 6,991,672 | 28,991 |
| 291 | 7,273,834 | 7,144,555 | 129,279 |
| 292 | 7,294,224 | 7,176,079 | 118,145 |
| 305 | 7,959,457 | 7,959,456 | 1 |
| 306 | 8,156,538 | 8,109,172 | 47,366 |
| 314 | 8,865,024 | 8,857,050 | 7,974 |
| 321 | 9,173,185 | 9,173,184 | 1 |
| 329 | 10,138,345 | 10,076,472 | 61,873 |
| 330 | 10,163,757 | 10,161,330 | 2,427 |
| 337 | 10,504,033 | 10,504,032 | 1 |
| 338 | 10,760,575 | 10,664,091 | 96,484 |
| 345 | 11,495,276 | 11,485,824 | 9,452 |
| 369 | 13,538,593 | 13,538,592 | 1 |
| 370 | 13,846,063 | 13,757,428 | 88,635 |
