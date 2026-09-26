# Finite exact constructions

All schemes compute square matrix multiplication over Q. The multiplication count is the number of rank-one bilinear terms, excluding additions and fixed-coefficient scaling.

The little-CW local identity has three matrix-vector legs. An exact rational completion writes it as q+2 products plus an all-positive error. In a tensor power, each retained sector has a specified zero mask. Auditing all cross-sector combinations ensures that the selected tensor is a direct sum of independent matrix products. The grading eliminates the unwanted all-positive error. Surviving pure subsets determine the exact packet rank; their counts are recomputed from the support witness.

A square outer bilinear scheme asks for a fixed number of inner square products. Divide these calls among the independent products in the CW packet. Further packets, ordinary recursive schemes or a partially filled final packet supply the remainder. A partially filled packet means its unused inputs are zero. Tensor product and padding preserve exactness over Q.

Other entries retain the fixed sectors of a finite involution. The six free sector pairs yield six full matrix products, while the fixed sectors split into products in rational eigenspaces. An unequal Strassen composition places those products in the seventh outer call and completes its remaining rectangular pieces with explicit ordinary recipes. Rectangular subproblems are an implementation ingredient; all delivered target problems are square.

Rational compression removes a term when its pair factor lies in the span of retained pair factors. If `u_j ⊗ v_j = Σ_b C[b,j] u_b ⊗ v_b`, transfer its output factor into those retained terms. The reduction preserves the tensor exactly. Symmetry permits checking the relations on finite multiplicity blocks rather than expanding every labelled column. Refined schemes split selected label orbits before repeating these relations. Successive schemes change factor pairs and reconstruct the corresponding updated Grams.

The source reader implements the resulting coefficient maps, not the search that selected a support or deletion set. Proof files store the specific support, rational coefficients and composition recipes. The exact checker reconstructs the relevant relations without choosing new supports or solving an optimization problem.

The largest new support has power 21 and 5,032 independent size-q^7 products. Its exact packet rank is

```
P(q) = (q+1)^21 + 21(q+1)^20 + 210(q+1)^19 + 1330(q+1)^18
       + 5985(q+1)^17 + 15583(q+1)^16 + 2500(q+1)^15 + 9(q+1)^14.
```

At q=16, the original published size-22 outer scheme has 5,155 calls. One packet and 123 recursive size-16^7 calls give

```
P(16) = 229066385161157064226434567
R(16^7) = 69980774018591548106130
R(22*16^7) <= P(16) + 123*R(16^7)
           = 237674020365443824643488557.
```

## Z4-quotient packets

The entries between sizes 1,024 and 65,536 use a power-9 support whose 28 rows are closed under a coordinate permutation g of order 4, with cycles (0 1 2 3)(4 5 6 7) and 8 fixed. The rows form 7 free orbits. The support is not an induced matching: 36 sunflower kernels K (profile 8, 24, 4 by |K|) receive Möbius-inverted corrections. The packet tensor is

```
T = G_q^{⊗9}|supp + Σ_K β_K (−D)^{⊗K} ⊗ G_q^{⊗(9−|K|)}|supp = 28 copies of ⟨q^3⟩.
```

Write it in the joint triangular representation of the mode-permuted rows. Packet modes (0, 1, 2) read the original (B, C, A). Columns are patterns H ∈ {pure, ordinary, special}^9 with labels on the ordinary coordinates. Two relation steps follow, both fixed independently of q.

1. Remove the 1,816 patterns whose factor vanishes on the support. 1,452 patterns remain.
2. Apply a pair (0, 2) Gram elimination. Its free mode is the output C. Deleted patterns form whole Z4-orbits (300 patterns), so 1,152 patterns remain.

Because the input factors are never modified, they stay Z4-equivariant; the kept set is Z4-stable. On Z4-invariant inputs, the columns of one orbit therefore compute the same product. Sum their output factors and keep one multiplication per column orbit:

```
T_quot = Σ_O u_O(X) v_O(Y) Σ_{t∈O} w_t |_rep.
```

This computes the 7 orbit products independently from 7 independent pairs of inputs. It is a "unit" for 7 disjoint ⟨q³⟩ products. By Burnside's lemma, its rank is

```
R7(q) = (q^9 + 17q^8 + 122q^7 + 392q^6 + 493q^5 + 131q^4 + 2q^3 + 2q^2)/4.
```

Strassen's outer scheme gives R(2q³) ≤ R7(q): 198,683,936 at q = 8 and 74,352,484,826 at q = 17 (n = 9,826, exponent 2.723013). An outer scheme of rank R_k with 7j of its products grouped into units gives R(kq³) ≤ j·R7(q) + (R_k − 7j)·R(q³). The catalogue uses this at 10,976, 16,384 (Strassen², j = 7), 32,768 (⟨8;336⟩, j = 48) and 65,536 (⟨16;2236⟩, j = 319, three products by the size-4,096 entry).

The input and output dimensions and every coefficient are finite and specified. The constructions are not border-rank approximations and do not rely on numerical near-identities.
