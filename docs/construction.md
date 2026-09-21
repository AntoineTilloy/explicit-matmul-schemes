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

The input and output dimensions and every coefficient are finite and specified. The constructions are not border-rank approximations and do not rely on numerical near-identities.
