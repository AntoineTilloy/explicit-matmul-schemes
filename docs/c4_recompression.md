# Rational C4 recompression at size 1,024

The size-1,024 square now uses **198,679,192** products, **4,096** fewer than
the preceding rank 198,683,288. The support is unchanged: 28 balanced power-9
CW rows in seven free C4 orbits, with eight singleton, 24 pair and four triple
sunflower kernels. The improvement comes from the choice of rational pair
relations and retained patterns.

Zero-factor removal leaves 1,452 patterns. A pair-(0,2) compression deletes
304 complete group-invariant patterns, leaving 1,148. Its only modified
factor is packet mode 1, which is the physical output C under the relation
permutation `[1,2,0]`. Both physical input factors retain their original
coordinate equivariance.

The independently counted retained-pattern polynomials are

```text
4 R_even(q) = q^9 + 17q^8 + 122q^7 + 392q^6 + 493q^5 + 127q^4 + 2q^3 + 2q^2
4 R_odd(q)  = q^9 + 17q^8 + 122q^7 + 392q^6 + 493q^5 + 127q^4 - 2q^3 - 2q^2
```

The old q⁴ coefficient was 131. At q=8 the new odd unit uses 198,679,264
products. The same signed Strassen embedding as before still has 72 zero
input factors, giving the stated square count. Every deletion ordinal is
rebound to the new unit and replayed from its coefficient tables. The total
saving against the original even-character rank 198,683,936 is 4,744.

## Exact averaged relations

Write the packet as `T = sum_h P_h ⊗ F_h`, where `P_h` is the product of the
two unchanged factors. The kept set is group invariant, `g P_h = P_(gh)`,
and T is coordinate invariant. Replacing the modified factor by

```text
Fbar_h = (1/4) sum_g g^-1 F_(gh)
```

preserves T, because the resulting decomposition equals
`(1/4) sum_g g^-1 T = T`. An arbitrary Gram basis may therefore be used and
its free-factor update averaged afterwards.

In the centered-label representation, the group sends the standard component
I, row zero mask Z and pattern H to `gI, gZ, gH`. For a relation
`P_(I,J) = sum_B C_I(B,J) P_(I,B)`, the explicitly stored averaged coefficients
are

```text
Cbar_I(B,J) = (1/4) sum_g C_(gI)(gB,gJ).
```

Absent coefficients are zero. The basis is the union of the transported
retained bases. It may be linearly dependent; only the displayed exact
relations are needed. The certificate stores Cbar explicitly, so the public
reader uses the existing relation format and needs no implicit averaging step.
Coordinate and role transport maps the discovery rows onto the existing
support, with every coefficient row and column reordered exactly.

## Proof scope and replay

Separate exact rational witnesses were produced at q=2 and q=8. Each source
replay checks 3,213,436 original and 3,213,436 averaged Gram equations and the
equality between the explicit update and direct free-factor averaging. The
published q=8 object carries `certified_q: 8`; the q=2 companion carries
`certified_q: 2`. Neither is licensed by the count polynomial for other q.
Larger catalogue entries keep their historical relation objects and counts.

An independent size-16 check verifies every 16⁶ coefficient of the q=2 odd
square at rank 15,838 using a rigorous rational CRT bound. This check uses
the unpermuted, ungauged Strassen embedding. Exact coordinate transport,
the 64-coefficient rational outer MM2 identity and the separate physical
zero-factor replay certify the gauged embedding. Its six q=2 cancellations
give rank 15,832; this companion is not a size-16 catalogue record.

`cwschemes verify 1024 --exact` reconstructs and checks the published q=8
relations. The default verifier also rebuilds all 72 declared zero factors
exactly. Integrity, counts, complete relation replay, small full-tensor scope
and installed-wheel access are recorded separately in
[checks/small1024_recompression.json](../checks/small1024_recompression.json)
and its companion check files. No characteristic-two inference is used.
