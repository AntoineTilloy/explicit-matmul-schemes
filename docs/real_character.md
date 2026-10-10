# Rational real-character restriction of the Z4 CW packet

This row identification changes the two input embeddings and term unranking
of the existing certified packet, retaining its relation coefficients and
Schur-updated output factors.

Let `g=(0123)(4567)` and let the real character satisfy `chi(g)=-1`.
For a packet row `rho=g^j rep_h`, embed the two input matrix variables with
the scalar `(-1)^j` and the existing transported inner-index map. Keep the
output restriction at the seven representative rows. The independent packet
identity has only matching row sectors. Its surviving row at each output
representative has both input signs equal to one, so it becomes exactly the
seven inner products in the existing outer Strassen decomposition.

The input factors of the retained packet columns are coordinate equivariant.
The initial triangular transformation is coordinate natural; zero removal
changes no factor; subsequent Gram relation steps change only original mode C.
Hence the signed input factor of the translated column satisfies
`U(g*t)=chi(g) U(t)` and likewise `V`. Their product is invariant. Therefore
one multiplication per column orbit suffices, with the existing output sum.
No equivariance of the Schur-updated output factors is required.

A column orbit survives precisely when its stabilizer lies in `ker(chi)`.
For this Z4 action all length-two column orbits have stabilizer `{1,g^2}` and
survive. A length-one column is fixed by g, so both signed input factors equal
their own negatives and are zero over Q. All length-four orbits survive.
The fixed-column polynomial of the final retained pattern set is `q^3+q^2`.
Thus

`4 R_odd(q) = q^9 + 17q^8 + 122q^7 + 392q^6 + 493q^5 + 131q^4 - 2q^3 - 2q^2`.

At q=8, n=1024, this is rank 198,683,360, 576 fewer than the previous bound.
At q=2 the rank is 15,854. The q2 orbit enumeration has 12 length-one orbits,
34 length-two orbits, and 15,820 length-four orbits. Unranking excludes all
and only the twelve forbidden orbits, including fixed points of every subgroup.

`cwschemes._engine.z4_checks.full_exact_check` checks every 16^6 square coefficient over Q. It calculates
rational upper bounds B_h on the absolute factor entries from the actual
coefficient tables and the exact q2 centered-factor magnitude 1/2. D_h are
integer common factor denominators. The CRT modulus exceeds
`2 D_0 D_1 D_2 (R B_0 B_1 B_2 + 1)`, and also each input-equivariance bound.
Every modular contraction uses integer-valued floating-point dot products
strictly below 2^53, so arithmetic is exact. It checks all twelve vanishing
fixed factors, all thirty-four length-two input orbit sign identities, and
eight length-four representatives after the existing output Schur updates.

The result is rational. Its denominators and the implication `U=-U => U=0`
do not provide a characteristic-two scheme.

## Further square-only cancellation

At q=8, permuting the seven Strassen calls by `[1,0,5,2,3,4,6]` and applying input gauges
`[1,1,-1,-1,-1,-1,-1]` and `[1,-1,-1,-1,-1,-1,-1]` removes another 72 products.
The output gauges are the reciprocal input products, so the outer tensor remains exactly MM2.
The square rank is **198,683,288**, a total saving of **648** against the earlier even-character entry.

The cancellation occurs in three actual label classes. Pattern 1093 with all seven ordinary labels equal has a zero B factor (8 columns). Pattern 10087 with all seven ordinary labels equal has a zero A factor (8 columns); its equality partition `0000001`, with the label at physical coordinate 8 different from the other six ordinary labels, also has a zero A factor (56 columns). The bound uses all numerical columns in these classes, with original unit ordinals frozen in the certificate. The underlying direct-sum unit still has rank 198,683,360: these factors vanish only after the square embedding.

`GaugedZ4Square.audit_deletions()` binds the proof to the exact outer permutation, gauges, original indices, patterns and numerical labels. It reconstructs the physical factor from the public packet tables in four outer blocks and six physical digit slots, and checks that every rational polynomial coefficient is zero. Simultaneous label relabeling reduces the 72 checks to three formal polynomials. This formal zero test is valid even when centered label functions have additional dependencies.

The size-16 q=2 instance removes six further columns, leaving rank 15,848. An independent public check compares scalar coefficient queries to the vectorized factors and checks all 16⁶ tensor coefficients with a rigorous CRT bound. The wrong-gauge and wrong-deletion regression controls retain a valid outer MM2 identity but fail the actual zero-factor check. See [the public replay](../checks/small1024_gauge.json). No characteristic-two argument is used in this cancellation.
