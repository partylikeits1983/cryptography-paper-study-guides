# Soloist: Distributing R1CS SNARKs Without the Proof-Size Penalty

*A distributed SNARK for R1CS that keeps proof size, verification, and per-prover communication constant — even as you scale to 64+ machines.*

**Paper:** https://eprint.iacr.org/2025/557.pdf

## TL;DR

Generating SNARK proofs is the bottleneck for almost every real deployment, and the obvious fix is to throw more machines at it. The trouble is that existing distributed SNARKs for R1CS leak that distribution into the proof: DIZK has linear communication, and Hekaton's proof grows logarithmically in the number of sub-provers while its inter-machine traffic runs in the megabytes per proof. Pianist gets you constant proof size, but it's a Plonk scheme, and many circuits — SHA-256, ECDSA, verifiable FHE — are noticeably cheaper in R1CS than after a Plonk transpile.

Soloist closes that gap. It's a distributed SNARK for R1CS with constant proof size, constant verification, and constant amortized communication per sub-prover, on a universal and updatable setup. The construction combines a new constant-size inner-product PIOP, a "two-dimension sum-check" that lifts it to the distributed setting via a Lagrange trick borrowed from Pianist, a distributed lookup argument that handles online-determined tables, and a batch bivariate KZG that opens many polynomials at many points with O(n) group elements. Concretely, with 64 machines on a general (non-data-parallel) zkRollup circuit, Soloist is roughly 7× faster than Hekaton, has 100× less communication, and uses about 3× less memory than Pianist on the same workload.

## Math you'll need

### R1CS and the PIOP-plus-PCS pipeline

R1CS — Rank-One Constraint System — is the granddaddy SNARK intermediate representation. A satisfying assignment is a witness vector w ∈ F^n together with three public sparse matrices P_a, P_b, P_c ∈ F^{n×n} such that (P_a · w) ∘ (P_b · w) = P_c · w, where ∘ is element-wise multiplication. Every multiplication gate in your circuit becomes one row; addition gates are free. The whole thing is a giant pile of "this product equals that sum."

Modern universal SNARKs follow a two-layer recipe. A polynomial interactive oracle proof (PIOP) reduces the relation to a small set of algebraic claims about polynomials, and a polynomial commitment scheme (PCS) like KZG lets the prover commit to those polynomials and selectively reveal evaluations. The PIOP gives you soundness; the PCS gives you succinctness. Soloist uses bivariate KZG for the PCS — same trapdoor structure as the univariate version, but the SRS is `{[σ^i · τ^j]_1}` and you can commit to f(X, Y) directly. Once you have that machinery, "distributed SNARK" really means "distributed PIOP plus distributed PCS."

### Univariate sum-check

The univariate sum-check is the workhorse subroutine here. Given a multiplicative subgroup H ⊂ F of order m and a polynomial f, you want to prove that Σ_{x ∈ H} f(x) = T. The key algebraic fact: that sum equals m times the constant term of f when f has degree less than m, because every non-constant monomial X^k vanishes when summed over a subgroup of order coprime to k. So the prover commits to f and sends two polynomials g and h such that f(X) = X · g(X) + T/m + Z_H(X) · h(X), where Z_H is the vanishing polynomial of H. The verifier picks a random α (think λ ≈ 128-bit security parameter), queries the three polynomials at α, and checks the identity. Soundness is roughly m/|F|.

A concrete example: to prove that ⟨a, b⟩ = T for two length-m vectors, you build polynomials f_a and f_b from a and b, then express the inner product as a sum-check claim over H. The exact way you build f_a and f_b — evaluation form (interpolation through H) or coefficient form (just dump the entries as coefficients) — matters a lot for prover work. Aurora and Marlin use evaluation form. Soloist uses coefficient form, which saves two FFTs because you skip the interpolation step.

### The Lagrange "Y-dimension" trick

This one is the heart of the distribution story and worth slowing down for. Suppose you have ℓ sub-provers and prover i holds a polynomial f_i(X). You want to operate on the entire collection without sending all ℓ polynomials to a verifier.

Pick a multiplicative subgroup L of order ℓ with generator η, so L = {1, η, η², …, η^{ℓ-1}}. Let L_i(Y) be the i-th Lagrange polynomial of L — that is, L_i(η^{i-1}) = 1 and L_i(η^{j-1}) = 0 for j ≠ i. Now define the bivariate polynomial f(X, Y) = Σ_i f_i(X) · L_i(Y). Evaluating at Y = η^{i-1} gives you exactly f_i(X), because all other Lagrange terms vanish. And critically, Σ_{y ∈ L} f(X, y) = Σ_i f_i(X), which is the sum you wanted in the first place.

So a bivariate polynomial whose Y-domain is the Lagrange basis encodes the entire collection of sub-prover polynomials, and a Y-direction sum over L collapses them. This is the structural move that lets the master prover treat ℓ sub-claims as a single bivariate object, run one sum-check over Y to aggregate them, and produce a proof whose size doesn't depend on ℓ.

### Lookup arguments (specifically logup)

Soloist uses a lookup argument as a subroutine for an R1CS-specific check, so you need the high-level shape. The claim is "every element of the lookup vector a appears somewhere in the table vector b." Haböck's logup observes that this is equivalent to a rational identity: Σ_i 1/(X + a_i) = Σ_j m_j / (X + b_j), where m_j is the multiplicity of b_j in a. Both sides are rational functions in a formal variable X; equality of these functions is a polynomial identity you can check with one random evaluation. Implementations clear denominators and commit to the multiplicities and the "fraction" polynomials, then check identities via sum-check. Soloist extends this to the bivariate, distributed setting and adds machinery for tables whose contents depend on online challenges — which standard lookup arguments don't support, because they assume the table is preprocessed by a trusted indexer.

## The construction

### Step 1: Reduce R1CS to inner products

The classical Aurora trick. An R1CS check `P_v · w = v` (one linear constraint, for v in {a, b, c}) can be batched with a random vector r = (1, r, r², …, r^{n-1}). Define p = r^⊤ P_v; then the matrix-vector equation collapses to a single inner product check: ⟨p, w⟩ = ⟨r, v⟩. The Hadamard check a ∘ b = c becomes a separate inner product check after another linear combination. R1CS reduces to a handful of inner-product claims, and you only need a really good distributed inner-product PIOP.

### Step 2: A leaner inner-product PIOP

The paper's first technical contribution is reducing the number of FFTs and oracles needed for a constant-size inner-product PIOP. Setup: the prover holds vectors f, s ∈ F^m and wants to prove ⟨f, s⟩ = T. Build coefficient-form polynomials f(X) = Σ f_{i-1} · X^{i-1} and s(X) similarly. The key observation is:

> **Lemma 1.** For an order-m multiplicative subgroup H, ⟨f, s⟩ = T iff Σ_{x ∈ H} f(x) · s(x⁻¹) = m · T.

The right-hand side is *almost* a sum-check, except for the s(x⁻¹) evaluated at inverses. The fix is to define s'(X) = X^m · Σ s_{i-1} · X^{-i+1} − s_0 · X^m + s_0, a degree-(m-1) polynomial that agrees with s(X⁻¹) on every x ∈ H (using the identity x^m = 1). After that, you have a clean sum-check, and a trick from Basilisk lets you eliminate the auxiliary low-degree test on g(X) by multiplying through by (X − u) for a fixed u outside H. Net result: 4 oracles, 3 FFTs, 4 field-element proof. Marlin's version is 5 / 5 / 5. Dark's is 5 / 3 / 6.

### Step 3: Distribute via two-dimension sum-check

Now suppose f, s ∈ F^{mℓ} are split across ℓ sub-provers, with P_i holding (f_i, s_i) of length m. You want to prove Σ_i ⟨f_i, s_i⟩ = T while keeping the proof and verifier work constant in ℓ.

The naive thing is to run ℓ parallel sum-checks and then aggregate, which gives a proof linear in ℓ. The fix is to bring in the Lagrange trick: build bivariate polynomials f(X, Y) = Σ_i f_i(X) · L_i(Y) and s'(X, Y) = Σ_i s'_i(X) · L_i(Y), and then the target identity becomes a two-dimension sum-check over H × L.

```rust
// Target identity: prove this evaluates to m·T,
// given oracles to f(X, Y) and s'(X, Y) over H × L.
fn target(f: &Bivariate, s: &Bivariate, H: &Subgroup, L: &Subgroup) -> Field {
    H.iter().map(|x| {
        L.iter().map(|y| f.eval(x, y) * s.eval(x, y)).sum::<Field>()
    }).sum()
}
```

The prover runs an X-direction sum-check first, producing auxiliary polynomials g_1, h_1 ∈ F[X] of size m and reducing the question to checking an identity at X = α for random α. Each sub-prover computes its slice (g_{1,i}, h_{1,i}) locally over its sub-vector. The master prover aggregates by summing — and this is correct because polynomial division by Z_H is unique, so the local quotient pieces literally add up to the global quotient. No coordination is required between sub-provers.

After the X sum-check fixes α, the remaining claim is Σ_{y ∈ L} f(α, y) · s'(α, y) = T_2 for some derived T_2. That's a univariate sum-check over Y of size ℓ, which the master prover can run locally because f(α, Y) is a degree-(ℓ-1) polynomial it can assemble from the {f_i(α)} evaluations the sub-provers send up.

Each sub-prover sends its constant-size evaluation packet — basically {f_i(α), s'_i(α), g_{1,i}(α), h_{1,i}(α)} — plus a couple of polynomial commitments. The master prover does O(ℓ log ℓ) work to assemble the Y-direction proof. Final proof size: constant.

### Step 4: Distributed preprocessing for the public matrices

The verifier needs to be O(1), which means it can't compute the polynomial f_p(α, β) corresponding to p = r^⊤ P_v itself — that's linear in the circuit size. The standard fix from Marlin is to delegate the description of P_v to a trusted indexer via three sparse-encoding polynomials row(X), col(X), val(X) over a subgroup K indexing the non-zero entries.

Marlin's encoding doesn't distribute cleanly. The paper's reformulation uses the coefficient-based encoding: each entry of p contributes a term of the form val(x) · r^{row(x)} · α^{col(x)}, and aggregating over all sub-matrices gives

```
f_p(α, β) = Σ_{x ∈ M, y ∈ L} val(x, y) · A(x, y) · B(x, y) · L(β, y)
```

where A and B are auxiliary bivariate polynomials defined by A(x, y) = r^{row(x, y)} and B(x, y) = α^{col(x, y)}. The catch: these exponential expressions can't be directly recognized as polynomials of bounded degree, so the verifier needs the prover to commit to A and B and then *prove they're consistent with the indexer's row and col polynomials*.

That consistency check is a lookup: prove that every pair (col(x, y), B(x, y)) sits in the table {(j-1, α^{j-1}) : j ∈ [m]}, and similarly for row and A. The table itself depends on the online challenge α, so it can't be preprocessed by the indexer.

> **The subtle part.** Existing lookup arguments assume the table is built up-front by a trusted party. Soloist's tables (α^{j-1} and r^{k-1}) only exist after the verifier sends challenges, so the prover has to prove the table is well-formed as part of the lookup. The paper builds T_col(X) such that T_col(ω^{j-1}) = α^{j-1}, then proves the recurrence T_col(ωX) = α · T_col(X) on H using a single polynomial identity that bundles in the base case T_col(1) = 1 by mixing with the m-th Lagrange polynomial L_{H,m}. The recurrence polynomial p(X) = T_col(ωX) − α · T_col(X) + L_{H,m}(X) · (α^m − 1) has at most m roots, and the construction supplies exactly m of them, so it must be identically zero. No extra quotient commitment is needed.

For the row lookup the table is size mℓ, which would blow up the sub-prover work to O(mℓ log mℓ). The fix is table decomposition: write row(x, y) = row_high(x, y) · √(mℓ) + row_low(x, y) and look up against two tables of size √(mℓ) each. Since √(mℓ) < m in practice, sub-prover work stays at O(m log m).

### Step 5: Batch bivariate KZG with shared Y-coordinate

The PIOP wants to open many bivariate polynomials at many points. Boneh-Drake-Fisch-Gabizon gave a clean batch scheme for univariate KZG with constant proof size; the paper extends it to bivariate.

The single-polynomial, Cartesian-product case: to prove that f(X, Y) takes the right values on A × B for size-t sets A and B, the prover computes a public r(X, Y) of bidegree (t-1, t-1) matching those values, then shows that f(X, Y) − r(X, Y) factors as p(X, Y) · Z_A(X) + q(X, Y) · Z_B(Y) for some quotient polynomials. This is the bivariate analog of "f(α) is the right value iff f(X) − f(α) is divisible by (X − α)."

The multi-polynomial case is the hard one, because each polynomial has its own A_k and B_k, and you can't naively common-denominator the quotients across different polynomials. The paper introduces an auxiliary polynomial s_k(X, Y) that agrees with f_k(X, Y) on A_k × everywhere and with r_k(X, Y) on everywhere × B_k. This splits the seam between the two dimensions and lets you batch p_k's and q_k's separately with a single random linear combination.

For Soloist specifically, all evaluation points share the same Y-coordinate β, because the verifier picks a single β for the Y-direction sum-check. That collapses the bivariate batch to something much closer to univariate batching, and the proof size drops to 4 |G| + (n+1) |F|. For n = 10 polynomials on BN254 that's 4 · 32 + 11 · 32 = 480 bytes for that piece of the proof, versus roughly 2 KB for the general bivariate batch.

### Step 6: Distribute the PCS

The final piece. Pianist showed how to distribute single-point bivariate KZG by choosing the quotient polynomials so that one of them depends only on Y (and so can be computed by the master prover alone) and the other splits cleanly over the Lagrange basis. Soloist generalizes this to the batch case using the same s_k machinery:

```rust
// Sketch of distributed batch opening for ℓ sub-provers.
// Each P_i has SRS slice {[σ^j · L_i(τ)]_1} and holds {f_{k,i}(X)}_k.
fn sub_prover_open(
    f_ki: &[Univariate],     // P_i's slice of each f_k(X, Y), one per polynomial
    a_sets: &[Subset],       // X-coordinates being opened, one set per polynomial
    gamma: Field,            // batch challenge from the verifier
    srs_slice: &SrsSlice,    // L_i(τ)-scaled monomials, owned by P_i
) -> SubProof {
    // Build s_{k,i}(X) such that s_{k,i}(α) = f_{k,i}(α) for α in A_k.
    // Invariant: s_{k,i} has X-degree < t even though f_{k,i} has degree < m.
    let s_ki: Vec<Univariate> = f_ki.iter().zip(a_sets)
        .map(|(f, a)| interpolate_through(f, a))
        .collect();

    // Local quotient piece. The global Z_A(X) divides
    //   Σ γ^{k-1} · (f_{k,i}(X) − s_{k,i}(X)) · Z_{A \ A_k}(X),
    // which the master prover will assemble by summing across i.
    let p_i = compute_local_quotient(f_ki, &s_ki, a_sets, gamma);

    // Commitments are pre-scaled by L_i(τ) via the SRS slice, so the
    // master prover gets the global commitment by group addition alone.
    SubProof {
        cm_s: s_ki.iter().map(|s| commit_scaled(s, srs_slice)).collect(),
        cm_p: commit_scaled(&p_i, srs_slice),
    }
}
```

The master prover sums commitments to assemble the global U_1, computes the Y-only quotient q(X, Y) locally — its size depends only on ℓ and t, not on m — and emits the final proof. The per-prover SRS is only O(m) elements; the full O(mℓ) SRS lives in aggregate across the cluster. That's a meaningful operational property. A 64-machine cluster proving a circuit of size 2^25 doesn't need any single machine to hold all 2^25 SRS elements.

## Why it's secure

Soundness follows the standard PIOP+PCS recipe: if the PIOP has soundness ε_p and the PCS has knowledge soundness, the compiled argument has knowledge soundness with error ε_p plus the PCS extractor's failure probability.

The PIOP's soundness has three layers. The reduction from R1CS to inner products is Freivalds's algorithm, with error mℓ/|F|. The two-dimension sum-check inherits soundness from the standard univariate sum-check on each axis, with Schwartz-Zippel bounds at each random challenge — total error roughly (2m + 3ℓ)/|F| for the inner-product piece. The lookup argument adds logup's soundness (essentially a polynomial identity check) plus the online table-validity check, whose correctness rests on the m-root argument for the recurrence polynomial mentioned earlier.

The KZG-based PCS gets knowledge soundness in the algebraic group model. If you could open a commitment to two different polynomials, you would have produced a polynomial g(X, Y) such that g(σ, τ) = 0 even though g is nonzero — a discrete-log relation in disguise, which is the standard AGM hardness assumption.

> **The subtle part of the security proof.** The lookup argument with online-determined tables doesn't fit the standard logup soundness analysis, because the table itself is an output of the prover rather than a fixed indexer artifact. The paper has to argue that any prover who satisfies the online table-validity check on a *wrong* table T'_col either fails the recurrence check at random δ (with probability 1 − m/|F|) or produces a polynomial with too many roots. The whole induction "T_col(ω^{j-1}) = α^{j-1} for all j" is encoded in a single algebraic identity rather than checked point-by-point, and that's what makes it cheap enough to embed inside the main proof.

For the distributed PCS, AGM extraction still works because the master prover's final commitment is a linear combination of sub-prover commitments under the Lagrange basis, and the algebraic adversary's representation of those sub-commitments lets the extractor recover each sub-polynomial individually. That's also what makes accountability cheap: each sub-prover's contribution is, with minor modifications, an independent KZG proof that the master can verify before aggregating. A malicious sub-prover is caught locally rather than corrupting the global proof.

## Performance

The numbers below are extracted from the paper's experiments on BN254 (BLS12-381 for the Marlin baseline), with each machine running 16 cores and 256 GB of RAM. Here n is the R1CS constraint count, ℓ is the number of sub-provers, and m ≈ n / ℓ is the per-prover slice size.

| Scheme | IR | Proof size | Verifier | Per-prover work | Amortized comm. | Per-prover SRS |
|---|---|---|---|---|---|---|
| Groth16 (non-distributed) | R1CS | ~192 B | ~1.5 ms | O(n log n) | — | O(n) |
| Marlin (non-distributed) | R1CS | ~600 B | ~6.5 ms | O(n log n) | — | O(n) |
| Pianist | Plonk | ~2.2 KB | ~2.8 ms | O(m log m) | ~2 KB | O(m) |
| HyperPianist | Plonk (multilinear) | O(log n) ≈ 15 KB | ~4 ms | O(m log m) | O(log n) | O(m) |
| Hekaton | R1CS | O(log ℓ), ~15–25 KB | ~83 ms | O(m log m) | ~10 MB | O(m) |
| DIZK | R1CS | ~192 B | ~1.5 ms | O(m log² m) | O(n) | O(n) |
| **Soloist** | **R1CS** | **~10 KB** | **~3 ms** | **O(m log m)** | **~7 KB** | **O(m)** |

Soloist's prover scales linearly with the number of machines on both data-parallel and general circuits, and the latter is the key result. Hekaton's prover time actually gets *worse* on general circuits as you add provers, because more sub-provers means more shared wires across sub-circuits and the aggregation overhead grows. At ℓ = 64 on a general zkRollup workload, Soloist is about 7× faster than Hekaton and roughly comparable on nearly-data-parallel Merkle-tree circuits. Against Pianist on the same R1CS workload, Soloist is 1.8× faster in proving, 2.8× faster in preprocessing, and uses about 3× less memory — much of that comes from skipping the R1CS-to-Plonk transpile, which inflates the constraint count roughly 10× on this circuit.

The crossover point against Hekaton sits around ℓ = 8 to 16. Below that, Hekaton wins on raw prover time because its sub-proofs use Groth16-style machinery and its aggregation overhead is small. Above that, the aggregation tree starts to dominate and Soloist pulls ahead. The crossover against HyperPianist runs the other direction. HyperPianist's multilinear PIOP gives it about a 4× prover speedup, but its proof size, communication, and verifier time all grow logarithmically. For applications where the verifier is on-chain — zkRollups, bridges, any setting where someone pays gas per proof — Soloist's flat numbers win even though raw proving is slower.

The proof size is the honest weakness. 10 KB is fine for off-chain settings and fine for an L2 rollup that posts proofs to L1 every few minutes, but if you're trying to fit a proof in a 1 KB envelope (some bridge protocols, some Bitcoin-anchored settings) it's too big. The paper points out that proof size is amortizable: Soloist can prove much larger batches in a single proof than Pianist — 320 transactions per prover versus 16 — because it doesn't hit the BN254 multiplicative-subgroup size limit. The per-transaction proof bytes can therefore be competitive even though the per-proof bytes are not.

## Implementation notes

**Existing crates.** The authors implemented Soloist in roughly 10,000 lines of Rust on top of `arkworks` and say they plan to open-source. As of the paper, nothing is public. If you wanted to build this today, you'd start from `ark-poly-commit` for KZG, `ark-bn254` or `ark-bls12-381` for the curve, and `ark-poly` for the polynomial machinery. The distributed PIOP needs a coordinator that owns the network — `tokio` plus `tonic` for gRPC between sub-provers and the master is the obvious choice. Pianist's released code demonstrates the Lagrange-SRS trick for the PCS layer; that's the closest existing reference.

**SRS structure matters a lot.** The paper uses a Lagrange-basis SRS where sub-prover i gets `{[σ^j · L_i(τ)]_1}` for j ∈ [m], not the monomial basis. This is what enables each sub-prover to produce a commitment whose group sum across sub-provers yields the correct global commitment. Generating this SRS is non-trivial: you need a trusted ceremony for the standard monomial SRS and then a one-time conversion. The conversion is itself an O(mℓ log(mℓ)) FFT in the group, which is annoying but happens once per (m, ℓ) configuration.

**Constant-time concerns.** R1CS provers don't typically run on user secrets in a side-channel-sensitive way — the witness is usually known to the prover and the proof is published. The constant-time discipline is mostly about not leaking which sub-circuit is harder to prove, which matters only if your distribution algorithm reveals workload skew. The paper's "making-even" algorithm is a public preprocessing step over the matrix structure and doesn't depend on the witness, which is the right design.

**Domain separation.** The Fiat-Shamir transform has to absorb commitments to many polynomials across multiple rounds. The natural mistake is to forget to absorb the *bivariate* structure — if your transcript only hashes the commitment to f(X, Y) without also binding the (m, ℓ) parameters and the Lagrange basis choice, an adversary can replay a transcript across different parameter settings and the soundness argument breaks. Bind everything: (m, ℓ, n_polys, point structure, indexer commitments) before the first challenge.

**Serialization.** Bivariate KZG commitments are single group elements, so on-chain verification cost is dominated by the four pairings in the batch opening check. Make sure your pairing precompile handles G2 element compression — BN254's G2 elements are 64 bytes compressed and 128 uncompressed. For the (n+1) field elements in the shared-Y batch scheme (32 bytes each on BN254), you'll want a stable little-endian wire format with explicit length prefixes, because the proof is no longer fixed-shape in n.

## Limitations

The concretely large proof size is the headline limitation. 10 KB is two orders of magnitude bigger than Groth16's 192 bytes and a few times bigger than Pianist's 2 KB. The paper traces this to the number of polynomial oracles in the preprocessing PIOP — the lookup machinery and the bivariate encoding of R1CS matrices both contribute. A future construction with a more compact preprocessing PIOP, perhaps one not based on online lookup, would close this gap.

The trade-off against multilinear schemes is real. HyperPianist and Cirrus have roughly 4× faster provers on the same hardware because multilinear PIOPs avoid FFTs over large smooth subgroups. If your application doesn't care about constant proof size or constant verifier work — for instance, if you're proving recursively and the inner verifier is itself implemented in a SNARK where O(log n) verification is fine — those schemes are probably the better fit. Soloist's argument is exactly that constant verifier work matters when the verifier sits on a chain that charges per gas, and that's the regime where it wins.

The "making-even" algorithm is heuristic. It produces R1CS sub-matrices where the maximum non-zero count is close to the average, which is what you need for full linear speedup, but it isn't provably optimal — multiway number partitioning is NP-hard. For pathological circuits with extremely skewed sparsity, the speedup will be sub-linear, and the paper offers no a-priori bound on how bad it can get.

Finally, the setup is universal and updatable (the same posture as Marlin and Plonk), but it isn't transparent. If you can't run a ceremony or trust an existing one, you need a transparent SNARK like Brakedown, Ligero, or one of the FRI-based schemes, and you'll pay for it in either proof size or verifier time.

Reach for Soloist when you have R1CS circuits, a cluster of provers, and an on-chain verifier you want to keep cheap; reach for Hekaton when your circuits are genuinely data-parallel with very few shared wires and Groth16-sized proofs matter; reach for HyperPianist or Cirrus when prover speed is everything and you can absorb logarithmic verifier work downstream.
