# Cyclo: Lattice Folding That Stops Checking the Accumulator

*A new framework for post-quantum folding schemes that trades unbounded recursion for ~30 KB proofs and a 3.5× faster prover than LatticeFold+.*

**Paper:** https://eprint.iacr.org/2026/359.pdf

## TL;DR

Lattice-based folding schemes — the post-quantum cousins of Nova — are getting practical, but they're held back by one problem: every time you fold a witness into an accumulator, the witness gets longer, and you have to keep proving the accumulator's norm is still bounded. That norm-check is expensive. Cyclo's central observation is that you can skip the accumulator norm-check entirely if you accept a bound on how many times you fold, because without the check the accumulator's norm grows additively rather than multiplicatively. The result is a folding scheme with ~30 KB proofs (down from ~100 KB in LatticeFold+), a prover that's about 3.5× faster on the dominant cost, and a much cleaner derivation when the underlying constraint system lives over a finite field like an R1CS instance over F_q.

## Math you'll need

### Cyclotomic rings and Ajtai commitments

Lattice folding lives in a ring R = Z[X]/⟨Φ_f(X)⟩, where Φ_f is the f-th cyclotomic polynomial. For intuition, take f = 256 and the ring becomes Z[X]/⟨X^128 + 1⟩ — polynomials of degree less than φ = 128 (φ is Euler's totient; gloss it as "the degree of the ring") with integer coefficients, where X^128 collapses to -1. Modular versions R_q = R/qR replace integer coefficients with coefficients in F_q.

The reason to work over R_q at all is the Ajtai commitment: pick a uniformly random matrix A ∈ R_q^{a×m}, and commit to a message w ∈ R_q^m by computing y = A·w mod q. This is binding if and only if w is *short* — meaning ‖w‖_∞ ≤ B for some small B — and the underlying hardness assumption is SIS (Short Integer Solutions): finding any nonzero short z with A·z = 0 mod q is hard. If you've used Pedersen commitments, the structure is the same except the secret has to be small. That smallness requirement is what makes the whole protocol design tricky.

### Folding as a reduction of knowledge

A folding scheme reduces (instance₁, witness₁) and (instance₂, witness₂) to a single (instance_folded, witness_folded). Soundness is captured by saying it's a *reduction of knowledge*: if a prover can convince the verifier of the folded instance, an extractor can rewind it and pull out witnesses for the original two. You fold many times to amortize one expensive SNARK over a long computation, which is the whole point of incrementally verifiable computation (IVC).

In Cyclo (and LatticeFold+) the input relation is what the paper calls the *principal linear relation* Ξ^lin: given an Ajtai matrix A, a list of constraint matrices (M_i), and evaluation points (r_i), find a short w such that A·w = y_top, and ⟨tensor(r_i), M_i·w⟩ = y_i for each i. The tensor here is the standard multilinear-extension trick: tensor(r) is the vector eq(j, r) over the Boolean cube, so ⟨tensor(r), v⟩ = MLE[v](r). Most things you'd want to prove — R1CS, CCS, sum-check residues — reduce to instances of this relation.

### The norm-growth problem

This is the technical headache that every lattice folding paper has to confront. When you fold w_acc and w_input with a short challenge c, the new witness is w_acc + c·w_input. The operator norm of c is some γ ≥ 1, so naively ‖new‖ ≤ ‖w_acc‖ + γ·‖w_input‖. After T folds, ‖w_acc‖ grows by a factor of γ^T — exponential in the recursion depth. Ajtai's commitment stops being binding well before that.

LatticeFold+ deals with this by decomposing the folded witness back into low-norm chunks after each round and running a norm-check on both the input and the accumulator. The decomposition is cheap (only 2 chunks), but the norm-check on the accumulator is what they call the *double commitment*, and it's the dominant cost — you have to compute φ separate Ajtai commitments and then re-commit. Cyclo's contribution is noticing that if you simply don't norm-check the accumulator, the extraction analysis still goes through, because the accumulator's norm grows only by an additive +γb per round (where b is the bound on the input's norm). After T folds, ‖w_acc‖ ≤ ‖initial‖ + T·γb. For T up to 2^20 and reasonable γ, b, this stays well within SIS-binding range.

### Sum-check and multilinear extensions

You probably know sum-check from STARK-adjacent constructions: it's an interactive protocol that reduces "does Σ_{z ∈ {0,1}^ℓ} f(z) = s?" to "does f(u) = v?" for a random u, in ℓ rounds with one univariate polynomial per round. Soundness error is ℓ·deg(f)/|F|. Cyclo uses sum-check twice: once inside the range test (to certify that every coefficient of a witness lies in [-b, b]), and once to unify the evaluation points of the accumulator with those of the input before the actual linear-combination step.

A useful subtlety: sum-check over the ring R_q can be re-cast as φ/e parallel sum-checks over F_{q^e}, where e is the multiplicative order of q mod f. So even though the high-level protocol "lives in R_q", the prover and verifier actually do field arithmetic, which is much cheaper.

## The construction

Cyclo is built from two reductions of knowledge, stacked: an *extension commitment* that decomposes the witness vertically into low-norm chunks, and a *range test* that certifies the chunks are short. Both produce new linear instances that get appended to the input, so by the time you're done folding, you still have a single principal linear relation — just with a few more rows.

### Extension commitment

The extension commitment is what replaces LatticeFold+'s double commitment. The idea is that you've got a witness w of norm at most B, and you want to expose a witness of norm at most b ≪ B so the downstream range test is cheap (range test cost is linear in b).

```rust
// Reduces a witness of norm B to a witness of norm b,
// by base-(2b) decomposition committed under a fresh matrix R.
fn extension_commitment(
    pp: &PublicParams,                  // contains secondary Ajtai matrix R
    instance: &LinearInstance,          // y = A·w mod q,  ‖w‖_∞ ≤ B
    witness: &Witness,                  // w ∈ Rq^m
    challenge_set: &StrongSamplingSet,  // C ⊆ Fq^e, differences invertible
) -> (LinearInstance, Witness) {
    let ell = log_base(2 * pp.b, 2 * pp.B);   // number of digits

    // Prover: decompose w into ℓ chunks of norm < b, concat into v.
    // Invariant: w = Σ_i (2b)^i · v_i,  with ‖v_i‖_∞ < b.
    let v: Vec<Rq> = base_2b_decompose(&witness.w, pp.b, ell);

    // Prover: commit to the *concatenated* v with a fresh matrix.
    // R has its own rank a'; t is small (~13 ring elements in practice).
    let t: Vec<Rq> = pp.R.matmul(&v);

    // Verifier: sample log_2(a) challenges from C, build the tensor.
    let c_hat = challenge_set.sample(log2(pp.a));
    let c = tensor(&c_hat);  // c ∈ Rq^a, multilinear in c_hat

    // The new relation carries one extra row enforcing that v
    // really is a base-(2b) decomposition of *some* w with A·w = y.
    // The row is f^T·v = ⟨c, y⟩, where
    //   f^T = c^T · ((2b)^0, (2b)^1, ..., (2b)^{ell-1}) ⊗ A
    // (a single inner product after the tensor flattens A across digits).
    let new_instance = append_decomposition_row(instance, &t, &c);
    (new_instance, Witness { w: v })   // dimension grew from m to m·ell
}
```

Two things are doing real work here. First, the decomposition is *vertical*: v stays as the witness of a single linear relation, just with m·ℓ rows instead of m. Previous lattice folding schemes decomposed horizontally — splitting one relation into ℓ parallel relations — which inflates the proof. Second, the row f^T·v = ⟨c, y⟩ is what stops the prover from sending a t that's unrelated to the real w. The challenge c is random and tensor-structured, so the Schwartz-Zippel / strong-sampling-set bound gives knowledge error log(a)/|C|.

### Range test

Once the witness has norm at most b, you have to prove it. The trick — older than Cyclo, but used here without batching, which makes it cheap — is that ‖w‖_∞ ≤ b is equivalent to "every coefficient of w lies in [-b, b]", which is equivalent to "the polynomial Π_{i=-b}^{b}(f(X) - i) vanishes on the Boolean cube", where f is the multilinear extension of w's coefficient vector.

```rust
// Proves ‖w‖_∞ ≤ b by sum-checking a product polynomial over F_{q^e}.
fn range_test(
    instance: &LinearInstance,
    witness: &Witness,    // claim: ‖w‖_∞ ≤ b
    b: u64,
) -> (LinearInstance, Witness) {
    let ell = log2(witness.w.len() * pp.phi);

    // f is the multilinear extension of the coefficient vector of w.
    // For z ∈ {0,1}^ell, f(z) is a coefficient; we want each in [-b,b].
    let f = mle_of_coefficients(&witness.w);

    // Verifier samples a "separator" η, prover uses eq(X, η) to absorb
    // the cube into a single sum-check claim.
    let eta = sample_extension(ell);
    let omega = eq_polynomial(&eta);

    // The polynomial whose cube-sum is zero iff every coeff lies in [-b, b].
    // Degree per variable is 2b + 2.
    let f_hat = product_range(&f, b) * &omega;

    // Reduce  Σ_{z ∈ {0,1}^ell} f_hat(z) = 0  to  f_hat(u) = s.
    let (proof, u, s) = sumcheck::prove(&f_hat, ell);

    // Verifier checks Π_{j=-b}^{b}(t - j) · eq(u, η) == s,
    // where t = f(u) — but f(u) is itself an evaluation claim on w,
    // which gets appended to the linear relation rather than checked now.
    let new_instance = append_evaluation_row(instance, &u);
    (new_instance, witness.clone())
}
```

The "appended row" is the load-bearing part. Instead of having the verifier check f(u) = t directly (which would mean reading w), the protocol turns it into another linear constraint on the witness: ⟨tensor(u), cf(w)⟩ = t, expressed in ring-land via a dual-basis trick. The verifier already trusts those constraints, so the range test costs essentially one sum-check transcript and one ring element.

### Folding step

With extension commitment and range test applied to the input witness, the actual fold is almost boring:

```rust
fn fold(
    acc: (LinearInstance, Witness),
    inputs: Vec<(LinearInstance, Witness)>,   // L inputs
    challenges: &ShortChallengeSet,            // D ⊆ Rq, operator norm ≤ γ
) -> (LinearInstance, Witness) {
    // 1. Apply extension commitment + range test to each input
    //    (skipped entirely when the input came from R1CS over F_q with b ≥ k).
    let inputs: Vec<_> = inputs.into_iter()
        .map(|(inst, w)| range_test(extension_commitment(inst, w)))
        .collect();

    // 2. Unify evaluation points across acc and inputs via one batched
    //    sum-check over F_{q^e}. The output: every relation now refers to
    //    the same random points (r̂, b̂).
    let (acc, inputs) = unify_via_sumcheck(acc, inputs);

    // 3. Linear-combine. The accumulator is NOT multiplied by a challenge.
    let s: Vec<Rq> = challenges.sample(inputs.len());   // L short challenges
    let folded_w  = acc.witness + sum(s_j * inputs[j].witness  for j);
    let folded_y  = acc.instance.y + sum(s_j * inputs[j].instance.y for j);

    (LinearInstance { y: folded_y, ..acc.instance }, Witness { w: folded_w })
}
```

The asymmetry in step 3 is what makes the norm grow additively instead of multiplicatively. The accumulator's witness is added in directly, not scaled, so its norm contribution to the next round is exactly ‖w_acc‖ — not γ·‖w_acc‖. The inputs contribute γ·b each, since they're scaled by short challenges and we just range-checked them. After T folds, ‖w_acc^{(T)}‖ ≤ ‖w_acc^{(0)}‖ + T·L·γ·b.

### Folding R1CS over F_q without decomposition

The last move is the cleanest part of the paper. Most R1CS instances live over a prime field F_q, not over a cyclotomic ring. Prior work either packed multiple field elements into one ring element via the NTT (which constrains q to split Φ_f) or did some ad-hoc encoding (which is what Neo does).

Cyclo's framing is that there's a *module homomorphism* θ_k : R_q → F_q given by θ_k(f(X)) = f(k) mod q for some small integer k. It's not a ring homomorphism in general, but it's F_q-linear, which is enough. To encode c ∈ F_q as a ring element, write c in base k and stuff the digits into a polynomial: p_c(X) = c_0 + c_1·X + c_2·X^2 + .... Then θ_k(p_c) = c, and ‖p_c‖_∞ < k by construction. If you pick k ≤ b, the input witness for the folding scheme is *already* short enough, so the extension commitment step can be skipped entirely. The whole expensive part of the protocol disappears when folding R1CS over F_q.

The R1CS → linear relation reduction then runs sum-check over F_{q^e} on the original constraint (M_0·z) ∘ (M_1·z) = M_2·z, using the fact that θ_k(MLE[v](x)) = MLE[θ_k(v)](x) (because θ_k is F_q-linear, and the MLE coefficients are F_q-valued). The prover sends three ring elements d_i that hide the field-side evaluations, and the verifier checks θ_k(d_0)·θ_k(d_1) - θ_k(d_2) matches the sum-check claim. Everything heavy stays in the field; only the commitment lives in R_q.

## Why it's secure

Soundness reduces to SIS in the standard way. If a malicious prover could convince the verifier of the folded instance without holding witnesses for the inputs, the extractor rewinds it on different challenges and uses the difference of two accepting transcripts to either reconstruct the input witnesses (the good case) or produce a nonzero short z with R·z = 0 — a SIS solution under the secondary commitment matrix R. The strong sampling set condition (differences of two challenges are invertible in R_q) is what makes the rewinding clean.

> The genuinely subtle part is the asymmetric treatment of the accumulator during extraction. The folded witness is ŵ = w_acc + Σ c_j · w_input_j. To extract the inputs you forge two transcripts that differ only in c_j, and the difference (ŵ - ŵ') / (c_j - c'_j) yields w_input_j — but this division has *slack*: the extracted witness comes with a short denominator. The paper handles this by defining a "slacked" version of the linear relation that carries the denominator explicitly, and only "discharges" the slack at extraction time. The accumulator gets reconstructed as ŵ - Σ c_j · w_input_j, which inherits an additive norm penalty per fold but no multiplicative blowup. This is the heart of why bounded recursion is the trade-off you're making.

The range test inherits soundness from sum-check over F_{q^e}, with knowledge error roughly ℓ(2b+2)/q^e. For ℓ ≈ 30 (witness with 2^30 coefficients), b = 1, and q^e ≈ 2^100, this is well below 2^-80. The extension commitment's knowledge error is log(a)/|C|; with C a subfield F_{q^2}, |C| ≈ 2^100 and a ≈ 13, this is also negligible. The dominant loss in the union bound is L/|D| (the challenge set for the actual fold), and the paper sizes D so this stays at the target security level.

> The other subtlety worth flagging: the "approximate" strong sampling set Cyclo uses for D — ternary {-1, 0, 1} coefficients, biased to make zero more likely — is *heuristic*. It's not proven that any two such samples differ in an invertible element; the paper gives a numerical estimate that the failure probability is ~2^-94 for their parameters. If you cared about provable security only, you'd use the LS18 construction at a constant-factor cost.

## Performance

The headline numbers, taken from §6.1 and Appendix C of the paper, for 128-bit security and a witness with 2^27 F_q coefficients (a fairly large instance, comparable to LatticeFold+'s benchmark):

| Metric | LatticeFold+ | Cyclo |
|---|---|---|
| Proof size | ~100 KB | ~31.8 KB |
| Dominant commitment cost | 129.4 s | 36.7 s (3.5× faster) |
| Ring degree φ | 64 | 128 |
| Modulus q | ~2^128 | ~2^50 |
| Witness size m | 2^21 | 2^20 |
| Folding depth | unbounded | 2^6 (configurable up to 2^20) |
| Prover memory (typical) | ~16 GB | ~1.56 GB |
| Underlying assumption | SIS | SIS |

The proof-size improvement is the easiest win to understand: LatticeFold+ has to keep two accumulator relations in flight (the input pair) and double-commits each, while Cyclo carries a single accumulator and never decomposes it. The 30 KB number is for a single fold; growing the fold depth from 2^6 to 2^20 only adds about 10 KB, because the dependence on depth is logarithmic.

The 3.5× prover speedup is more interesting because it's measured on just the dominant subroutine — the commitment computation — under matched arithmetic settings (both using AVX-512 IFMA via Intel HEXL, both at q ≈ 2^50). The full pipeline including sum-check should land at a similar ratio, since Cyclo's sum-check work is bounded by ~12 R_q multiplications per witness element (versus ~13 for the commitment), and LatticeFold+'s sum-check load is at least comparable.

The crossover point that matters in practice is the folding-depth cap. If your IVC chain is bounded — sequential rollup batches, ML inference layers, long-running services with periodic snapshots — Cyclo is straightforwardly better. If you need genuinely unbounded recursion, the paper suggests running Cyclo for 2^k rounds and then refreshing the accumulator with one expensive LatticeFold+ step, amortizing the cost. That hybrid still wins overall because the LatticeFold+ step is rare.

## Implementation notes

**No production crate exists.** The reference code at the paper's GitHub repo (https://github.com/osdnk/cyclo) is benchmark-grade, not library-grade. If you're building on this, your foundation is `intel-hexl` (or a Rust wrapper around it) for NTT-domain arithmetic, plus your own implementation of sum-check over F_{q^e} and the extension commitment. The closest reusable crate for the lattice primitives is `lattirust`, but it doesn't have Cyclo's specific decomposition routines.

**Constant-time arithmetic.** The folding prover sees the secret witness throughout. Every operation on w needs to be constant-time with respect to the secret coefficients — which means NTT, base-(2b) decomposition, and the inner products that build f all need scrubbing. The decomposition is the most dangerous: a naive `while x > b { digit = x % b; x /= b; }` leaks the magnitude of each coefficient. Use a fixed-width digit extraction. The sum-check rounds operate on the witness via inner products with public tensors, so those are easier — the public tensor schedules the access pattern.

**Domain separation in the random oracle.** Cyclo's extractor does coordinate-wise forking: it rewinds the prover and reprograms one coordinate of the challenge vector s at a time. This only works if each coordinate is generated by an independent random-oracle query. In practice that means tagging the transcript hash with the coordinate index (`H("s_j" || j || transcript)`) when deriving challenges, not just slicing one long output. Getting this wrong silently breaks soundness — the protocol still runs, but the extraction reduction stops being valid.

**Serialization and field-vs-ring boundaries.** A proof carries elements of three different algebraic structures: R_q (for commitments), R_{q^e} (for unified evaluation claims), and F_{q^e} (for sum-check messages). You'll want a clean encoding that distinguishes them — a tagged-union format rather than raw byte concatenation. The verifier reconstructs challenges via Fiat-Shamir, and any ambiguity in how a ring element vs. a field element gets hashed creates a malleability hole.

**Side channels in the strong sampling set.** Cyclo samples the fold challenge s from a biased ternary distribution. Generating that sample needs to be constant-time with respect to the rejection-sampling rolls — leaking the rejection rate leaks information about the seed. This is solved territory (Kyber and Dilithium both have constant-time samplers you can lift), but it's easy to get wrong.

## Limitations

The bounded fold count is the structural one. The paper presents 2^20 folds as "enough for almost any application," and for many it is — but for unbounded IVC chains or proof-carrying data with adversarial scheduling, you need the hybrid scheme with periodic LatticeFold+ refreshes, and the analysis of that hybrid isn't fully worked out in the paper.

The approximate strong sampling set is heuristic. The paper provides experimental evidence and a numerical estimate of the failure probability, but no proof. If your threat model includes adversaries who might exploit ring-element structure to make non-invertible differences more likely, you'd want the exact strong sampling set from LS18 — at the cost of restricting q and losing some NTT speedup.

Memory is still 1.5 GB for a 2^27-coefficient witness. Cyclo doesn't claim to be streaming-friendly on the prover side; it's just less catastrophic than LatticeFold+. If you're folding on hardware where memory is genuinely scarce, you need something like LatticeFold (which uses the sparser monomial-decomposition technique and has been estimated to fit in ~256 MB with some computational overhead).

Recursive verification — the in-circuit cost of verifying one Cyclo fold inside another Cyclo fold, which is what you'd need for true IVC — is sketched in Remark 1 of the paper but not benchmarked. The claim is that the verifier circuit is smaller than LatticeFold+'s because there's no double-commitment to check, but until somebody actually wires up the recursion, the practical IVC story is unfinished.

Reach for Cyclo when you have a bounded-depth folding chain over R1CS or CCS, want post-quantum security with a small modulus, and care more about proof size and prover wall-clock than about unbounded recursion; reach for LatticeFold+ when you need unbounded folding or your security audit won't accept a heuristic sampling set, and for STARK-style constructions when you want hash-based security without the lattice baggage.
