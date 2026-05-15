# Rackle: Three-Round Threshold Schnorr from the Discrete Log Assumption Alone

*The first adaptively-secure threshold Schnorr signature scheme to hit three rounds without strengthening any cryptographic assumption.*

**Paper:** https://eprint.iacr.org/2025/1941.pdf

## TL;DR

Threshold Schnorr lets N signers jointly produce a single, standard-format Schnorr signature when any T of them cooperate. Until now you had to pick two of three: low round count, adaptive corruption resistance, or DL-only security. Sparkle and FROST achieve few rounds but lean on AOMDL or AGM; Gargos and Glacius get adaptive security in three rounds but require DDH; Crackle hits DL but needs five rounds. Rackle, in this Eurocrypt 2026 paper by Niot, Reichle, and Takemure, closes the gap: three rounds (one preprocessable), full adaptive security, no secure erasures, no authenticated channels, and security under DL in the random oracle model. The mechanism is a new "product equivocation" property for the BCJ commitment scheme paired with a simulation-extractable randomized Fischlin NIZK, glued together with a masking trick that forces consistency among signers without an explicit confirmation round.

## Math you'll need

### Adaptive vs static corruption

In the static model, the adversary names up to T−1 corrupted signers before the protocol starts. That's a polite assumption: in reality, an attacker compromises machines opportunistically, often after watching signatures fly past. Adaptive security says the adversary can pick whom to corrupt at any point, including mid-protocol, having already seen partial transcripts. The paper goes further with what it calls "no secure erasures": when signer *i* is corrupted, the adversary recovers every random coin *ρ* that *i* ever sampled, not just its current state. This is the strongest reasonable corruption model — anything weaker assumes the signer's process had time to scrub memory, which in practice it usually didn't.

This matters for the reduction. The simulator has to be ready, at any moment, to hand the adversary a plausible secret key share and a complete history of randomness for any honest signer. If the simulator committed itself to specific values during signing, late corruption can force it into a contradiction.

### Equivocal commitments

A commitment hides a value and binds you to it: you publish a commitment *c* to message *m*, and later you can open it to reveal *m* but not anything else. An *equivocal* commitment has a back door: with a trapdoor, you can produce a commitment that opens to anything you want. Pedersen commitments are equivocal over field elements: *c = g^m · h^r*, and if you know logₕ *g* you can open *c* to any *m* you like by adjusting *r*.

For threshold Schnorr you need to commit to group elements (the nonces *R = G^r*), and the natural choices have a tension: hash-based commitments (just publish H(R)) are equivocal in the ROM but require revealing R to open, which contradicts adaptive security. Structured commitments can be set up either equivocable *or* extractable but not both simultaneously without help.

The paper uses BCJ commitments [Bagherzandi-Cheon-Jarecki 2008], which commit to a group element *M* under key (G, H, Y, Z) as the pair (G^s · H^t, Y^s · Z^t · M) for random (s, t). The notation is (C₀, C₁) ← BCJ.Commit(ck, M; (s,t)). BCJ is binding under DL, perfectly hiding, and crucially *homomorphic*: the product of commitments to M₁ and M₂ is a commitment to M₁·M₂ under randomness s₁+s₂.

### Straight-line simulation-extractability

A NIZK is *extractable* if there's an extractor that, given an accepting proof, can recover the witness. Most ROM-based extractors work by rewinding the prover, which is fine for a security proof but useless when you need to extract on-the-fly during a live protocol with many concurrent sessions. *Straight-line extraction* (SLE) means the extractor reads proofs as they arrive, without rewinding.

*Simulation-extractability* (SIM-EXT) is the stronger property the paper needs: even while the extractor is simulating proofs for honest parties (using a trapdoor), it can still extract witnesses from any proof an adversary produces. The paper proves the randomized Fischlin transform [Kondi-shelat 2022] satisfies SIM-EXT with two extra refinements: proofs carry a *label* (so a simulated proof for label 0 doesn't poison extraction for label 1), and simulated proofs are *explainable* — given a witness, the simulator can produce well-distributed randomness *ρ* such that running the honest prover on *ρ* would have produced the simulated proof. Explainability is what lets the simulator survive adaptive corruption without erasing the original coins.

### Zero shares (the masking trick)

If T signers each compute a value vᵢ and you want only their sum Σvᵢ to be visible, the textbook approach is additive secret sharing in reverse: each signer broadcasts vᵢ + Δᵢ where Σᵢ Δᵢ = 0 and the Δᵢ look individually random. Pairwise PRG seeds make this cheap. Each pair (i, j) of signers shares seeds seedᵢⱼ and seedⱼᵢ, and signer *i* sets

Δᵢ = Σⱼ (H(seedᵢⱼ, input) − H(seedⱼᵢ, input))

The terms cancel pairwise in the sum, so Σᵢ Δᵢ = 0. The paper notates this as Δᵢ = ZeroShare(seedsᵢ[SS], input).

The critical hygiene rule is that *input* must be fresh per session. If you reuse input, the same Δᵢ shows up twice, and an adversary subtracting two masked values learns vᵢ¹ − vᵢ². Rackle ensures freshness by including the full set of commitments in *input*. This is also how the protocol implicitly enforces that all honest signers agree on what they're signing: if any honest signer's view of "the set of commitments" differs by even one bit, their Δᵢ values stop summing to zero, the aggregate breaks, and no signature emerges.

## The construction

The verification key is a standard Schnorr public key *X = G^x*, with *x* distributed via T-of-N Shamir secret sharing so that signer *i* holds *xᵢ* with *x = Σ Lᵢ · xᵢ* over any T-subset. Each signer also holds pairwise seeds with every other party, and a BCJ commitment key *ck* generated during setup.

### Round 1 (offline, preprocessable)

Signer *i* samples a Schnorr nonce, commits to it under BCJ, and attaches a NIZK proving knowledge of the opening:

```rust
struct Round1Message {
    commitment: BcjCommitment,  // (C0, C1) ∈ G²
    proof: FischlinProof,       // proves knowledge of (r_i, s_i)
}

fn sign_round1(sk_i: &SecretShare, ck: &CommitmentKey) -> (Round1Message, State1) {
    let r_i = Scalar::random();
    let R_i = G * r_i;                        // the actual nonce
    let s_i = (Scalar::random(), Scalar::random());
    let C_i = bcj_commit(ck, R_i, s_i);       // hides R_i

    // Prove knowledge of (r_i, s_i) such that C_i opens to G^r_i.
    // Label is the signer's own index — this matters for the security proof.
    let stmt = Statement { ck, commitment: C_i };
    let witness = Witness { r: r_i, s: s_i };
    let proof = fischlin_prove(&stmt, &witness, /* label = */ i);

    let msg = Round1Message { commitment: C_i, proof };
    let state = State1 { C_i, r_i, s_i };
    (msg, state)
}
```

The message and signer set are not yet known. This entire round can be batched and precomputed.

### Round 2 (online): masked opening

After receiving everyone's commitments, signer *i* verifies the NIZKs, fixes its view of the session context, and broadcasts a masked opening:

```rust
fn sign_round2(
    sk_i: &SecretShare,
    state1: State1,
    SS: &SignerSet,
    msg: &[u8],
    received: &[Round1Message],
) -> (Round2Message, State2) {
    // Verify everyone's proof of knowledge. Any failure aborts the session.
    for (j, r1) in SS.iter().zip(received) {
        assert!(fischlin_verify(&Statement { ck, commitment: r1.commitment },
                                &r1.proof, /* label = */ j));
    }

    // Context binds the signer set, the message, and *every* commitment.
    // This is the input that has to be fresh for zero-shares to mask safely.
    let ctnt_tilde = encode(b"open", SS, msg, received.commitments());

    let delta_tilde = zero_share(&sk_i.seeds, SS, &ctnt_tilde);  // ∈ Z_p²

    let s_tilde_i = state1.s_i + delta_tilde;  // masked opening

    let state = State2 { SS, msg, commitments: received.commitments(), 
                         r_i: state1.r_i, s_tilde_i };
    (Round2Message { s_tilde: s_tilde_i }, state)
}
```

Two things are happening here. First, the masking means no individual *sᵢ* is leaked — only the sum across all honest signers is, and that sum is exactly what's needed to open the aggregate commitment. Second, the input to ZeroShare commits to *everyone's* commitment Cⱼ, so any divergence in views among the honest signers prevents the masks from canceling.

### Round 3 (online): masked response

Each signer computes the aggregate, recovers the aggregated nonce R, derives the Schnorr challenge *c*, and broadcasts a masked response:

```rust
fn sign_round3(
    sk_i: &SecretShare,
    state2: State2,
    received: &[Round2Message],
) -> Round3Message {
    let s_agg: (Scalar, Scalar) = received.iter().map(|m| m.s_tilde).sum();
    let C_agg = state2.commitments.product();    // ∏ C_j, by homomorphism

    // Recover the committed value: R is uniquely determined by C_agg and s_agg.
    let R = bcj_recover(ck, C_agg, s_agg);
    assert_eq!(C_agg, bcj_commit(ck, R, s_agg)); // binding check

    let c = hash_to_scalar(b"chal", &vk, state2.msg, R);

    // Second zero-share, different domain tag, so it's a fresh input.
    let ctnt = encode(b"sign", state2.SS, state2.msg, state2.commitments, R);
    let delta = zero_share(&sk_i.seeds, state2.SS, &ctnt);  // ∈ Z_p

    let L_i = lagrange(&state2.SS, i);
    let z_i = c * L_i * sk_i.x_i + state2.r_i;
    let z_tilde_i = z_i + delta;

    Round3Message { z_tilde: z_tilde_i }
}
```

Aggregation is trivial: *z = Σ z̃ᵢ*, and *σ = (c, z)* is a standard Schnorr signature, verifying with R = G^z · X^(−c) as usual. The Δᵢ values cancel, leaving *z = c·x + Σ rᵢ*.

### Why the rounds collapse

Compared to Crackle's five rounds, the trick is that BCJ's homomorphic structure lets the commitments do double duty. There's no explicit "share your view, sign that you agree" round, because the act of opening the aggregate commitment is implicitly a consistency check — if any honest signer disagreed about the commitment set, their Δ̃ᵢ wouldn't combine into a valid opening, and the binding check in Round 3 would fail. The first round is offline because nothing in it depends on the message or the signer set.

## Why it's secure

The end goal is a reduction to discrete log: if an adversary can produce a forgery, the reduction can compute *x* from *X*. The reduction goes through SelfTargetDL [Bellare-Dai 2021], an interactive variant of DL where the adversary must produce (M, c, z) such that H(X, G^z · X^(−c), M) = c. SelfTargetDL is loss-Q tight to DL, which is where the (Q_H + Q_S) factor in the security bound comes from.

The reduction has to simulate signing without knowing *x*. The standard playbook is: simulate the Schnorr signature using HVZK (pick *c* and *z* first, derive R), program the random oracle so that H(vk, M, R) = c, and use this simulated signature in place of an honest one. The tension in threshold Schnorr is that you can only program the oracle on the *aggregate* R, but the adversary helps construct R from its own (possibly malicious) commitments. You need to extract the adversary's contribution before you can compute R, but you also need to commit to your own honest contribution before learning the adversary's.

> The subtle bit is the order of operations across fourteen hybrids. The simulator first extracts the adversary's nonces from its NIZK proofs (SIM-EXT lets it do this on-the-fly while simulating honest proofs). Then, when the *last* honest signer in a session takes its turn, the simulator knows everyone's nonces, computes the aggregate R, and programs H. For all the earlier honest signers, the simulator publishes uniformly random masked values — these look correct because the masking is information-theoretically hiding until the last share lands. Only the last share is computed for consistency. Adaptive corruption is survived by retroactively "explaining" the random tape via the NIZK's explainability and by reprogramming the zero-share oracle to back-fit the Δᵢ.

The product equivocation property of BCJ is what makes the commitments themselves simulatable. Ordinary equivocation lets you open a single commitment to anything; product equivocation says you can simulate N commitments, equivocate the *product* of any subset to a chosen value, *and* still equivocate up to N−1 of the individual commitments to chosen values, as long as at least one is left unopened individually. That last commitment absorbs the randomness needed to make the algebra consistent.

> Why does at least one have to stay unopened? The aggregate opening reveals z = β·x + Σαₜ, where the αₜ are the individual messages. If you open every αₜ individually, you've over-determined the system: z is now a linear function of public values, but in the simulation the αₜ for the unopened commitments are what carries the entropy that makes z look random. With N−1 unopened commitments at the start and one remaining at the end, you always have at least one source of fresh randomness available. The proof uses this to construct a bijection between the real and ideal joint distributions; without the constraint, the bijection fails to be full-rank.

The final reduction guesses which of the adversary's H queries corresponds to the forgery, programs that specific query through the SelfTargetDL challenger, and outputs the forged (c, z) as the SelfTargetDL solution. The guess succeeds with probability 1/(Q_H + Q_S), which is where the dominant tightness loss comes from.

## Performance

Concrete numbers for 128-bit security, 256-bit prime field, 32-byte group elements:

| Scheme | Rounds | Online comm/party | Offline comm/party | Assumption | Adaptive |
|---|---|---|---|---|---|
| FROST | 1+1 | 32 B | 64 B | AOMDL | half (T/2) |
| FROST (full) | 1+1 | 32 B | 64 B | AOMDL+LDVR+AGM | full |
| Sparkle | 2+1 | ~128 B | 32 B | AOMDL | half |
| Crackle | 4+1 | ~128 B | 32 B | DL | full |
| Gargos | 3+0 | ~256 B | — | DDH | full |
| **Rackle** | **2+1** | **96 B** | **~3 KB + 64 B** | **DL** | **full** |

Online communication is excellent: three field elements per party, which is actually less than Sparkle or Crackle. The cost lands almost entirely on the offline round, dominated by the Fischlin SLE proof at roughly 3 KB. For comparison, an Ed25519 signature is 64 bytes total, and a single FROST run is essentially three signatures' worth of traffic per party.

The crossover with FROST is where the security model bites. FROST is faster and smaller, but its full-adaptive security needs three assumptions stacked together (AOMDL, LDVR, and the AGM) and even then only via [CKK+25] in the strongest version. Sparkle is similar in shape to Rackle's online phase but only proven half-adaptive (the adversary corrupts ≤ T/2). If you actually need the strong guarantees Rackle provides — adaptive corruption of up to T−1 parties, no erasures, no AGM — there is no faster scheme.

The 3 KB offline blob is the elephant in the room. In throughput terms it's not bad: at typical signing rates and modern network speeds, transmitting 3 KB during preprocessing is invisible. In latency terms it's also fine because Round 1 doesn't depend on the message or signer set, so it can be run hours before signing. The pain shows up in storage and bookkeeping: you have to keep the preprocessing state until you know which signers and message you'll commit to, and if you preprocess for a session that never happens, you've burnt 3 KB of bandwidth.

## Implementation notes

**No production crate exists.** The closest existing building blocks in Rust are `arkworks` for elliptic curve and pairing-free group operations, and `frost-secp256k1` / `frost-ed25519` for reference threshold Schnorr implementations you'd want to compare APIs against. The BCJ commitment is a half-page of group arithmetic; the randomized Fischlin transform is harder to get right and there is, to my knowledge, no audited Rust implementation. A real deployment would have to build the Fischlin compiler from scratch on top of the Σ-protocol shown in the paper's Appendix C.

**Constant-time discipline matters in places you might not expect.** The masking values Δᵢ are derived by hashing seeds, so if your scalar arithmetic for *zᵢ = c · Lᵢ · xᵢ + rᵢ + Δᵢ* leaks via timing, you leak the share *xᵢ*. The Lagrange coefficients themselves are public, but their multiplication with *xᵢ* must be constant-time. Don't reuse a scalar library that branches on operand bits.

**Domain separation is load-bearing.** The two zero-share invocations use different domain tags ("open" vs "sign" in my pseudocode, encoded as a leading bit in the paper), and the random oracles H_mask, H_c, H_aux are all formally separate. If you collapse them into a single hash with naive prefixing, you can break the proof. The paper's input encoding ctnt = "open"‖SS‖M‖(Cⱼ)ⱼ is also order-sensitive: signers must canonicalize the commitment list (sort by signer index) before hashing, or two honest signers with different iteration orders produce different Δᵢ values and the session fails.

**The Fischlin proof needs a hash with explicit output-length control.** The transform runs a proof-of-work where you repeatedly sample challenges and accept only when H(...) outputs b bits of zero. The paper sets b=8 and r=⌈λ/b⌉=16, giving each honest prover an expected ~16 × 256 = 4096 hash invocations. SHA-256 with output truncation is fine, but you need the API to expose "give me the first 8 bits" cleanly.

**The NIZK label is the signer index.** This is critical and easy to miss when porting code. The label ensures that if a commitment Cᵢ from signer *i* is reused or replayed by a malicious signer *j*, the proof π associated with Cⱼ uses label *j* and is still independently extractable. Drop the label and the security proof's "impersonation" handling collapses.

**Serialization gotchas around point compression.** BCJ commitments are pairs of group elements. If your serialization is variable-length (e.g., some libraries omit the leading byte for compressed points on certain curves), the hash inputs across signers can diverge silently. Pin a fixed-width encoding and test cross-implementation.

## Limitations

Rackle does not provide identifiable abort. If signing fails — a malicious participant sends garbage, a NIZK doesn't verify, the binding check rejects — you know the session failed but not which signer to blame. The paper conjectures the [PKN+25] interactive abort-identification protocol can be bolted on, but that's future work, and in a real deployment with N=20 or N=100 signers, "someone misbehaved, retry without them all" is a denial-of-service vector.

The 3 KB offline cost will dominate any deployment where signers run on bandwidth-constrained links or where preprocessing can't be amortized — IoT, embedded HSMs, mobile-to-mobile MPC. The communication is asymptotic in the security parameter (the Fischlin proof grows linearly in λ), so there's no easy way to shrink it without changing the underlying SLE proof technique. If SLE proofs get cheaper — and there's active work on this — Rackle inherits the improvement directly.

The DL assumption is minimal, but the proof goes through SelfTargetDL with a loss of √(Q_H · ε_DL), which is the standard Pointcheval-Stern style loss. For 128-bit security against an adversary with 2⁶⁰ hash queries, you want a curve where DL is 188-bit hard — secp256k1 is on the edge, Ristretto/Curve25519 is fine, the NIST P-curves are fine. This isn't a Rackle-specific concern (all single-Schnorr ROM proofs have it) but it's worth knowing.

The "no secure erasures" guarantee is strong, but the threat model still assumes signers don't leak intermediate state during the session itself — corruption is modeled as a discrete event between protocol messages, not as a continuous side channel. If your threat model includes a powerful local attacker watching memory while a signer is mid-computation, you need orthogonal defenses.

If you're building threshold Schnorr today and you can tolerate the AGM or the algebraic-OMDL assumption, FROST in two rounds is faster, smaller, well-deployed, and increasingly well-understood; reach for Rackle when you specifically need the combination of three rounds, full adaptive corruption, and a security argument that won't fall apart if AOMDL turns out to be subtly broken.
