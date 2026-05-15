# Aggregating Falcon Signatures with LaBRADOR
*How to squash 1000 post-quantum signatures into 120 KB while preserving NIST compatibility*

**Paper:** https://hal.science/hal-04700114/document

## TL;DR

Falcon is one of NIST's three post-quantum signature standards, and a single compressed signature is around 666 bytes. If you're running something like a consensus layer that bundles thousands of signatures into a block, that adds up fast — 1000 signatures is 666 KB of pure signature overhead. This paper shows how to use LaBRADOR, a recent lattice SNARK with quasi-optimal proof sizes, to aggregate N Falcon signatures into a single compact proof: 120 KB for N=1000, 417 KB for N=8192. Along the way the authors had to invent a new soundness framework called *Predicate Special Soundness* because LaBRADOR's structure doesn't fit cleanly into existing Fiat-Shamir analyses, and they had to make several careful surgical modifications to LaBRADOR itself so that Falcon's verification equation could be expressed in its native constraint language. The result is provably secure under standard lattice assumptions (Module-SIS) and works with the actual NIST-standardized Falcon, not a custom variant.

## Math you'll need

### Polynomial rings and the splitting trick

Falcon, LaBRADOR, and basically every modern lattice scheme work over polynomial rings rather than plain integers mod q. The ring is $R_q = \mathbb{Z}_q[X]/(X^d + 1)$, where d is a power of two (512 or 1024 for Falcon). An element of $R_q$ is a polynomial of degree less than d with coefficients in $\mathbb{Z}_q$.

The reason this beats working over $\mathbb{Z}_q$ directly is twofold. First, multiplication in $R_q$ via NTT is O(d log d), which makes polynomial dot products fast. Second, by the Chinese Remainder Theorem, $R_q$ factors into a product of smaller rings depending on how $X^d + 1$ splits modulo q. If $X^d + 1$ factors into l irreducibles, we say $R_q$ is "l-splitting." Fully-splitting (l = d) gives the fastest NTTs but the smallest factors; two-splitting (l = 2) is the opposite extreme. Falcon uses q = 12289 specifically because $X^d + 1$ factors very nicely there.

For our purposes, just remember: a polynomial ring is what you commit to when you do lattice ZK, and how it splits affects both performance and the soundness analysis. The paper picks two-splitting for the proof system because it gives an exponentially large challenge space where every nonzero element of small norm is invertible — a property the soundness proof leans on heavily.

### Module-SIS

Module-SIS (M-SIS) is the lattice problem that underpins binding for the commitments in LaBRADOR. Given a random matrix $A \in R_q^{n \times m}$, find a short vector $\vec{x} \in R_q^m$ with $A\vec{x} = \vec{0}$ and $0 < \|\vec{x}\|_2 \leq \beta$ (where β is the "shortness" bound, in this paper typically a few thousand to a few million depending on context).

If you've seen plain SIS, this is the structured version that lives over a polynomial ring. Hardness reduces (in the worst case) to finding short bases in module lattices, which is believed quantum-hard. Concretely, breaking M-SIS for the parameters in this paper requires running BKZ with block sizes well above 400 — comparable to the assumed security of Kyber and Dilithium.

Ajtai commitments are just M-SIS instances: to commit to a short vector $\vec{w}$, publish $\vec{v} = A\vec{w}$. The commitment is binding because two distinct short openings to the same $\vec{v}$ would give a short M-SIS solution $\vec{w}_1 - \vec{w}_2$. It's not hiding by default, but that doesn't matter for SNARKs.

### Falcon's verification equation

A Falcon signature on a message m, under public key $h \in R_q$, is a tuple $(r, s_1, s_2)$ where r is a random 320-bit salt. Verification has exactly two checks:

$$s_1 + h \cdot s_2 = H(r, m) \pmod{q}$$
$$\|(s_1, s_2)\|_2 \leq \beta$$

That's it. The first check is a single linear equation over $R_q$. The second is a Euclidean norm bound where β is around 5800 for Falcon-512. Compressed signatures omit $s_1$ since it can be recovered from $s_2$ and $H(r,m)$, but for aggregation we'll need both.

The beauty here is that this is essentially the same shape as LaBRADOR's native constraint language — one linear equation plus one norm bound. That's not an accident. The authors picked Falcon partly because the fit is so close.

### LaBRADOR's relation

LaBRADOR proves statements of the form "I know witness vectors $\vec{w}_1, \ldots, \vec{w}_r \in R_q^n$ satisfying a system of dot product constraints and a global norm bound." A dot product constraint looks like:

$$f(\vec{w}_1, \ldots, \vec{w}_r) = \sum_{i,j} a_{i,j} \langle \vec{w}_i, \vec{w}_j \rangle + \sum_i \langle \vec{\varphi}_i, \vec{w}_i \rangle - b = 0$$

with the global check $\sum_i \|\vec{w}_i\|_2^2 \leq \beta^2$. The parameter r is called the *multiplicity* and n the *rank*. Performance depends on the ratio between them — LaBRADOR works best when they're roughly balanced, which is the source of one of the adaptations later.

The protocol is recursive. Each iteration replaces a statement-witness pair with a smaller statement-witness pair, and the prover keeps recursing until further folding stops helping. The final proof is a logarithmic stack of intermediate messages plus a small base-case opening. Proof sizes scale roughly as O(polylog) in the witness, not log — but in practice the constant matters more than the asymptotic.

## The construction

### What's in an aggregate signature

An aggregator takes N raw Falcon signatures $\{(r_i, s_{i,1}, s_{i,2})\}_{i=1}^N$ on messages $m_i$ under public keys $h_i$, and produces a single compact proof $\pi$. The aggregate signature is $(\pi, r_1, \ldots, r_N)$ — yes, the salts have to ride along, because the verifier needs them to compute the hash targets $t_i = H(r_i, m_i)$ locally.

```rust
struct AggregateSignature {
    proof: LabradorProof,
    salts: Vec<[u8; 40]>, // 320-bit salts, one per aggregated sig
}

fn aggregate(
    sigs: &[(PublicKey, Message, FalconSignature)],
) -> AggregateSignature {
    let salts: Vec<_> = sigs.iter().map(|(_, _, s)| s.salt).collect();
    let targets: Vec<RingElem> = sigs.iter()
        .map(|(_, m, s)| hash_to_ring(&s.salt, m))
        .collect();

    // Witness: the (s_1, s_2) pairs we want to prove knowledge of.
    let witness: Vec<(RingElem, RingElem)> = sigs.iter()
        .map(|(_, _, s)| (s.s1.clone(), s.s2.clone()))
        .collect();

    // Statement: public keys, targets, norm bound. Invariant for each i:
    //   s1_i + h_i * s2_i = t_i (mod q)   AND   ||(s1_i, s2_i)||_2 <= beta.
    let statement = build_labrador_statement(&sigs, &targets);

    let proof = labrador_prove(&statement, &witness);
    AggregateSignature { proof, salts }
}
```

Note the salts make the size technically linear in N, so the construction isn't asymptotically succinct. In practice the 40-byte salts are dwarfed by the rest of the signature data for any meaningful N, and the authors point out that deterministic or "synchronized" Falcon variants would let you drop them entirely. For now, just accept that there's a per-signature 40-byte tax.

### Adaptation 1: a bigger modulus

Falcon uses q = 12289. LaBRADOR can't. The reason is the Johnson-Lindenstrauss projection inside LaBRADOR's norm check: to prove a witness has small Euclidean norm, the verifier sends a random ±1-valued projection matrix and the prover replies with the projected vector. Soundness requires $\sqrt{\lambda} \cdot \beta \leq q / C_1$ for a constant $C_1 \approx 120$ at λ = 128 (security parameter). Falcon's q is way too small to satisfy this.

The fix is to use a larger modulus q' for the proof system (40–60 bits depending on N), and to embed Falcon's verification equation by lifting it to the integers:

$$s_{i,1} + h_i \cdot s_{i,2} + q v_i - t_i = 0 \in \mathbb{Z}[X]/(X^d + 1)$$

The new witness element $v_i$ captures the modular wraparound that originally happened mod q. To enforce that the equation actually holds over the integers (not just mod q'), the prover proves an infinity-norm bound tight enough that no wraparound mod q' could have occurred. This is itself done with a Johnson-Lindenstrauss projection, since the infinity norm is upper-bounded by the Euclidean norm.

### Adaptation 2: exact norm proofs via four squares

LaBRADOR's built-in norm check has slack: if you prove $\|\vec{w}\|_2 \leq \beta$, the extractor only recovers a witness of norm up to roughly $\sqrt{\lambda / C_2} \cdot \beta \approx 2\beta$. For internal recursion that's fine, but for our application we need the extracted witness to be a *valid* Falcon signature, which means norm exactly bounded by β. Otherwise the security reduction can't claim to have produced a real forgery.

The trick is Lagrange's four-square theorem: every non-negative integer is a sum of four squares. So if $\beta^2 - \|(s_{i,1}, s_{i,2})\|_2^2 \geq 0$, there exist integers $\epsilon_{i,0}, \ldots, \epsilon_{i,3}$ with:

$$\beta^2 - \|s_{i,1}\|_2^2 - \|s_{i,2}\|_2^2 = \epsilon_{i,0}^2 + \epsilon_{i,1}^2 + \epsilon_{i,2}^2 + \epsilon_{i,3}^2$$

The prover adds the epsilons to the witness as coefficients of a degree-3 polynomial, and the four-square identity becomes a quadratic constraint expressible in LaBRADOR's dot-product language. This is the [GHL22] four-square trick, originally from non-interactive verifiable secret sharing. The catch: the constraint only holds *modulo q'*, so we again need an infinity-norm bound on the witness coefficients to ensure no wraparound could have made the sum spuriously equal $\beta^2$ over the reals.

### Adaptation 3: reshape the witness

A naive encoding gives r = 2N witness vectors of rank n = 1 — radically unbalanced, and LaBRADOR doesn't like that. Performance and proof size are sensitive to the ratio between multiplicity and rank.

The authors reshape: pack the $s_{i,j}$ into vectors of length N, getting roughly $r = O(\sqrt{N})$ witness vectors of rank n = N. This is mostly bookkeeping for the linear constraints, but the four-square constraints are quadratic and require a clever padding scheme so that $\|s_i\|_2^2$ still comes out of an inner product (the construction places each $s_i$ at a fixed offset and its conjugation-twisted partner at a matching offset). Appendix F has the gory details; the upshot is that the prover's runtime drops from $O(N^2)$ to $O(N^{1.5})$ and proof sizes shrink modestly.

### Adaptation 4: work over a subring

Falcon operates over a ring of degree d = 512 or 1024, but the proof sizes shrink considerably if LaBRADOR operates over a subring S of degree d' = d/c for some small c. There's a norm-preserving bijection $\phi: R_q^n \to S_q^{cn}$ that maps polynomial vectors to longer vectors over the smaller ring without changing the underlying constraints. After applying this, the proof sizes drop by at least 2×.

The choice of c interacts with the splitting behavior of S mod q'. The paper experiments with several configurations and concludes that c=8 (giving d'=64 for Falcon-512) combined with two-splitting and a separate "lifting" trick — multiplying in a larger NTT-friendly ring and reducing back — gives both the smallest proofs and the fastest arithmetic.

### Putting it together

```rust
// Roughly what the aggregator does, after all four adaptations.
fn labrador_prove(stmt: &Stmt, w: &Witness) -> LabradorProof {
    // Witness now lives in the subring S_{q'} with degree d' = d/c.
    let (witness_subring, statement_subring) = lift_to_subring(w, stmt);

    // Augment witness with v_i (wraparound), epsilon_i (four-squares),
    // and the conjugation-twisted copies needed to express ||.||_2 as
    // an inner product inside LaBRADOR's language.
    let augmented = augment_witness(witness_subring);

    // Build the constraint system: Falcon verification per signature,
    // four-square equation per signature, padding-consistency checks,
    // and infinity-norm projections for both wraparound bounds.
    let constraints = build_constraints(&augmented, &statement_subring);

    // Run base LaBRADOR, then recurse log(log N) times until folding
    // no longer helps. Each iteration is commit -> project -> aggregate
    // -> amortize, with a Fiat-Shamir transcript binding the messages.
    labrador_base(constraints, augmented).recurse_until_base_case()
}
```

## Why it's secure

Unforgeability of the aggregate signature reduces to two things: unforgeability of plain Falcon, and knowledge soundness of the underlying argument. The former is NIST-blessed. The latter is where the paper does its real work.

The high-level reduction is standard: a forgery on the aggregate means the prover produced a valid proof for a statement that includes some message $m^*$ never signed by the legitimate key. Extract from the prover; one of the extracted "signatures" must be a forgery on $m^*$. Appendix C handles a subtle gap noted by Fiore-Nitulescu — that a powerful signing oracle could in principle interfere with extraction — by proving that hash-then-sign signing oracles don't, as long as the SNARK's knowledge soundness holds against an adversary given auxiliary input consisting of random preimage-image pairs. But the real devil is in extraction itself, which is where the new PSS framework comes in.

### Why standard special soundness isn't enough

The standard tool for analyzing Fiat-Shamir SNARKs is the AFK22 result on K-tree special soundness: if your interactive protocol admits extraction from any K-tree of accepting transcripts, then its Fiat-Shamir transform has knowledge error roughly Q·κ where Q is the prover's RO query budget and κ is the interactive knowledge error.

LaBRADOR breaks this framework in three places. First, its challenge differences aren't always invertible in $R_{q'}$ — in a high-splitting ring you only get invertibility with overwhelming probability. Second, extraction doesn't immediately produce a witness; it produces a *candidate* witness that the extractor then has to verify is short via additional probabilistic checks. Third, the amortization phase uses challenge *vectors*, and extraction needs vectors that differ in only one coordinate at a time.

> The second issue is the most painful. You'd like to set k (the number of distinct transcripts you collect at one tree level) large enough that with overwhelming probability a fresh independent challenge would land in the "good" set for whatever invariant you're enforcing — say, that the witness is short. But the Johnson-Lindenstrauss soundness means the fraction of bad challenges is $2^{-\lambda}$, and to be sure of avoiding them by sampling you'd need k exponential in λ. Then the extractor isn't polynomial-time anymore. The way out is to extract the candidate witness from a small tree, then run *one extra* fresh probabilistic check against it. AFK22 doesn't model this.

### Predicate Special Soundness

The PSS framework attaches to each level of the transcript tree two kinds of predicates: *challenge predicates*, which enforce conditions on the challenges sampled at that level (e.g. "this challenge difference is invertible"), and *commitment predicates*, which enforce conditions on the values extracted from sub-trees (e.g. "the extracted witness is short, given an independent projection challenge"). Each predicate comes with a *failure density* — an upper bound on the fraction of challenges that let a cheating prover pass without satisfying the predicate.

The headline result is that the Fiat-Shamir knowledge error decomposes neatly as a sum over levels of failure densities, scaled by Q. For LaBRADOR concretely, the per-iteration knowledge error is roughly:

$$2(Q+1) \cdot \left(2^{-\lambda} + q^{-d/l} + q^{-\lceil \lambda/\log q \rceil} + (5 + 2l) \cdot r \cdot B\right)$$

where l is the splitting factor, r is the witness multiplicity, and B is the "well-spreadness" of the challenge space — essentially the probability that a random challenge takes any specific value in any CRT slot. When you set B to be $2^{-\lambda}$ via the challenge-space design, every term is negligible. Across the log-log many recursive iterations, the error sums additively.

> The other subtle bit is *coordinate-wise* extraction in the amortization round. The prover replies with $\vec{z} = \sum_i c_i \vec{w}_i$ for a vector of challenges $(c_1, \ldots, c_r)$. To extract $\vec{w}_i$, the extractor needs two transcripts that share all challenges except $c_i$. That's a coordinate-wise tree of transcripts — a structure introduced for related reasons in recent work on coordinate-wise special soundness. PSS is extended to handle this; the proof in Appendix G is mostly bookkeeping but it's where you'd look if you wanted to adapt the framework to another protocol with amortization.

The binding part of the security argument is more conventional: any failure of a commitment predicate either implies the candidate witness was valid all along, or yields two distinct openings to the same Ajtai commitment, which is an M-SIS solution. The norms work out so that the M-SIS challenge has bound around $8 T_{op} (b+1) \beta'$, where $T_{op}$ is the operator norm of a challenge — still well within the hard regime.

## Performance

All numbers are for Falcon-512, targeting 121-bit security (Falcon-512 itself is 128-bit; the small loss comes from the SNARK's knowledge error).

| N | Trivial concat | Phoenix | Squirrel | Chipmunk | This paper (with salts) | This paper (no salts) |
|---|---|---|---|---|---|---|
| 500 | 333 KB | 3,616 KB | — | — | 93 KB | 73 KB |
| 1,000 | 666 KB | 7,444 KB | 572 KB | 118 KB | 122 KB | 81 KB |
| 2,000 | 1,332 KB | 15,319 KB | — | — | 165 KB | 85 KB |
| 4,096 | 2,729 KB | — | 693 KB | — | 252 KB | 88 KB |
| 8,192 | 5,458 KB | — | 762 KB | 160 KB | 417 KB | 89 KB |

A few crossover points to internalize. Against trivial concatenation, the aggregate proof becomes smaller than the naive bundle at around N ≈ 110 signatures. Below that, just send the signatures. Against Phoenix — the only other lattice aggregate signature scheme with a comparable security model — this construction wins by one to two orders of magnitude at every size, and Phoenix actually gets *worse* beyond N = 4000 while LaBRADOR keeps scaling.

The interesting comparison is against Chipmunk, which is smaller for large N. Chipmunk wins on raw size but operates in the *synchronized* model (one signature per time period per key, bounded total life cycle) and uses its own Merkle-tree-based signature scheme rather than NIST-standardized Falcon. If you drop the salts — which you can do if you're willing to move to a synchronized or deterministic Falcon variant — this construction beats Chipmunk at every size: 81 KB vs 118 KB at N=1024, 89 KB vs 160 KB at N=8192.

The unmentioned cost is prover time. The paper doesn't publish wall-clock numbers for the full aggregation, only microbenchmarks for the polynomial arithmetic underneath. The two-splitting-with-lifting strategy runs polynomial multiplication in about 1,090 cycles for Falcon-512's d', which is roughly an order of magnitude slower than equivalent operations in fully-splitting NTT-friendly rings. Expect the aggregator to be the slowest part of any system that uses this.

## Implementation notes

**Existing code.** There is no production Rust implementation. The authors publish parameter estimation scripts and proof-of-concept benchmarks at github.com/dfaranha/aggregate-falcon, mostly in Python and C. Beullens and Seiler's original LaBRADOR reference is in C and not designed for production. If you're building this for real, you'd be starting from the LaBRADOR base protocol and porting it carefully — you can borrow NTT machinery from existing Falcon implementations like pqcrypto-falcon or pqclean, but the lattice SNARK layer is bespoke.

**Constant-time considerations.** The aggregator doesn't touch any secret keys — it only manipulates already-issued signatures, which are public from the perspective of the original signer. So constant-time discipline for the aggregator is less critical than for Falcon's signing routine. However, the Fiat-Shamir challenges depend on the witness via the commitments, and any timing leak in the polynomial arithmetic on the witness side could theoretically reveal something about the signatures being aggregated. In a setting where the set of aggregated signatures is itself sensitive — private mempools, anonymous credential bundles — implement the prover constant-time. Verification is on entirely public data and constant-time doesn't matter.

**Domain separation.** The construction uses two distinct random oracles: H for Falcon's hash-to-ring, and a separate transcript hash inside LaBRADOR. These must be domain-separated (different prefixes or different SHAKE personalizations) or you'll have transcript collisions across the layers. The salts must be hashed *together with the message* in Falcon's H — this matches Falcon's spec but it's easy to break when wiring things up. Test against the NIST KAT vectors before doing anything else.

**Serialization.** The proof contains a mix of polynomial commitments (vectors in $R_{q'}^{\kappa_1}$), short responses (vectors in $R_{q'}^n$ with provably small coefficients), and challenge seeds (32 bytes each, one per Fiat-Shamir round). Pack the short responses bit-tightly using the known norm bound — naively encoding each coefficient as a full $\log q' \approx 60$-bit value will double your proof size. Reuse Falcon's compress/decompress routines for inspiration; the math is similar.

**Side channels for the verifier.** Verification involves polynomial multiplications mod q'. Use Montgomery or Barrett reduction in constant time even though there's no secret here — it's good hygiene and makes the verifier robust if it's deployed somewhere weird like a TEE or an attestation enclave.

## Limitations

The salts ruin asymptotic succinctness. Each aggregated signature contributes 40 bytes of salt to the final aggregate, so the construction is technically linear in N. For N = 8192 the salts contribute around 320 KB out of 417 KB total — they're the dominant cost. The synchronized variant fixes this but loses standardization compatibility, which was the whole point of picking Falcon.

The modulus q' is much bigger than Falcon's q, which means the aggregator's memory footprint per signature is several times larger than the signature itself. For N = 10,000 you need a 61-bit q' for Falcon-512 and a 63-bit q' for Falcon-1024 — right at the edge of fitting in a single machine word. Beyond that you'd need to move to double-word arithmetic, which is roughly 7× slower per multiplication.

The security loss from the PSS analysis costs a few bits — the paper reports Falcon-512 aggregates at 121-bit security rather than Falcon-512's nominal 128. For Falcon-1024 the loss is symmetric, ending at 249-bit. This is fine for most applications but worth noting if you're targeting a specific security level for compliance reasons.

Finally, the heaviest open question is prover performance in absolute terms. The paper benchmarks polynomial arithmetic but doesn't publish end-to-end aggregation times for N in the thousands. Lattice SNARKs are generally slow to prove, and LaBRADOR is no exception. If you're aggregating signatures in latency-sensitive contexts — say, on the critical path of block production — measure carefully before committing.

Reach for this construction when you have hundreds to thousands of independently generated Falcon signatures, you need to ship them as a single blob, you can afford a heavy aggregator, and asymptotic succinctness isn't required; reach for Chipmunk if you control the signature scheme and can live in the synchronized model; reach for plain concatenation if N is under about 100.
