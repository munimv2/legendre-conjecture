# Legendre's conjecture in October 2026

**What is proved, what has been computed, and what is open.** This is an exposition with complete
proofs of the elementary reductions and a small, tested Python companion (`legendre.py`) that
re-derives the classical computational range from a prime-gap table, checks Oppermann's stronger form
by sieve, finds witnesses for any *n* up to 1.8 · 10¹², and reproduces the 2026 Mills-type constant
attached to the conjecture.

By MV2 of Munim, Inc. This article and the code were researched and written with AI agents working in
MV2's name. Every statement below is either proved here or cited to a source; section 9 separates
what is proved from what is open, and section 11 lists what is new here. Corrections are welcome by
issue or pull request.

Version 1.0, 2026-10-14. Text CC BY 4.0, code MIT.

Related problem records: Justin Sun Prize JSP-000012 (MathDB #316127), OEIS A014085, A144140.

---

## 1. The question

**Legendre's conjecture.** For every positive integer *n* there is a prime strictly between *n*² and
(*n* + 1)².

Legendre stated it in the second edition of his *Essai sur la théorie des nombres* (1808). At the 1912
International Congress of Mathematicians, Landau listed it among four problems on primes he called
"unattackable" with the methods of the time; the others are Goldbach's conjecture, the twin prime
conjecture and the infinitude of primes of the form *n*² + 1. All four are still open.

Writing *a*(*n*) for the number of primes in (*n*², (*n* + 1)²) (OEIS A014085), the conjecture says
*a*(*n*) ≥ 1 for all *n* ≥ 1. The data look nothing like a borderline case:

    n     0  1  2  3  4  5  6  7  8  9  10 ...  79  80
    a(n)  0  2  2  2  3  2  4  3  4  3  5  ...  23  13

The interval has length 2*n* + 1 and sits near *n*², where primes have density about 1/ln(*n*²), so
one expects about *n*/ln *n* primes in it. The difficulty is not that the conjecture is nearly false;
it is that no method controls primes in *every* interval as short as the square root of its location
(section 3).

Two strengthenings are often used:

- **Oppermann's conjecture:** for every *n* ≥ 2 there is a prime in (*n*², *n*(*n* + 1)) and a
  prime in (*n*(*n* + 1), (*n* + 1)²). It implies *a*(*n*) ≥ 2.
- **Real version** (Chamberland–Straub 2026): every interval [*x*², (*x* + 1)²] with real *x* ≥ 1
  contains a prime.

## 2. Three elementary reductions, with proofs

**Lemma 1 (gap tables settle ranges).** Let *X* ≥ 4 and *G* ≥ 1 be such that every prime *p* ≤ *X* is
followed by the next prime within *G*. Then (*n*², (*n* + 1)²) contains a prime for every *n* with
max(2, ⌈*G*/2⌉) ≤ *n* ≤ ⌊√*X*⌋.

*Proof.* Since *n* ≥ 2, there is a largest prime *p* ≤ *n*², and *p* ≤ *n*² ≤ *X*. The next prime *p*′
satisfies *p*′ > *n*² (by maximality of *p*) and *p*′ ≤ *p* + *G* ≤ *n*² + 2*n* < (*n* + 1)². ∎

`legendre_range_from_gaps(X, G)` returns this range. With the computation of all prime gaps up to
4 · 10¹⁸ by Oliveira e Silva, Herzog and Pardi (2014), whose largest gap is 1476 (after the prime
1425172824437699411, found by Oliveira e Silva in 2009), Lemma 1 gives 738 ≤ *n* ≤ 2 · 10⁹, and a sieve
covers *n* < 738 (the test suite checks *n* ≤ 2000). That is the range 2 · 10⁹ that Sorenson and Webster
cite as the state of verification before their work (section 6).

**Lemma 2 (Andrica ⇒ Legendre).** Andrica's conjecture says √*p*′ − √*p* < 1 for consecutive primes
*p* < *p*′. It implies Legendre's conjecture.

*Proof.* If (*n*², (*n* + 1)²) had no prime for some *n* ≥ 2, take the largest prime *p* ≤ *n*² and the
next prime *p*′ ≥ (*n* + 1)²; then √*p*′ − √*p* ≥ (*n* + 1) − *n* = 1. (For *n* = 1 the interval contains 2
and 3.) ∎

**Lemma 3 (exponents of gaps and of powers; Chamberland–Straub 2026, Prop. 1.3).** Fix 0 < θ < 1 and put
α = 1/(1 − θ). If *p*ₖ₊₁ − *p*ₖ < *p*ₖ^θ for all large *k*, then [*x*^α, (*x* + 1)^α] contains a prime for
all large real *x*.

*Proof.* If not, choose *k* with *p*ₖ < *x*^α < (*x* + 1)^α < *p*ₖ₊₁. Then
(*p*ₖ₊₁ − *p*ₖ)/*p*ₖ^θ > ((*x* + 1)^α − *x*^α)/*x*^{αθ} > α*x*^{α−1}/*x*^{αθ} = α > 1, since α(1 − θ) = 1. ∎

Legendre's conjecture is the case α = 2, which needs θ = 1/2 *with no logarithmic loss and a constant
below 2*. Conversely, Legendre for *n* and *n* + 1 gives gaps below 4√*p* + 4 (a prime *p* in
[*n*², (*n* + 1)²) has a successor below (*n* + 2)²). So the conjecture is essentially the statement
"prime gaps are O(√*p*)" with an explicit constant, and that is exactly the size where every known
method stops, even under the Riemann hypothesis (section 4).

## 3. What is proved unconditionally

**Gaps between consecutive primes** (*p*ₖ₊₁ − *p*ₖ < *p*ₖ^θ for large *k*):

| θ | Who | Consequence by Lemma 3 |
|---|---|---|
| 1 | Bertrand (conjecture 1845), Chebyshev (1852) | prime in (*n*, 2*n*) |
| 32999/33000 | Hoheisel (1930) | first θ < 1 |
| 5/8 + ε | Ingham (1937) | prime between *n*³ and (*n* + 1)³ for all large *n* |
| 0.535 | Baker–Harman (1996) | |
| 0.525 | Baker–Harman–Pintz (2001) | prime in [*x*^α, (*x* + 1)^α] for α ≥ 2.106 and large *x* |

Ingham's theorem is conditional on a bound ζ(1/2 + *it*) = O(*t*^*c*), giving any θ > (1 + 4*c*)/(2 + 4*c*);
with the classical *c* = 1/6 + ε this is θ > 5/8, and Mills (1947) built his prime-representing
constant on it. Baker, Harman and Pintz remark that their *x*₀ could be made explicit "with enough
effort"; nobody has done it. A 2023 preprint by Runbo Li (arXiv:2308.04458, several revisions)
claims θ = 0.52; we found no record of publication.

**Explicit consecutive powers, for every n ≥ 1.** A prime lies between *n*^α and (*n* + 1)^α for all
*n* ≥ 1 when α ≥ 5 · 10⁹ (Dudek 2016), α ≥ 1,438,989 (Mattner 2017), α ≥ 155 (Cully-Hugill 2023),
α ≥ 140 (Cully-Hugill–Johnston 2023) and **α ≥ 90** (Cully-Hugill–Johnston, Funct. Approx. 2025).

**Cubes for huge n.** A prime lies between *n*³ and (*n* + 1)³ for all *n* ≥ exp(exp(32.537))
(Cully-Hugill, PhD thesis 2023, Thm 4.3), improving Dudek's 2016 bound of the same shape.

**Almost all n.** Bazzanella (2000): for every ε > 0, all but O(*N*^{1/4+ε}) of the intervals
[*n*², (*n* + 1)²] ⊆ [1, *N*] contain a prime, and indeed the expected number of primes. Much shorter
intervals almost always contain primes: [*n*, *n* + *n*^θ] with θ = 1/20 + ε for all but o(*N*) of
*n* ∈ [*N*, 2*N*] (Jia 1996), and with θ = 2/15 + ε the expected number of primes (Guth–Maynard,
Annals of Mathematics, to appear 2026; arXiv:2405.20552).

**Almost primes between all squares.** Chen (1975): between sufficiently large consecutive squares
there is an integer with at most two prime factors. Removing "sufficiently large":

- Dudek–Johnston, *Almost primes between all squares*, J. Number Theory 278 (2026) 726–745
  (arXiv:2501.18048): for **every** *n* ≥ 1 there is an integer in (*n*², (*n* + 1)²) with at most
  **four** prime factors counted with multiplicity.
- P. J. Campbell, arXiv:2603.10356 (March 2026, revised May 2026): at most **three** prime factors,
  for every *n* ≥ 1, using the Sorenson–Webster computation and maximal gap bounds for *n*² ≤ 10³¹ and an
  explicit linear sieve with Richert's logarithmic weights above. The arXiv page lists a Journal of
  Number Theory DOI (10.1016/j.jnt.2026.07.014).
- Johnston–Sorenson–Thomas–Webster, arXiv:2601.15564 (2026): between **cubes**, an integer with at most
  **two** prime factors for every *n* ≥ 1.

**A prime-representing constant without the conjecture.** Matomäki (2010) proved that some *A* makes
⌊*A*^{2ⁿ}⌋ prime for every *n* ≥ 1, the conclusion Legendre's conjecture would give (section 7), by a
different construction that gives no value of *A*.

## 4. What is proved under hypotheses

- **Riemann hypothesis.** Cramér (1936): gaps are O(√*p* log *p*), a logarithm too weak for Lemma 1.
  The best explicit version is Carneiro–Milinovich–Soundararajan (2019): on RH there is a prime in
  [*x*, *x* + (22/25)√*x* log *x*] for all *x* ≥ 4. On RH, [*x*^α, (*x* + 1)^α] contains a prime for all
  large *x* whenever α > 2; Ingham's theorem shows the Lindelöf hypothesis already suffices.
- **Chamberland–Straub (arXiv:2602.22502, February 2026).** On RH, [*x*^{2+δ}, (*x* + 1)^{2+δ}] contains a
  prime for **all real x ≥ 1** when δ ≥ 1/4 (the method allows δ down to 0.2253), and for
  *x* ≥ [β(log β + log log β)]^β, β = 2/δ, when δ < 1/4. The proof joins the explicit RH interval above
  (covering *x* ≳ 6.55 · 10¹⁰ at δ = 1/4) with the Sorenson–Webster computation (covering
  *x* < 2.04 · 10¹²). They judge δ ≥ 1/5 out of reach with current theory and computation. Caldwell and
  Cheng (2005) had shown the case δ = 1 (cubes, every *n* ≥ 1) on RH.
- **Cramér's conjecture** (gaps O((log *p*)²)) would imply Legendre's conjecture for all large *n*; the
  gap computations show *p*ₖ₊₁ − *p*ₖ < (log *p*ₖ)² for 11 ≤ *p*ₖ ≤ 4 · 10¹⁸.

So even RH leaves Legendre's conjecture open: it gives exponent 2 + δ for any fixed δ > 0, not 2.

## 5. Smaller exponents: where the heuristics stop

For 1 < α < 2 one can still ask whether [*x*^α, (*x* + 1)^α] contains a prime for large *x*
(Chamberland–Straub, Conjecture 4.2, equivalent to gaps O(*p*^ε) for every ε > 0; open even on RH). For
α = 3/2 the integers *n* with no prime in [*n*^{3/2}, (*n* + 1)^{3/2}] appear to be exactly

    10, 20, 24, 27, 32, 65, 121, 139, 141, 187, 306, 321, 348, 1006, 1051   (OEIS A144140)

`three_halves_gaps(N)` tests this with integers only (*p* ≥ *n*^{3/2} ⟺ *p*² ≥ *n*³), and the test suite
checks the list for *n* ≤ 2000.

## 6. What has been computed

| Statement | Range | Source |
|---|---|---|
| All prime gaps; largest 1476 | primes up to 4 · 10¹⁸ | Oliveira e Silva, Herzog, Pardi (Math. Comp. 2014) |
| Legendre, via Lemma 1 | *n* ≤ 2 · 10⁹ | from the row above |
| **Oppermann** (both halves) | **2 ≤ *n* ≤ 7.05 · 10¹³**, i.e. *n*² ≤ 4.97 · 10²⁷ | Sorenson–Webster, Research in Number Theory 11 (2025); arXiv:2401.13753 (earlier version to 3.33 · 10¹³) |
| Prime between consecutive cubes | *n*³ ≤ 1.649 · 10⁴⁰ | Johnston–Sorenson–Thomas–Webster (2026) |

Sorenson and Webster search for a prime in each half of every interval with a sieve-and-test algorithm
whose running time, under a heuristic based on Cramér's conjecture, is O(*N* log *N* log log *N*); the
computation found 27-digit primes. By Chamberland–Straub's Proposition 3.4, their result also gives
the real version for 1 ≤ *x* < 7.05 · 10¹³.

`legendre.py witness n` finds the two Oppermann witnesses for any single *n* up to 1.8 · 10¹² with a
Miller–Rabin test that is deterministic below ψ₁₃ = 3317044064679887385961981 (the smallest strong
pseudoprime to the first 13 prime bases, found by the same Sorenson and Webster). This checks
individual cases; it is not a verification of a range.

## 7. A concrete face of the conjecture: Mills-type constants

**Lemma 4 (Niven's construction; Chamberland–Straub Lemma 5.1).** Let α > 1 and let *q*₁, *q*₂, ... be
primes with *q*ₖ^α ≤ *q*ₖ₊₁ < (*q*ₖ + 1)^α − 1. Then *A* = lim *q*ₖ^{α^{−k}} exists and ⌊*A*^{αᵏ}⌋ = *q*ₖ for
all *k*.

*Proof.* Raising *q*ₖ^α ≤ *q*ₖ₊₁ < *q*ₖ₊₁ + 1 < (*q*ₖ + 1)^α to the power α^{−(k+1)} gives

    q_k^(α^-k) ≤ q_{k+1}^(α^-(k+1)) < (q_{k+1} + 1)^(α^-(k+1)) < (q_k + 1)^(α^-k).

The left sequence increases and is bounded by the right one, which decreases, so the limit *A* exists
and *q*ₖ^{α^{−k}} ≤ *A* < (*q*ₖ + 1)^{α^{−k}} for every *k*, i.e. *q*ₖ ≤ *A*^{αᵏ} < *q*ₖ + 1. ∎

For α = 2 the condition is a prime in [*q*ₖ², (*q*ₖ + 1)² − 1), which Legendre's conjecture supplies.
Choosing *q*₁ = 2 and each *q*ₖ₊₁ as the least prime above *q*ₖ² gives (Chamberland–Straub, Example 5.3)

    q = 2, 5, 29, 853, 727613, 529420677791, ...      A = 1.524699960538094359923363575688...

and *A* is then the smallest constant with ⌊*A*^{2ᵏ}⌋ prime for all *k*, *if* the construction never
breaks. Whether it breaks is **Conjecture 5.4** of that paper: *q*ₖ₊₁ < (*q*ₖ + 1)² for all *k*. It is
weaker than Legendre's conjecture but, in the authors' words, "likely of a similar difficulty".

The excesses *q*ₖ₊₁ − *q*ₖ² are 1, 4, 12, 4, 22 for the listed terms, against an allowance of 2*q*ₖ.
`python legendre.py mills 7` computes the seventh term (about 2.8 · 10²³; deterministic primality
applies) and brackets *A* to within 10⁻²⁵ using Lemma 4; the tests check that the published 30 digits lie
in that bracket and that ⌊*A*^{2ᵏ}⌋ = *q*ₖ for *k* ≤ 6. For cubes the same code reproduces Mills'
setting: *q* = 2, 11, 1361, 2521008887 and the Caldwell–Cheng constant 1.3063778838630806904686144926...
(on RH the smallest Mills constant), bracketed to within 6.4 · 10⁻¹² as stated in the paper.

## 8. Claimed proofs

Proofs of Legendre's conjecture are announced regularly in preprints, for example arXiv:2403.08306
(2024) and several 2024 preprints on Oppermann's conjecture. We found none accepted by the research
community, and every 2026 paper cited above treats the conjecture as open. A claimed proof should be
read against section 2: it must, in effect, control primes in every interval of length 2√*x* around *x*,
which is beyond what RH is known to give.

## 9. Proved versus open, as of 2026-10-14

**Proved:**

1. Lemmas 1–4 above (elementary).
2. A prime between *n*^α and (*n* + 1)^α for every *n* ≥ 1 when α ≥ 90; between cubes for
   *n* ≥ exp(exp(32.537)); for α ≥ 2.106 and all large *x* (non-explicit).
3. For every *n* ≥ 1, (*n*², (*n* + 1)²) contains an integer with at most four prime factors
   (published); at most three (2026 preprint, see section 3).
4. All but O(*N*^{1/4+ε}) of the intervals up to *N* contain a prime.
5. Computation: Oppermann's conjecture, hence Legendre's with *a*(*n*) ≥ 2, for 2 ≤ *n* ≤ 7.05 · 10¹³.
6. On RH: primes between *x*^{2+δ} and (*x* + 1)^{2+δ} for all *x* ≥ 1 when δ ≥ 1/4.
7. Unconditionally, some *A* has ⌊*A*^{2ᵏ}⌋ prime for all *k* (Matomäki), value unknown.

**Open:**

1. Legendre's conjecture, Oppermann's conjecture and Andrica's conjecture.
2. Legendre's conjecture assuming RH.
3. Primes between consecutive cubes for every *n* (true for *n*³ ≤ 1.649 · 10⁴⁰ and for
   *n* ≥ exp(exp(32.537)); the gap in between is open unconditionally).
4. An integer with at most two prime factors between all consecutive squares (Chen's theorem for every
   *n*).
5. Chamberland–Straub Conjecture 5.4, and Conjecture 4.2 for any α < 2.
6. An explicit *x*₀ in the Baker–Harman–Pintz theorem.

## 10. Why this matters beyond the problem

Legendre's conjecture is the cleanest test of how well we understand primes in short intervals. Each
unconditional advance listed in section 3 came from a new input (zero-density estimates, sieve
weights, large-value estimates for Dirichlet polynomials), and the 2025–2026 results show the field
moving on the explicit side: computation to 7 · 10¹³ made "for all *n* ≥ 1" statements possible for
almost primes and conditional power results.

## 11. What is original here

- Lemma 1 written out with explicit constants and turned into code, re-deriving the pre-2025
  verification range *n* ≤ 2 · 10⁹ from the 4 · 10¹⁸ gap table plus a sieve.
- Short proofs of Andrica ⇒ Legendre and of the nested-interval Lemma 4, and a dated map of the 2026
  literature (Dudek–Johnston, Campbell, Johnston–Sorenson–Thomas–Webster, Chamberland–Straub,
  Guth–Maynard).
- `legendre.py`: integer-only sieve and deterministic Miller–Rabin (below ψ₁₃), A014085 and its cube
  analogue, Oppermann counts and single-*n* witnesses, Lemma 1, the Mills-type sequences and brackets
  for squares and cubes, and the α = 3/2 exceptions.
- `test_legendre.py`: tests against OEIS A014085, A003617 and A144140, π(10⁴) and π(10⁶), known strong
  pseudoprimes ψ₄, ψ₉ = ψ₁₀ = ψ₁₁, ψ₁₂ and ψ₁₃, the published Mills-type sequences and constants, and
  Oppermann for 2 ≤ *n* ≤ 2000.

Nothing here claims progress on the conjecture.

## 12. Using the companion code

    python legendre.py count 100                 # a(n), OEIS A014085
    python legendre.py oppermann 10000           # both halves by sieve, 2 <= n <= N
    python legendre.py witness 1000000000000     # Oppermann witnesses for one n (n <= 1.8e12)
    python legendre.py gaps 4000000000000000000 1476
    python legendre.py mills 7                   # squares; add "3" for cubes (at most 4 terms)
    python legendre.py three-halves 2000         # OEIS A144140
    python -m unittest test_legendre             # or: python test_legendre.py

Python 3.6+, no dependencies.

## Sources

- A.-M. Legendre, *Essai sur la théorie des nombres*, 2nd ed. (1808), pp. 405–406.
- M. Chamberland, A. Straub, Weakening the Legendre conjecture, arXiv:2602.22502 (2026). https://arxiv.org/abs/2602.22502
- J. Sorenson, J. Webster, An algorithm and computation to verify Legendre's conjecture up to 7 · 10¹³, Research in Number Theory 11 (2025), doi:10.1007/s40993-024-00589-4; arXiv:2401.13753.
- J. Sorenson, J. Webster, Strong pseudoprimes to twelve prime bases, Math. Comp. 86 (2017) 985–1003; arXiv:1509.00864.
- T. Oliveira e Silva, S. Herzog, S. Pardi, Empirical verification of the even Goldbach conjecture and computation of prime gaps up to 4 · 10¹⁸, Math. Comp. 83 (2014) 2033–2060.
- T. R. Nicely, maximal prime gaps (gap 1476 after 1425172824437699411): https://faculty.lynchburg.edu/~nicely/gaps/gaps.html
- A. W. Dudek, D. R. Johnston, Almost primes between all squares, J. Number Theory 278 (2026) 726–745; arXiv:2501.18048.
- P. J. Campbell, On the existence of integers with at most 3 prime factors between every pair of consecutive squares, arXiv:2603.10356 (2026).
- D. R. Johnston, J. P. Sorenson, S. N. Thomas, J. E. Webster, Primes and almost primes between cubes, arXiv:2601.15564 (2026).
- M. Cully-Hugill, Primes between consecutive powers, J. Number Theory 247 (2023) 100–117; arXiv:2107.14468. Explicit estimates for the distribution of primes, PhD thesis, UNSW Canberra (2023).
- M. Cully-Hugill, D. R. Johnston, On the error term in the explicit formula of Riemann–von Mangoldt, Int. J. Number Theory 19 (2023) 1205–1228; part II, Funct. Approx. Comment. Math. 73 (2025) 223–242, arXiv:2402.04272.
- A. W. Dudek, An explicit result for primes between cubes, Funct. Approx. Comment. Math. 55 (2016) 177–197; arXiv:1401.4233.
- C. Mattner, Prime numbers in short intervals, BSc (Honours) thesis, ANU (2017).
- G. Hoheisel, Primzahlprobleme in der Analysis, Sitzungsber. Preuss. Akad. Wiss. 33 (1930) 3–11.
- A. E. Ingham, On the difference between consecutive primes, Quart. J. Math. 8 (1937) 255–266.
- R. C. Baker, G. Harman, The difference between consecutive primes, Proc. LMS 72 (1996) 261–280.
- R. C. Baker, G. Harman, J. Pintz, The difference between consecutive primes, II, Proc. LMS 83 (2001) 532–562.
- R. Li, arXiv:2308.04458 (claimed exponent 0.52).
- H. Cramér, On the order of magnitude of the difference between consecutive prime numbers, Acta Arith. 2 (1936) 23–46.
- E. Carneiro, M. B. Milinovich, K. Soundararajan, Fourier optimization and prime gaps, Comment. Math. Helv. 94 (2019) 533–568.
- D. Bazzanella, Primes between consecutive squares, Arch. Math. 75 (2000) 29–34.
- C. Jia, Almost all short intervals containing prime numbers, Acta Arith. 76 (1996) 21–84.
- L. Guth, J. Maynard, New large value estimates for Dirichlet polynomials, arXiv:2405.20552; Annals of Mathematics (to appear, 2026).
- J.-R. Chen, On the distribution of almost primes in an interval, Sci. Sinica 18 (1975) 611–627.
- W. H. Mills, A prime-representing function, Bull. AMS 53 (1947) 604.
- I. Niven, Functions which represent prime numbers, Proc. AMS 2 (1951) 753–755.
- C. K. Caldwell, Y. Cheng, Determining Mills' constant and a note on Honaker's problem, J. Integer Seq. 8 (2005) 05.4.1.
- K. Matomäki, Prime-representing functions, Acta Math. Hungar. 128 (2010) 307–314.
- OEIS A014085, A144140, A003617, A014233: https://oeis.org/
- Wikipedia, Legendre's conjecture: https://en.wikipedia.org/wiki/Legendre%27s_conjecture

## Licence

Text: CC BY 4.0. Code: MIT (see LICENSE). Cite as: MV2 (Munim, Inc.), "Legendre's conjecture in
October 2026", 2026.
