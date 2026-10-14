import io
import os
import sys
import unittest
from contextlib import redirect_stdout
from decimal import Decimal, localcontext

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import legendre as L  # noqa: E402

# OEIS A014085, n = 0..80
A014085 = [
    0, 2, 2, 2, 3, 2, 4, 3, 4, 3,
    5, 4, 5, 5, 4, 6, 7, 5, 6, 6,
    7, 7, 7, 6, 9, 8, 7, 8, 9, 8,
    8, 10, 9, 10, 9, 10, 9, 9, 12, 11,
    12, 11, 9, 12, 11, 13, 10, 13, 15, 10,
    11, 15, 16, 12, 13, 11, 12, 17, 13, 16,
    16, 13, 17, 15, 14, 16, 15, 15, 17, 13,
    21, 15, 15, 17, 17, 18, 22, 14, 18, 23, 13,
]

# OEIS A003617: smallest d-digit prime, d = 1..20
A003617 = [
    2, 11, 101, 1009, 10007, 100003, 1000003, 10000019, 100000007, 1000000007,
    10000000019, 100000000003, 1000000000039, 10000000000037, 100000000000031,
    1000000000000037, 10000000000000061, 100000000000000003, 1000000000000000003,
    10000000000000000051,
]

# OEIS A144140 (and Chamberland-Straub 2026, Example 4.1)
A144140 = [10, 20, 24, 27, 32, 65, 121, 139, 141, 187, 306, 321, 348, 1006, 1051]

# Chamberland-Straub 2026, Example 5.3 (power 2) and equation (1) / Example 5.2 (power 3)
MILLS_2 = [2, 5, 29, 853, 727613, 529420677791]
A_2 = Decimal("1.524699960538094359923363575688")
MILLS_3 = [2, 11, 1361, 2521008887]
A_3 = Decimal("1.306377883863080690468614492602")

PSI_12 = 318665857834031151167461


class TestPrimality(unittest.TestCase):
    def test_flags_small(self):
        f = L.prime_flags(30)
        self.assertEqual([k for k in range(31) if f[k]],
                         [2, 3, 5, 7, 11, 13, 17, 19, 23, 29])
        self.assertEqual(list(L.prime_flags(0)), [0])
        self.assertEqual(list(L.prime_flags(1)), [0, 0])

    def test_is_prime_matches_sieve(self):
        f = L.prime_flags(20000)
        for n in range(-3, 20001):
            self.assertEqual(L.is_prime(n), n >= 0 and bool(f[n]), n)

    def test_strong_pseudoprimes(self):
        self.assertTrue(L.is_strong_probable_prime(2047, 2))
        self.assertEqual(23 * 89, 2047)
        self.assertFalse(L.is_prime(2047))

        self.assertEqual(151 * 751 * 28351, 3215031751)
        for a in (2, 3, 5, 7):
            self.assertTrue(L.is_strong_probable_prime(3215031751, a))
        self.assertFalse(L.is_prime(3215031751))

        psi_11 = 3825123056546413051
        self.assertEqual(149491 * 747451 * 34233211, psi_11)
        for a in L.FIRST_PRIMES[:11]:
            self.assertTrue(L.is_strong_probable_prime(psi_11, a))
        self.assertFalse(L.is_prime(psi_11))

        for a in L.FIRST_PRIMES[:12]:
            self.assertTrue(L.is_strong_probable_prime(PSI_12, a))
        self.assertFalse(L.is_prime(PSI_12))

        for a in L.FIRST_PRIMES:
            self.assertTrue(L.is_strong_probable_prime(L.PSI_13, a))
        with self.assertRaises(ValueError):
            L.is_prime(L.PSI_13)

    def test_smallest_d_digit_primes(self):
        for d, p in enumerate(A003617, 1):
            self.assertEqual(L.first_prime_after(10 ** (d - 1) - 1), p, d)


class TestCounts(unittest.TestCase):
    def test_A014085(self):
        self.assertEqual(L.primes_between_powers(80), A014085)

    def test_partial_sums_are_pi(self):
        a = L.primes_between_powers(999)
        self.assertEqual(sum(a[:100]), 1229)    # pi(10^4)
        self.assertEqual(sum(a), 78498)         # pi(10^6)
        self.assertGreaterEqual(min(a[1:]), 2)

    def test_cubes_small(self):
        self.assertEqual(L.primes_between_powers(4, 3), [0, 4, 5, 9, 12])

    def test_oppermann_to_2000(self):
        counts = L.oppermann_counts(2000)
        a = L.primes_between_powers(2000)
        self.assertEqual(counts[1], (0, 1))
        self.assertEqual(counts[2], (1, 1))
        self.assertEqual(counts[3], (1, 1))
        self.assertEqual(counts[4], (2, 1))
        for n in range(2, 2001):
            lower, upper = counts[n]
            self.assertGreaterEqual(min(lower, upper), 1, n)
            self.assertEqual(lower + upper, a[n], n)

    def test_witnesses_match_sieve(self):
        f = L.prime_flags(302 ** 2)

        def next_prime(m):
            k = m + 1
            while not f[k]:
                k += 1
            return k

        for n in range(2, 301):
            p, q = L.oppermann_witnesses(n)
            self.assertEqual(p, next_prime(n * n), n)
            self.assertEqual(q, next_prime(n * (n + 1)), n)
        self.assertEqual(L.oppermann_witnesses(1), (None, 3))

    def test_witness_large(self):
        n = 10 ** 12
        p, q = L.oppermann_witnesses(n)
        self.assertIsNotNone(p)
        self.assertIsNotNone(q)
        self.assertTrue(n * n < p < n * (n + 1) < q < (n + 1) ** 2)
        for s in range(2, 1000):
            self.assertNotEqual(p % s, 0)
            self.assertNotEqual(q % s, 0)


class TestGapLemma(unittest.TestCase):
    def test_oliveira_e_silva_table(self):
        self.assertEqual(L.legendre_range_from_gaps(4 * 10 ** 18, 1476), (738, 2000000000))
        self.assertEqual(L.legendre_range_from_gaps(100, 5), (3, 10))
        self.assertEqual(L.legendre_range_from_gaps(100, 1), (2, 10))

    def test_lemma_on_a_small_table(self):
        x = 10 ** 5
        f = L.prime_flags(x + 1000)
        primes = [k for k in range(x + 1001) if f[k]]
        gap = max(primes[i + 1] - primes[i] for i in range(len(primes) - 1) if primes[i] <= x)
        lo, hi = L.legendre_range_from_gaps(x, gap)
        a = L.primes_between_powers(hi)
        for n in range(lo, hi + 1):
            self.assertGreaterEqual(a[n], 1, n)


class TestMills(unittest.TestCase):
    def test_squares_sequence(self):
        q = L.mills_sequence(6)
        self.assertEqual(q, MILLS_2)
        self.assertEqual(L.excesses(q), [1, 4, 12, 4, 22])
        self.assertTrue(L.conjecture_5_4_holds(q))

    def test_squares_constant(self):
        q = L.mills_sequence(7)
        self.assertEqual(q[:6], MILLS_2)
        self.assertTrue(q[5] ** 2 < q[6] < (q[5] + 1) ** 2)
        self.assertTrue(L.conjecture_5_4_holds(q))
        lo, hi = L.mills_bounds(q)
        lo6, hi6 = L.mills_bounds(q[:6])
        with localcontext() as ctx:
            ctx.prec = 60
            tol = Decimal("1e-29")
            self.assertLess(hi - lo, Decimal("1e-25"))
            self.assertTrue(lo6 <= lo < hi <= hi6)
            self.assertLessEqual(lo - tol, A_2)
            self.assertLessEqual(A_2, hi + tol)

    def test_floor_of_powers(self):
        with localcontext() as ctx:
            ctx.prec = 60
            for k in range(1, 7):
                self.assertEqual(int(A_2 ** (2 ** k)), MILLS_2[k - 1], k)
            for k in range(1, 4):
                self.assertEqual(int(A_3 ** (3 ** k)), MILLS_3[k - 1], k)

    def test_cubes(self):
        q = L.mills_sequence(4, 3)
        self.assertEqual(q, MILLS_3)
        self.assertEqual(L.excesses(q, 3), [3, 30, 6])
        lo, hi = L.mills_bounds(q, 3)
        with localcontext() as ctx:
            ctx.prec = 60
            tol = Decimal("1e-29")
            self.assertLess(hi - lo, Decimal("6.4e-12"))
            self.assertLessEqual(lo - tol, A_3)
            self.assertLessEqual(A_3, hi + tol)
        with self.assertRaises(ValueError):
            L.mills_sequence(5, 3)


class TestThreeHalves(unittest.TestCase):
    def test_A144140(self):
        self.assertEqual(L.three_halves_gaps(2000), A144140)


class TestCli(unittest.TestCase):
    def run_cli(self, *argv):
        buf = io.StringIO()
        with redirect_stdout(buf):
            code = L.main(list(argv))
        return code, buf.getvalue()

    def test_commands(self):
        code, out = self.run_cli("count", "10")
        self.assertEqual(code, 0)
        self.assertIn("10 5", out)
        code, out = self.run_cli("oppermann", "500")
        self.assertIn("holds for 2 <= n <= 500", out)
        code, out = self.run_cli("gaps", "4000000000000000000", "1476")
        self.assertIn("738 <= n <= 2000000000", out)
        code, out = self.run_cli("mills", "6")
        self.assertIn("q_6 = 529420677791", out)
        code, out = self.run_cli("witness", "10")
        self.assertIn("101", out)
        code, out = self.run_cli("three-halves", "100")
        self.assertIn("[10, 20, 24, 27, 32, 65]", out)
        code, out = self.run_cli()
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
