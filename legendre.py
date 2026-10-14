"""Exact-integer companion to "Legendre's conjecture in October 2026".

Every decision is made with Python integers (Decimal only for the Mills-type constants), so
nothing near a boundary depends on floating point. Standard library only.

    python legendre.py count 100              # a(n) = #primes in (n^2, (n+1)^2), OEIS A014085
    python legendre.py oppermann 10000        # both halves (n^2, n(n+1)), (n(n+1), (n+1)^2) by sieve
    python legendre.py witness 1000000000000  # least primes past n^2 and n(n+1), n up to 1.8e12
    python legendre.py gaps 4000000000000000000 1476   # n settled by a table of prime gaps
    python legendre.py mills 7                # q1 = 2, q_{k+1} = least prime > q_k^2, bounds on A
    python legendre.py mills 4 3              # the same with cubes (Mills, Caldwell-Cheng)
    python legendre.py three-halves 2000      # n with no prime in [n^1.5, (n+1)^1.5], OEIS A144140
"""

import sys
from decimal import Decimal, localcontext
from math import isqrt

FIRST_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)

# psi_13 (Sorenson and Webster, Math. Comp. 86 (2017) 985-1003): the smallest strong
# pseudoprime to all of the first 13 prime bases. Below it, those 13 bases decide primality.
PSI_13 = 3317044064679887385961981


def is_strong_probable_prime(n, a):
    """One Miller-Rabin round: True if the odd number n > 2 passes base a."""
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def is_prime(n):
    """Deterministic primality for n < PSI_13; raises ValueError at or above it."""
    if n < 2:
        return False
    for p in FIRST_PRIMES:
        if n % p == 0:
            return n == p
    if n >= PSI_13:
        raise ValueError("n >= psi_13: the 13-base Miller-Rabin test is not proved there")
    return all(is_strong_probable_prime(n, a) for a in FIRST_PRIMES)


def first_prime_after(m):
    """The least prime > m."""
    k = m + 1
    while not is_prime(k):
        k += 1
    return k


def prime_flags(limit):
    """bytearray f with f[k] == 1 exactly when k is prime, for 0 <= k <= limit."""
    flags = bytearray([1]) * (limit + 1)
    for k in range(min(2, limit + 1)):
        flags[k] = 0
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p::p] = bytes(len(range(p * p, limit + 1, p)))
    return flags


def primes_between_powers(n_max, k=2):
    """[a(0), ..., a(n_max)] with a(n) = number of primes strictly between n^k and (n+1)^k.

    k = 2 is OEIS A014085 (Legendre: a(n) >= 1 for n >= 1)."""
    flags = prime_flags((n_max + 1) ** k)
    return [sum(flags[n ** k + 1:(n + 1) ** k]) for n in range(n_max + 1)]


def oppermann_counts(n_max):
    """{n: (#primes in (n^2, n(n+1)), #primes in (n(n+1), (n+1)^2))} for 1 <= n <= n_max."""
    flags = prime_flags((n_max + 1) ** 2)
    out = {}
    for n in range(1, n_max + 1):
        mid = n * (n + 1)
        out[n] = (sum(flags[n * n + 1:mid]), sum(flags[mid + 1:(n + 1) ** 2]))
    return out


def oppermann_witnesses(n):
    """(p, q): the least prime above n^2 and the least prime above n(n+1).

    p is replaced by None if it is not below n(n+1), q by None if it is not below (n+1)^2;
    a None for n >= 2 would be a counterexample to Oppermann's conjecture."""
    mid = n * (n + 1)
    p = first_prime_after(n * n)
    q = first_prime_after(mid)
    return (p if p < mid else None, q if q < (n + 1) ** 2 else None)


def legendre_range_from_gaps(x_max, max_gap):
    """Lemma 1 of the README.

    If every prime p <= x_max is followed by the next prime within max_gap, then (n^2, (n+1)^2)
    contains a prime for every n with lo <= n <= hi; returns (lo, hi)."""
    return max(2, (max_gap + 1) // 2), isqrt(x_max)


def mills_sequence(terms, power=2):
    """q_1 = 2, q_{k+1} = the least prime > q_k^power (Chamberland-Straub 2026, Ex. 5.2, 5.3)."""
    q = [2]
    while len(q) < terms:
        q.append(first_prime_after(q[-1] ** power))
    return q


def excesses(q, power=2):
    """q_{k+1} - q_k^power for consecutive terms."""
    return [q[i + 1] - q[i] ** power for i in range(len(q) - 1)]


def conjecture_5_4_holds(q, power=2):
    """True if q_{k+1} < (q_k + 1)^power for every consecutive pair (Conjecture 5.4 for power 2)."""
    return all(q[i + 1] < (q[i] + 1) ** power for i in range(len(q) - 1))


def mills_bounds(q, power=2, digits=60):
    """Bracket from Lemma 5.1: q_k^(power^-k) <= A < (q_k + 1)^(power^-k), k = len(q)."""
    k = len(q)
    with localcontext() as ctx:
        ctx.prec = digits
        e = Decimal(1) / Decimal(power ** k)
        return Decimal(q[-1]) ** e, Decimal(q[-1] + 1) ** e


def three_halves_gaps(n_max):
    """n <= n_max such that no prime p satisfies n^(3/2) <= p <= (n+1)^(3/2) (OEIS A144140)."""
    flags = prime_flags(isqrt((n_max + 1) ** 3))
    out = []
    for n in range(1, n_max + 1):
        lo = isqrt(n ** 3)
        if lo * lo < n ** 3:
            lo += 1
        hi = isqrt((n + 1) ** 3)
        if not any(flags[lo:hi + 1]):
            out.append(n)
    return out


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    cmd, args = argv[0], [int(a) for a in argv[1:]]
    if cmd == "count":
        a = primes_between_powers(args[0])
        for n, c in enumerate(a):
            print(n, c)
        if args[0] >= 1:
            print("min over 1 <= n <= %d: %d" % (args[0], min(a[1:])))
    elif cmd == "oppermann":
        counts = oppermann_counts(args[0])
        bad = [n for n in range(2, args[0] + 1) if min(counts[n]) == 0]
        if bad:
            print("Oppermann fails at", bad)
        else:
            print("Oppermann (hence Legendre) holds for 2 <= n <= %d" % args[0])
    elif cmd == "witness":
        n = args[0]
        p, q = oppermann_witnesses(n)
        print("n =", n)
        print("least prime in (n^2, n(n+1)):     ", p)
        print("least prime in (n(n+1), (n+1)^2): ", q)
    elif cmd == "gaps":
        lo, hi = legendre_range_from_gaps(args[0], args[1])
        print("Legendre holds for %d <= n <= %d (check n < %d separately)" % (lo, hi, lo))
    elif cmd == "mills":
        power = args[1] if len(args) > 1 else 2
        q = mills_sequence(args[0], power)
        for i, v in enumerate(q, 1):
            print("q_%d = %d" % (i, v))
        print("q_{k+1} - q_k^%d:" % power, excesses(q, power))
        print("q_{k+1} < (q_k + 1)^%d for all listed terms:" % power, conjecture_5_4_holds(q, power))
        lo, hi = mills_bounds(q, power)
        print("A in [%s, %s)" % (lo, hi))
    elif cmd == "three-halves":
        print(three_halves_gaps(args[0]))
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
