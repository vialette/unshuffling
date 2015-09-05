import itertools as it

def binary_strings_with_composition(p, q):
    """Generate all binary strings with p 0's and q 1's.
    """
    for pos0 in it.combinations(range(p + q), p):
        yield "".join(['0' if pos in pos0 else '1' for pos in range(p + q)])

def binary_strings(n):
    """Generate all binary strings of length n.
    """
    for p in range(n+1):
        for s in binary_strings_with_composition(p, n-p):
            yield s

def is_shuffle_of(u, v, w, cache=set()):
    """Return true if w is in the shuffle of u with v.
    """
    if (u, v) in cache:
        return False

    if len(u) + len(v) != len(w):
        return False

    if not u or not v or not w:
        if u + v == w:
            return True
        else:
            return False

    if u[0] != w[0] and v[0] != w[0]:
        return False

    if u[0] == w[0] and is_shuffle_of(u[1:], v, w[1:], cache):
            return True

    if v[0] == w[0] and is_shuffle_of(u, v[1:], w[1:], cache):
            return True

    cache.add((u, v))
    return False

def shuffle_square_roots(u):
    """Return all square roots of u.
    """
    n = len(u)

    if n % 2:
        # even length, no square roots
        return []

    # dynalic programming
    T = {}
    T[(1, 0)] = set([u[0]])
    for i in range(2, n + 1):
        for j in range(max(0, i - n//2), i//2 + 1):
            T[(i, j)] = set()
            if j > 0:
                for v in T[(i - 1, j - 1)]:
                    if v[j - 1] == u[i - 1]:
                        T[(i, j)].add(v)
            if j <= min((i-1)/2, n/2):
                for v in T[(i - 1, j)]:
                    T[(i, j)].add(v + u[i-1])

    # done, all square roots are in T[(n, n//2)])
    return (v for v in T[(n, n//2)])

def is_square(u):
    """Return True iff u is a square.
    """
    return shuffle_square_roots(u) != []

def count_squares(n):
    """Return the number of binary squares of length n.
    """
    # no odd length binary squares
    if n % 2 == 1:
        return 0

    # enumerate and test
    return len([u
                for p in range(0, n+1, 2)
                for u in binary_strings_with_composition(p, n-p)
                if is_square(u)])

def count_subsequences(u, v):
    """Return the number of occurrences of v in u as a subsequence.
    """
    return sum([1 for w in it.combinations(u, len(v)) if v == "".join([c for c in w])])

def count_subsequences_0110_1001(u):
    """Return the number of occurrences of 0110 or 1001 in u as a subsequence.
    """
    return count_subsequences(u, "0110") + count_subsequences(u, "1001")

def fact(n):
    """Return factorial n.
    """
    res = 1;
    for i in range(2, n + 1):
        res *= i
    return res

def binomial(n ,k):
    """Return n choose k.
    """
    result = 1
    for i in range(1, k+1):
        result = result * (n-i+1) / i
    return result

def lower_bound(n, k):
    if k <= (n - k):
        return (2*n - 2*k - 1) * binomial(k + 1, 2)
    else:
        return (2*k - 1) * binomial(n - k + 1, 2)

if __name__ == '__main__':
    n = 6
    h = {}

    total = 0
    sq_total = 0
    non_sq_total = 0

    for k in range(1, n//2 + 1):
        lb = lower_bound(n, k)
        for s in binary_strings_with_composition(2*k, 2*n - 2*k):
            square_roots = [r for r in shuffle_square_roots(s)]
            l = count_subsequences_0110_1001(s)
            total += 1
            if len(square_roots) > 0:
                sq_total += 1
            else:
                non_sq_total += 1
            if l < lb:
                if len(square_roots) == 0:
                    print("%s is not a square but has lower bound %d" % (s, l))
                    import sys
                    sys.exit(1)
                if l not in h:
                    h[l] = [s]
                else:
                    h[l].append(s)

    for l in sorted(h.keys()):
        print("%d (%d): %s" % (l, len(h[l]), [s for s in h[l]]))
    print
    print("total: %d" % (total,))
    print("squares/non squares: %d/%d" % (sq_total, non_sq_total))
    ll = sum(len(h[key]) for key in h.keys())
    print("squares below bound: %d (ratio %f)" % (ll, (1.0 * ll) / sq_total))
    print
    print

    h = {}
    for k in range(1, n/2 + 1):
        min_value = None
        min_s     = None
        max_value = None
        max_s     = None
        squares = []

        for s in binary_strings_with_composition(2*k, 2*n - 2*k):
            square_roots = [r for r in shuffle_square_roots(s)]
            squares.append((s, square_roots, count_subsequences_0110_1001(s), lower_bound(n, k)))
            if count_subsequences_0110_1001(s) < lower_bound(n, k):
                if len(square_roots) == 0:
                    print("%s. count_subsequences_0110_1001 = %d, bound = %d" % (s, count_subsequences_0110_1001(s), lower_bound(n, k)))
            if len(square_roots) == 0:
                if min_value is None:
                    min_value = count_subsequences_0110_1001(s)
                    min_s = [s]
                elif min_value > count_subsequences_0110_1001(s):
                    min_value = count_subsequences_0110_1001(s)
                    min_s = [s]
                elif min_value == count_subsequences_0110_1001(s):
                    min_s.append(s)
                if max_value is None:
                    max_value = count_subsequences_0110_1001(s)
                    max_s = [s]
                elif max_value < count_subsequences_0110_1001(s):
                    max_value = count_subsequences_0110_1001(s)
                    max_s = [s]
                elif max_value == count_subsequences_0110_1001(s):
                    max_s.append(s)
        print("n=%d, k=%d" % (n, k))
        print("min_value / bound = %d / %d" % (min_value, lower_bound(n, k)))
        print(min_s)
        print("max_value / bound = %d / %d" % (max_value, lower_bound(n, k)))
        print(max_s)
        print
