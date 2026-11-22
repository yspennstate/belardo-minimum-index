"""Numerical checks for the paper on connected graphs of fixed order and size with minimum index.

    python check_min_index.py            # quadratic-form identity, length inequality, closed forms, constructed minimizers
    geng -c N E:E | python check_min_index.py exhaustive N S    # E = N + S, graph6 from nauty's geng

Requires Python 3.8 or later and numpy. Floating point throughout; the proofs do not depend on it.
"""
import itertools
import math
import random
import sys

import numpy as np


def split(n, s):
    """n + s = 3 s a + r with 0 <= r < 3 s."""
    return divmod(n + s, 3 * s)


def phi(n, s, t):
    a, r = split(n, s)
    return (3 * s - r) * math.tanh(a * t / 2) + r * math.tanh((a + 1) * t / 2)


def tstar(n, s):
    lo, hi = 1e-12, 60.0
    for _ in range(300):
        mid = (lo + hi) / 2
        if phi(n, s, mid) - s / math.tanh(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def lam(n, s):
    return 2 * math.cosh(tstar(n, s))


def beh(n, s):
    """Minimum over the weighted relaxation of Ben-Shimon, Eisner and Hoory (their Theorems 1.4 and 1.5)."""
    e = n + s
    d_lo, n_hi = divmod(2 * e, n)
    if n_hi == 0:
        return float(d_lo)
    n_lo, d_hi = n - n_hi, d_lo + 1
    if n_lo * d_lo <= n_hi * d_hi:
        d1, n1, d2, n2 = d_lo, n_lo, d_hi, n_hi
    else:
        d1, n1, d2, n2 = d_hi, n_hi, d_lo, n_lo
    x = n1 / n2
    return 0.5 * (d2 - d1 * x + math.sqrt(4 * d1 * d1 * x + (d2 - d1 * x) ** 2))


def subdivide(kernel_edges, lengths):
    """Replace kernel edge (u, v) (u == v for a loop) by a path with the given number of edges."""
    nk = 1 + max(max(e) for e in kernel_edges)
    edges, nxt = [], nk
    paths = []
    for (u, v), l in zip(kernel_edges, lengths):
        seq = [u] + list(range(nxt, nxt + l - 1)) + [v]
        nxt += l - 1
        paths.append(seq)
        edges += list(zip(seq, seq[1:]))
    A = np.zeros((nxt, nxt))
    for x, y in edges:
        if x == y or A[x, y]:
            return None, None, None  # not simple
        A[x, y] = A[y, x] = 1
    return A, paths, nk


def trial_vector(nv, nk, paths, t):
    f = np.zeros(nv)
    f[:nk] = 1.0
    for seq in paths:
        l = len(seq) - 1
        for k in range(1, l):
            f[seq[k]] = math.cosh((l / 2 - k) * t) / math.cosh(l * t / 2)
    return f


def random_cubic_multigraph(m, rng):
    """Configuration model on m vertices of degree 3 (m even); retried until connected."""
    while True:
        stubs = [v for v in range(m) for _ in range(3)]
        rng.shuffle(stubs)
        E = [(stubs[2 * i], stubs[2 * i + 1]) for i in range(len(stubs) // 2)]
        seen, todo = {0}, [0]
        while todo:
            x = todo.pop()
            for u, v in E:
                for p, q in ((u, v), (v, u)):
                    if p == x and q not in seen:
                        seen.add(q)
                        todo.append(q)
        if len(seen) == m:
            return E


def perfect_matchings(E, m):
    idx = [i for i, (u, v) in enumerate(E) if u != v]
    for M in itertools.combinations(idx, m // 2):
        cov = sorted(x for i in M for x in E[i])
        if cov == list(range(m)):
            yield set(M)


def check_identity(rng, trials=400):
    """Quadratic-form identity of the test-vector lemma: <f,(A - 2cosh t)f> = 2 s cosh t - 2 sinh t sum tanh(l_i t / 2), any kernel, any lengths."""
    worst = 0.0
    done = 0
    while done < trials:
        s = rng.randint(1, 5)
        E = random_cubic_multigraph(2 * s, rng)
        # also kernels with a vertex of degree > 3: merge two vertices sometimes
        if s >= 2 and rng.random() < 0.3:
            E = [(0 if u == 1 else u, 0 if v == 1 else v) for u, v in E]
            E = [(u - (u > 1), v - (v > 1)) for u, v in E]
        lengths = [rng.randint(1, 6) for _ in E]
        A, paths, nk = subdivide(E, lengths)
        if A is None:
            continue
        t = rng.uniform(0.05, 3.0)
        f = trial_vector(len(A), nk, paths, t)
        lhs = f @ A @ f - 2 * math.cosh(t) * (f @ f)
        ex = len(E) - nk
        rhs = 2 * ex * math.cosh(t) - 2 * math.sinh(t) * sum(math.tanh(l * t / 2) for l in lengths)
        worst = max(worst, abs(lhs - rhs))
        done += 1
    return worst


def check_lengths(rng, trials=20000):
    """Length inequality (the lemma on path lengths) on random inputs: sum tanh(l_i t/2) <= Phi(t)."""
    worst = -1.0
    for _ in range(trials):
        s = rng.randint(1, 6)
        n = rng.randint(2 * s, 12 * s)
        E = n + s
        q = rng.randint(1, 3 * s)
        cuts = sorted(rng.sample(range(1, E), q - 1)) if q > 1 else []
        parts = [b - a for a, b in zip([0] + cuts, cuts + [E])]
        drop = rng.randint(0, max(0, min(parts) - 1))
        parts[0] -= drop
        t = rng.uniform(0.01, 5.0)
        worst = max(worst, sum(math.tanh(l * t / 2) for l in parts) - phi(n, s, t))
    return worst


def check_minimizers(rng):
    """Balanced subdivisions of random cubic kernels: spectral radius equals 2 cosh t*."""
    worst, count = 0.0, 0
    for s in range(1, 6):
        for ratio in range(2, 13):
            n = ratio * s
            if n < 4:
                continue
            a, r = split(n, s)
            j = r // s
            for _ in range(6):
                E = random_cubic_multigraph(2 * s, rng)
                if j == 0:
                    lengths = [a] * len(E)
                else:
                    Ms = list(perfect_matchings(E, 2 * s))
                    if not Ms:
                        continue
                    M = rng.choice(Ms)
                    long_on_M = (j == 1)
                    lengths = [a + 1 if ((i in M) == long_on_M) else a for i in range(len(E))]
                A, paths, nk = subdivide(E, lengths)
                if A is None or len(A) != n:
                    continue
                rho = np.linalg.eigvalsh(A)[-1]
                worst = max(worst, abs(rho - lam(n, s)))
                count += 1
    return worst, count


def g6(line):
    b = [c - 63 for c in line]
    m = b[0]
    A = np.zeros((m, m))
    bits = [(x >> k) & 1 for x in b[1:] for k in range(5, -1, -1)]
    i = 0
    for j in range(1, m):
        for k in range(j):
            if bits[i]:
                A[k, j] = A[j, k] = 1
            i += 1
    return A


def structure(A, a):
    """Degree set, number of branch vertices, path lengths, set of long incidences per branch vertex."""
    deg = A.sum(1).astype(int)
    adj = [list(np.nonzero(A[v])[0]) for v in range(len(deg))]
    br = [v for v in range(len(deg)) if deg[v] >= 3]
    used, out = set(), []
    for v in br:
        for w in adj[v]:
            if (v, w) in used:
                continue
            prev, cur, L = v, w, 1
            used.add((v, w))
            while deg[cur] == 2:
                nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
                prev, cur, L = cur, nxt, L + 1
            used.add((cur, prev))
            out.append((v, cur, L))
    longs = {v: 0 for v in br}
    for v, w, L in out:
        if L == a + 1:
            longs[v] += 1
            longs[w] += 1
    return sorted(set(deg.tolist())), len(br), sorted(L for _, _, L in out), sorted(set(longs.values()))


def exhaustive(n, s):
    a, r = split(n, s)
    bound = lam(n, s)
    count, mn, cand = 0, float('inf'), []
    for line in sys.stdin.buffer:
        line = line.strip()
        if not line:
            continue
        count += 1
        rho = np.linalg.eigvalsh(g6(line))[-1]
        if rho < mn - 1e-9:
            mn, cand = rho, [(rho, line)]
        elif rho < mn + 1e-9:
            cand.append((rho, line))
            mn = min(mn, rho)
    mins = [c for c in cand if c[0] < mn + 1e-9]
    print(f"n={n} s={s} graphs={count} a={a} r={r} s|n={n % s == 0} bound={bound:.15f} "
          f"min={mn:.15f} min-bound={mn - bound:.3e} minimizers={len(mins)}")
    for _, gl in mins:
        degs, nb, lens, le = structure(g6(gl), a)
        print(f"  {gl.decode()} degrees={degs} branch={nb} lengths={lens} long_ends={le}")


def verify(path):
    """Read the output of run_exhaustive.sh and check each row against the main theorem."""
    import re
    head = re.compile(r'n=(\d+) s=(\d+) graphs=(\d+) a=(\d+) r=(\d+) s\|n=(True|False) bound=(\S+) min=(\S+) '
                      r'min-bound=(\S+) minimizers=(\d+)$')
    body = re.compile(r'\s+\S+ degrees=\[([\d, ]+)\] branch=(\d+) lengths=\[([\d, ]*)\] long_ends=\[([\d, ]*)\]$')
    rows, cur, failures = [], None, 0
    for line in open(path, encoding='ascii'):
        line = line.rstrip('\n')
        m = head.match(line)
        if m:
            n, s, g, a, r = (int(m.group(i)) for i in range(1, 6))
            cur = dict(n=n, s=s, a=a, r=r, div=m.group(6) == 'True', gap=float(m.group(9)), k=int(m.group(10)), seen=0)
            rows.append(cur)
            continue
        m = body.match(line)
        if m and cur is not None:
            cur['seen'] += 1
            if not cur['div']:
                continue
            degs = [int(x) for x in m.group(1).split(',')]
            lens = [int(x) for x in m.group(3).split(',') if x.strip()]
            ends = [int(x) for x in m.group(4).split(',') if x.strip()]
            s, a, r = cur['s'], cur['a'], cur['r']
            good = degs == ([3] if cur['n'] == 2 * s else [2, 3]) and int(m.group(2)) == 2 * s
            good = good and len(lens) == 3 * s and set(lens) <= {a, a + 1} and ends == [r // s]
            if not good:
                failures += 1
                print('  not as predicted:', cur['n'], cur['s'], line.strip())
    for c in rows:
        if c['seen'] != c['k']:
            failures += 1
            print(f"  n={c['n']} s={c['s']}: {c['seen']} minimizer lines for {c['k']} minimizers")
        if c['div'] and abs(c['gap']) > 1e-14 or not c['div'] and c['gap'] <= 0:
            failures += 1
            print(f"  n={c['n']} s={c['s']}: min - bound = {c['gap']:.3e}")
    checked = sum(c['seen'] for c in rows if c['div'])
    print(f'rows {len(rows)}, minimizers checked against the theorem {checked}, failures {failures}')


def main():
    if len(sys.argv) == 4 and sys.argv[1] == 'exhaustive':
        exhaustive(int(sys.argv[2]), int(sys.argv[3]))
        return
    if len(sys.argv) == 3 and sys.argv[1] == 'verify':
        verify(sys.argv[2])
        return
    rng = random.Random(20261121)
    print('quadratic-form identity, largest deviation:', '%.2e' % check_identity(rng))
    print('length inequality, largest value of sum - Phi (must be <= 0 up to rounding):', '%.2e' % check_lengths(rng))
    w, c = check_minimizers(rng)
    print('balanced subdivisions checked:', c, ' largest |rho - 2cosh t*|:', '%.2e' % w)
    closed = {2: 3.0, 3: 1 + math.sqrt(3), 4: (1 + math.sqrt(17)) / 2, 5: math.sqrt(6), 8: (1 + math.sqrt(13)) / 2}
    print('ratio n/s, 2cosh t*, relaxation of Ben-Shimon-Eisner-Hoory, closed form')
    for ratio in range(2, 13):
        L, B = lam(ratio, 1) if ratio >= 4 else lam(2 * ratio, 2), beh(6 * ratio, 6)
        cf = closed.get(ratio)
        print(f'  {ratio:2d}  {L:.12f}  {B:.12f}  ' + (f'{cf:.12f}' if cf else '-'))
    for n, s in [(5, 2), (7, 2), (7, 3), (9, 4), (11, 2), (11, 3)]:
        print(f'  n={n} s={s}: 2cosh t* = {lam(n, s):.12f}')
    agree, larger, smaller = 0, 0, 0
    for s in range(1, 13):
        for n in range(2 * s, 12 * s + 1):
            d = lam(n, s) - beh(n, s)
            if abs(d) < 1e-9:
                agree += 1
                if not (n <= 5 * s or n == 8 * s):
                    print(f'  agreement outside n <= 5s, n = 8s: n={n} s={s}')
            elif d > 0:
                larger += 1
                if n <= 5 * s or n == 8 * s:
                    print(f'  2cosh t* larger at n={n} s={s} although n <= 5s or n = 8s')
            else:
                smaller += 1
                print(f'  2cosh t* smaller than the relaxation at n={n} s={s}')
    print(f'pairs s <= 12, 2s <= n <= 12s: agree to 1e-9 {agree}, 2cosh t* larger {larger}, smaller {smaller}')


if __name__ == '__main__':
    main()
