"""Check of Theorem 6 against the computed minimizers, for s even and n an odd multiple of s/2 (n+s = 3sa+r, r/s = j+1/2).

Builds every graph B_L with L split (Definition 5: sets X, Y of s vertices, a perfect matching of short edges in X and
of long edges in Y, every vertex an end of two X-Y edges of which exactly j are long), keeps the simple ones, and compares
their indices and isomorphism classes with the minimizers listed by check_min_index.py. It also prints
kappa(n,s) = 2cosh(tau), tau the zero of the larger eigenvalue of the matrix Z(t) of the paper.
Usage: python split_family.py exhaustive_output.txt half_output.txt   (requires numpy and networkx)
"""
import itertools
import sys

import networkx as nx
import numpy as np


def perfect_matchings(vs):
    if not vs:
        yield []
        return
    v = vs[0]
    for i in range(1, len(vs)):
        rest = vs[1:i] + vs[i + 1:]
        for m in perfect_matchings(rest):
            yield [(v, vs[i])] + m


def crossing_multigraphs(s):
    # s x s nonnegative integer matrices with all row and column sums 2
    rows = [r for r in itertools.product(range(3), repeat=s) if sum(r) == 2]
    out = []

    def rec(prefix, colsum):
        if len(prefix) == s:
            if all(c == 2 for c in colsum):
                out.append(prefix)
            return
        for r in rows:
            cs = [c + x for c, x in zip(colsum, r)]
            if max(cs) <= 2:
                rec(prefix + [r], cs)

    rec([], [0] * s)
    return out


def subdivide(edges, lengths, nk):
    g = nx.Graph()
    g.add_nodes_from(range(nk))
    nxt = nk
    for (u, v), l in zip(edges, lengths):
        path = [u] + list(range(nxt, nxt + l - 1)) + [v]
        nxt += l - 1
        for x, y in zip(path, path[1:]):
            if x == y or g.has_edge(x, y):
                return None
            g.add_edge(x, y)
    return g


def kappa(n, s):
    a, r = divmod(n + s, 3 * s)
    j = (2 * r // s - 1) // 2

    def zeta(t):
        u = lambda l: np.tanh(l * t / 2)
        sg = lambda l: 1 / np.sinh(l * t)
        eta = u(a + 1) - u(a)
        b = 1 / np.tanh(t) - 3 * u(a) - (j + 0.5) * eta
        om = (2 - j) * sg(a) + j * sg(a + 1)
        return b - om + np.sqrt(om * om + eta * eta / 4)

    lo, hi = 1e-6, np.log((3 + np.sqrt(5)) / 2)
    for _ in range(200):
        mid = (lo + hi) / 2
        if zeta(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 2 * np.cosh((lo + hi) / 2)


def rho(g):
    return max(np.linalg.eigvalsh(nx.to_numpy_array(g)))


def family(n, s):
    E = n + s
    a, r = divmod(E, 3 * s)
    assert 2 * r % s == 0 and (2 * r // s) % 2 == 1, (n, s, a, r)
    j = (2 * r // s - 1) // 2
    A = list(range(s))
    B = list(range(s, 2 * s))
    graphs = []
    for MA in perfect_matchings(A):
        for MB in perfect_matchings(B):
            for X in crossing_multigraphs(s):
                cross = [(A[i], B[k]) for i in range(s) for k in range(s) for _ in range(X[i][k])]
                for lens in itertools.product((a, a + 1), repeat=len(cross)):
                    ok = True
                    for v in A + B:
                        if sum(1 for (e, l) in zip(cross, lens) if v in e and l == a + 1) != j:
                            ok = False
                            break
                    if not ok:
                        continue
                    edges = MA + MB + cross
                    L = [a] * len(MA) + [a + 1] * len(MB) + list(lens)
                    assert sum(L) == E and len(edges) == 3 * s
                    g = subdivide(edges, L, 2 * s)
                    if g is not None and nx.is_connected(g):
                        assert g.number_of_nodes() == n and g.number_of_edges() == E
                        graphs.append(g)
    classes = []
    for g in graphs:
        if not any(nx.is_isomorphic(g, h) for h in classes):
            classes.append(g)
    return a, r, j, classes


def exhaustive_rows(path):
    rows = {}
    cur = None
    for line in open(path, encoding="utf-8"):
        if line.startswith("n="):
            f = dict(x.split("=", 1) for x in line.split() if "=" in x and not x.startswith("s|n"))
            cur = (int(f["n"]), int(f["s"]), int(f.get("graphs", 0)))
            mn = f["min"].replace("np.float64(", "").rstrip(")")
            rows[cur] = {"min": float(mn), "count": int(f["minimizers"]), "g6": []}
        elif line.startswith("  ") and cur is not None and not line.strip().startswith("degrees"):
            rows[cur]["g6"].append(line.split()[0])
        elif not line.startswith(" "):
            cur = None
    return rows


def main():
    rows = {}
    for path in sys.argv[1:]:
        rows.update(exhaustive_rows(path))
    checked = 0
    bad = 0
    for (n, s, ng) in sorted(rows):
        E = n + s
        a, r = divmod(E, 3 * s)
        if s % 2 or (2 * r) % s or (2 * r // s) % 2 == 0 or a < 1:
            continue
        a, r, j, classes = family(n, s)
        row = rows[(n, s, ng)]
        mins = [nx.from_graph6_bytes(x.encode()) for x in row["g6"]]
        rv = sorted(rho(g) for g in classes)
        same_value = rv and abs(rv[0] - row["min"]) < 1e-9 and abs(rv[-1] - row["min"]) < 1e-9
        matched = all(any(nx.is_isomorphic(g, h) for h in mins) for g in classes)
        covered = all(any(nx.is_isomorphic(g, h) for h in classes) for g in mins)
        kp = kappa(n, s)
        ok = (same_value and matched and covered and len(classes) == row["count"] == len(mins)
              and abs(kp - row["min"]) < 1e-12)
        checked += 1
        bad += not ok
        print(f"n={n} s={s} graphs={ng} a={a} r={r} j={j} family={len(classes)} computed_minimizers={row['count']} "
              f"family_rho=[{rv[0] if rv else None:.12f},{rv[-1] if rv else None:.12f}] computed_min={row['min']:.12f} "
              f"kappa={kp:.12f} "
              f"{'OK' if ok else 'MISMATCH'}")
    print(f"rows checked {checked}, mismatches {bad}")


if __name__ == "__main__":
    main()
