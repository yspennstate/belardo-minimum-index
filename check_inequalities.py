"""High-precision check of the inequalities used in the proof of Theorem 6 (odd multiples of s/2).

For j in {0,1,2} and a >= 1, tau is the zero of the larger eigenvalue of Z(t), with
Z = [[b+eta/2-omega, omega], [omega, b-eta/2-omega]], u_l = tanh(lt/2), sigma_l = 1/sinh(lt),
eta = u_{a+1}-u_a, b = coth t - 3u_a - (j+1/2)eta, omega = (2-j)sigma_a + j sigma_{a+1}.
At t = tau the script checks: e^tau in (sqrt 2, (3+sqrt 5)/2); R > 1 and R - 1/R = eta/omega; 0 < -b < eta^2/(8 omega);
gamma_l (R-1) < 2 eta_l for l >= 1; R^2 < eta_a/eta_{a+1}; R^2 < eta_{a-1}/eta_a for a >= 2; and
m_7 = min(eta_a, eta_{a-1}-eta_a, eta_a-eta_{a+1}) > s eta^2/(16 sigma_{a+1}) for s in {2,4,6}.
The working precision is 60+a decimal digits, since eta is of order e^{-a tau}.
Usage: python check_inequalities.py   (requires mpmath)
"""
from mpmath import mp, mpf, tanh, sinh, cosh, coth, sqrt, exp, log


def u(l, t):
    return mpf(0) if l == 0 else tanh(l * t / 2)


def sigma(l, t):
    return 1 / sinh(l * t)


def eta_l(l, t):
    # u_{l+1} - u_l in closed form, without cancellation at large l
    return sinh(t / 2) / (cosh((l + 1) * t / 2) * cosh(l * t / 2))


def gamma_l(l, t):
    return sinh(t) / (sinh(l * t) * sinh((l + 1) * t))


def parts(a, j, t):
    eta = eta_l(a, t)
    b = coth(t) - 3 * u(a, t) - (j + mpf(1) / 2) * eta
    omega = (2 - j) * sigma(a, t) + j * sigma(a + 1, t)
    return eta, b, omega


def zeta(a, j, t):
    # larger eigenvalue of Z(t)
    eta, b, omega = parts(a, j, t)
    return b - omega + sqrt(eta * eta / 4 + omega * omega)


def root(a, j):
    lo, hi = mpf("1e-6"), log((3 + sqrt(5)) / 2)
    assert zeta(a, j, lo) > 0 and zeta(a, j, hi) < 0, (a, j)
    for _ in range(int(mp.prec) + 20):
        mid = (lo + hi) / 2
        if zeta(a, j, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


n_checks = 0
fails = []
for j in (0, 1, 2):
    for a in list(range(1, 60)) + [100, 300, 1000]:
        mp.dps = 60 + a
        t = root(a, j)
        eta, b, omega = parts(a, j, t)
        R = (omega + eta / 2 - b) / omega
        conds = {
            "e^tau range": sqrt(2) < exp(t) < (3 + sqrt(5)) / 2,
            "R > 1": R > 1,
            "R - 1/R = eta/omega": abs(R - 1 / R - eta / omega) < mpf(10) ** (-40),
            "0 < -b < eta^2/(8 omega)": 0 < -b < eta * eta / (8 * omega),
            "gamma_l (R-1) < 2 eta_l": all(gamma_l(l, t) * (R - 1) < 2 * eta_l(l, t) for l in range(1, a + 5)),
            "R^2 < eta_a/eta_(a+1)": R * R < eta_l(a, t) / eta_l(a + 1, t),
            "R^2 < eta_(a-1)/eta_a": a < 2 or R * R < eta_l(a - 1, t) / eta_l(a, t),
        }
        m7 = min(eta_l(a, t), eta_l(a - 1, t) - eta_l(a, t), eta_l(a, t) - eta_l(a + 1, t))
        for s in (2, 4, 6):
            conds[f"m_7 > s eta^2/(16 sigma_(a+1)), s={s}"] = m7 > s * eta * eta / (16 * sigma(a + 1, t))
        for k, v in conds.items():
            n_checks += 1
            if not v:
                fails.append((j, a, k))
        if a in (1, 2, 10, 1000):
            print(f"j={j} a={a} kappa=2cosh(tau)={float(2 * cosh(t)):.12f} R={float(R):.6f} "
                  f"m_7/(6 eta^2/(16 sigma_(a+1)))={float(m7 / (6 * eta * eta / (16 * sigma(a + 1, t)))):.4f}")
print("checks", n_checks, "failures", len(fails), fails[:10])
# control: the last inequality fails for larger s, so s <= 6 is used
mp.dps = 60
t = root(1, 0)
eta, b, omega = parts(1, 0, t)
m7 = min(eta_l(1, t), eta_l(0, t) - eta_l(1, t), eta_l(1, t) - eta_l(2, t))
print("control: smallest even s with m_7 <= s eta^2/(16 sigma_2) at a=1, j=0:",
      next(s for s in range(2, 200, 2) if not m7 > s * eta * eta / (16 * sigma(2, t))))
