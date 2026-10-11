"""Nonlinear 2D PDE via fixed-point residual proofs.

Solves u_xx + u_yy + u^2 = f(x,y) with zero Dirichlet BCs on [-1,1]^2,
using float Newton on the Chebyshev tau discretization.
"""

from fractions import Fraction
import math


def pde2d_nl_coefficients(f_str, deg, max_iter=50):
    """Solve u_xx + u_yy + u^2 = f via float Newton.

    Returns list of lists c[m][n] as floats, shape (N+1)x(N+1).
    """
    N = deg
    if N % 2 == 0:
        raise ValueError("deg must be odd")
    n_vars = (N + 1) * (N + 1)

    def idx(m, n):
        return m * (N + 1) + n

    # Chebyshev D2
    def cheb_diff_mat(K):
        D = [[0.0] * (K + 1) for _ in range(K + 1)]
        if K >= 1:
            D[K - 1][K] = 2.0 * K
        for nn in range(K - 2, -1, -1):
            cbar = 2.0 if nn == 0 else 1.0
            for col in range(K + 1):
                D[nn][col] += D[nn + 2][col] / cbar
            D[nn][nn + 1] += 2.0 * (nn + 1) / cbar
        return D
    D1 = cheb_diff_mat(N)
    D2 = [[0.0] * (N + 1) for _ in range(N + 1)]
    for i in range(N + 1):
        for j in range(N + 1):
            s = 0.0
            for k in range(N + 1):
                s += D1[i][k] * D1[k][j]
            D2[i][j] = s

    # f in Chebyshev basis (float)
    def x_pow_cheb(ex):
        if ex == 0:
            return [1.0]
        p = [1.0]
        for _ in range(ex):
            q = [0.0] * (len(p) + 1)
            for k, pk in enumerate(p):
                if pk == 0:
                    continue
                q[abs(1 - k)] += pk / 2
                q[k + 1] += pk / 2
            p = q
        return p
    import re
    f_dict = {}
    s = f_str.replace(" ", "").replace("-", "+-")
    for term in [t for t in s.split("+") if t]:
        sign = 1.0
        if term.startswith("-"):
            sign = -1.0
            term = term[1:]
        if not term:
            continue
        coeff = 1.0
        ex, ey = 0, 0
        for fac in term.split("*"):
            if not fac:
                continue
            m = re.match(r'^x\^(\d+)$', fac)
            if m:
                ex = int(m.group(1))
                continue
            m = re.match(r'^y\^(\d+)$', fac)
            if m:
                ey = int(m.group(1))
                continue
            if fac == "x":
                ex = 1
                continue
            if fac == "y":
                ey = 1
                continue
            coeff *= float(Fraction(fac))
        f_dict[(ex, ey)] = f_dict.get((ex, ey), 0.0) + sign * coeff
    f_cheb = [[0.0] * (N + 1) for _ in range(N + 1)]
    for (ex, ey), coeff in f_dict.items():
        px = x_pow_cheb(ex)
        py = x_pow_cheb(ey)
        for m in range(min(len(px), N + 1)):
            for n_ in range(min(len(py), N + 1)):
                f_cheb[m][n_] += coeff * px[m] * py[n_]

    def residual(c_flat):
        # c_flat: list of n_vars floats
        # Returns list of n_vars floats
        c = [[c_flat[idx(m, n)] for n in range(N + 1)] for m in range(N + 1)]
        F = [0.0] * n_vars
        # Tau: (u_xx + u_yy + u^2 - f)_{m,n} = 0 for m,n=0..N-2
        # These occupy rows 0..(N-1)^2-1
        row = 0
        for m in range(N - 1):
            for n_ in range(N - 1):
                acc = 0.0
                # u_xx + u_yy
                for a in range(N + 1):
                    acc += D2[m][a] * c[a][n_]
                for b in range(N + 1):
                    acc += D2[n_][b] * c[m][b]
                # u^2: (u^2)_{m,n}
                u2 = 0.0
                for a in range(N + 1):
                    for b in range(N + 1):
                        for cc in range(N + 1):
                            for d in range(N + 1):
                                kx = 0.0
                                if abs(a - cc) == m:
                                    kx += 0.5
                                if a + cc == m:
                                    kx += 0.5
                                if kx == 0:
                                    continue
                                ky = 0.0
                                if abs(b - d) == n_:
                                    ky += 0.5
                                if b + d == n_:
                                    ky += 0.5
                                if ky == 0:
                                    continue
                                u2 += kx * ky * c[a][b] * c[cc][d]
                acc += u2
                acc -= f_cheb[m][n_]
                F[row] = acc
                row += 1
        for n_ in range(N + 1):
            s = sum(((-1) ** mm) * c[mm][n_] for mm in range(N + 1))
            F[row] = s
            row += 1
        for n_ in range(N + 1):
            s = sum(c[mm][n_] for mm in range(N + 1))
            F[row] = s
            row += 1
        for mm in range(1, N):
            s = sum(((-1) ** n_) * c[mm][n_] for n_ in range(N + 1))
            F[row] = s
            row += 1
        for mm in range(1, N):
            s = sum(c[mm][n_] for n_ in range(N + 1))
            F[row] = s
            row += 1
        return F

    # Initial guess: zeros (Newton will converge from here for small nonlinearity)
    # Note: do NOT use Galerkin linear solution as initial guess;
    # this solver uses raw Chebyshev basis (matching _build_pde),
    # while pde2d_poisson_coefficients uses Galerkin basis.
    c_flat = [0.0] * n_vars
    # Newton
    for it in range(max_iter):
        F = residual(c_flat)
        norm = math.sqrt(sum(x * x for x in F))
        if norm < 1e-10:
            break
        # Numerical Jacobian
        J = [[0.0] * n_vars for _ in range(n_vars)]
        eps = 1e-7
        for j in range(n_vars):
            c2 = c_flat[:]
            c2[j] += eps
            F2 = residual(c2)
            for i in range(n_vars):
                J[i][j] = (F2[i] - F[i]) / eps
        # Solve J*delta = -F
        M = [J[i][:] + [-F[i]] for i in range(n_vars)]
        for col in range(n_vars):
            piv = max(range(col, n_vars), key=lambda r: abs(M[r][col]))
            if abs(M[piv][col]) < 1e-12:
                raise ValueError("singular Jacobian in Newton")
            M[col], M[piv] = M[piv], M[col]
            pivval = M[col][col]
            for c2 in range(col, n_vars + 1):
                M[col][c2] /= pivval
            for r in range(n_vars):
                if r != col:
                    f = M[r][col]
                    for c2 in range(col, n_vars + 1):
                        M[r][c2] -= f * M[col][c2]
        delta = [M[i][n_vars] for i in range(n_vars)]
        c_flat = [c_flat[i] + delta[i] for i in range(n_vars)]
    else:
        raise ValueError("Newton did not converge for nonlinear PDE")
    # Reshape
    c = [[c_flat[idx(m, n)] for n in range(N + 1)] for m in range(N + 1)]
    return c
