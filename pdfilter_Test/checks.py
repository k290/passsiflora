#!/usr/bin/env python3
# checks.py -- exact-integer checks of the statements about resultants, discriminants and the two
# Euclidean chains that are proved in pdfilter_chains.tex (Appendix B lists what is checked).
# Independent of chain_trace.cpp: own polynomial arithmetic, own models of chain_block and
# chain_general, resultants as Sylvester determinants (a quoted standard result), checked against
# the root definition for split polynomials.  Python 3 standard library only, fixed seed.
#
# Usage:   python3 checks.py [quick] [--gp FILE]
#   quick      a tenth of the cases (a few seconds; the full run takes about a minute)
#   --gp FILE  also write a PARI/GP script that compares 300 of the discriminants with poldisc():
#              gp -q FILE < /dev/null   prints "gp: disc samples checked" and nothing else if all agree
# Exit status 1 if any check fails.  Section C8 prints the numbers of the worked example (q = 41);
# section C9 checks Section 9 of the document (zeroing the x^8 coefficient) with its own random generator.
import random, sys
from math import comb

RNG = random.Random(20261024)
QUICK = "quick" in sys.argv[1:]
GP_FILE = sys.argv[sys.argv.index("--gp") + 1] if "--gp" in sys.argv[1:-1] else None
FAIL = []

def check(cond, what):
    if not cond:
        FAIL.append(what)
        if len(FAIL) < 20:
            print("FAILED:", what)

# ---------- polynomials: ascending coefficient lists, trimmed (zero polynomial = [])
def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p
def red(p, q): return trim([c % q for c in p])
def deg(p): return len(p) - 1
def sym(x, q):
    x %= q
    return x - q if x > (q - 1) // 2 else x
def chi(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1
def inv(a, q): return pow(a % q, q - 2, q)
def powm(a, e, q):  # e may be negative (a != 0)
    return pow(a % q, e, q) if e >= 0 else pow(inv(a, q), -e, q)
def add(a, b, q):
    n = max(len(a), len(b))
    return red([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)], q)
def sub(a, b, q):
    n = max(len(a), len(b))
    return red([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)], q)
def scal(a, c, q): return red([c * x for x in a], q)
def shift(a, k): return [0] * k + list(a) if a else []
def mul(a, b, q):
    if not a or not b: return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return red(r, q)
def mul_int(a, b):
    if not a or not b: return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return trim(r)
def pdivmod(a, b, q):
    a = red(a, q); b = red(b, q)
    assert b
    quo = [0] * max(1, len(a) - len(b) + 1)
    il = inv(b[-1], q)
    while a and deg(a) >= deg(b):
        c = a[-1] * il % q
        k = deg(a) - deg(b)
        quo[k] = c
        a = sub(a, shift(scal(b, c, q), k), q)
    return red(quo, q), a
def pmod(a, b, q): return pdivmod(a, b, q)[1]
def deriv_int(f): return trim([k * f[k] for k in range(1, len(f))])
def evalp(p, x, q):
    r = 0
    for c in reversed(p):
        r = (r * x + c) % q
    return r
def translate_int(f, t):  # f(x+t), exact integers
    n = deg(f)
    g = [0] * (n + 1)
    for k, a in enumerate(f):
        for j in range(k + 1):
            g[j] += a * comb(k, j) * t ** (k - j)
    return trim(g)
def translate_mod(f, t, q): return red(translate_int([c % q for c in f], t % q), q)
def from_roots_int(roots, lead=1):
    p = [lead]
    for r in roots:
        p = mul_int(p, [-r, 1])
    return p
def show(p, q):  # descending symmetric residues, for printing
    return [sym(c, q) for c in reversed(p)] if p else [0]

# ---------- determinants and resultants
def det_mod(M, q):
    M = [[x % q for x in row] for row in M]
    n = len(M); d = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c]), None)
        if piv is None: return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]; d = -d
        d = d * M[c][c] % q
        ic = inv(M[c][c], q)
        for r in range(c + 1, n):
            fct = M[r][c] * ic % q
            if fct:
                for k in range(c, n):
                    M[r][k] = (M[r][k] - fct * M[c][k]) % q
    return d % q
def det_int(M):  # Bareiss, exact
    M = [row[:] for row in M]
    n = len(M)
    if n == 0: return 1
    sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((r for r in range(k + 1, n) if M[r][k]), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]
def sylvester(A, B):
    m, n = deg(A), deg(B); N = m + n
    rows = []
    for i in range(n):
        row = [0] * N
        for j in range(m + 1): row[i + j] = A[m - j]
        rows.append(row)
    for i in range(m):
        row = [0] * N
        for j in range(n + 1): row[i + j] = B[n - j]
        rows.append(row)
    return rows
def res(V, P, q):
    """Res(V,P) = lc(V)^deg(P) * prod_{V(b)=0} P(b), via det Syl(V,P) mod q."""
    V = red(V, q); P = red(P, q)
    assert V
    if not P:
        assert deg(V) >= 1
        return 0
    if deg(V) + deg(P) == 0: return 1
    return det_mod(sylvester(V, P), q)
def res_int(V, P):
    if deg(V) + deg(P) == 0: return 1
    return det_int(sylvester(V, P))
def res_rec(V, P, q):
    """The recursion of the lemma "The resultant is computed inside F" (Section 3)."""
    V = red(V, q); P = red(P, q)
    if deg(V) == 0: return powm(V[0], deg(P), q)
    if deg(P) == 0: return powm(P[0], deg(V), q)
    if deg(P) < deg(V):
        return (-1) ** (deg(V) * deg(P)) * res_rec(P, V, q) % q
    R = pmod(P, V, q)
    if not R: return 0
    return powm(V[-1], deg(P) - deg(R), q) * res_rec(V, R, q) % q

def rand_poly(d, q, lead=None):
    p = [RNG.randrange(q) for _ in range(d)] + [lead if lead is not None else RNG.randrange(1, q)]
    return red(p, q)
def rand_pair(q):
    """U, V nonzero with deg U >= deg V >= 1, including degenerate constructions."""
    kind = RNG.randrange(6)
    dv = RNG.randrange(1, 7)
    V = rand_poly(dv, q)
    if kind == 0:
        U = rand_poly(RNG.randrange(dv, 12), q)
    elif kind == 1:   # U = Q V + small R: big degree drops
        Qp = rand_poly(RNG.randrange(0, 5), q)
        U = add(mul(Qp, V, q), rand_poly(RNG.randrange(0, max(1, dv - 2)), q), q)
    elif kind == 2:   # U a multiple of V: zero remainder
        U = mul(rand_poly(RNG.randrange(0, 5), q), V, q)
    elif kind == 3:   # common factor
        G = rand_poly(1, q)
        V = mul(G, rand_poly(max(0, dv - 1), q), q)
        U = mul(G, rand_poly(RNG.randrange(max(0, deg(V) - 1), 9), q), q)
    elif kind == 4:   # repeated roots in V
        G = rand_poly(1, q)
        V = mul(mul(G, G, q), rand_poly(RNG.randrange(0, 3), q), q)
        U = rand_poly(RNG.randrange(deg(V), 10), q)
    else:             # sparse
        U = red([RNG.choice([0, 0, 0, RNG.randrange(q)]) for _ in range(RNG.randrange(dv, 11))] + [RNG.randrange(1, q)], q)
        V = red([RNG.choice([0, 0, RNG.randrange(q)]) for _ in range(dv)] + [RNG.randrange(1, q)], q)
    if not U or deg(U) < deg(V) or deg(V) < 1:
        return rand_pair(q)
    return U, V

PR1 = [13, 17, 29, 37, 41, 53, 73, 89, 97, 113, 761, 1009]   # 1 mod 4
PR3 = [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]         # 3 mod 4
N_REP = 150 if QUICK else 1500

# ---------- C0: Sylvester determinant against the root definition (split V)
cnt = 0
for q in PR1 + PR3:
    for _ in range(N_REP // 5):
        roots = [RNG.randrange(q) for _ in range(RNG.randrange(1, 7))]
        if RNG.random() < 0.3: roots.append(roots[0])
        c = RNG.randrange(1, q)
        V = red(from_roots_int(roots, c), q)
        P = rand_poly(RNG.randrange(0, 8), q)
        if not P: continue
        rootdef = powm(c, deg(P), q)
        for b in roots: rootdef = rootdef * evalp(P, b, q) % q
        check(res(V, P, q) == rootdef, f"C0 root definition q={q} V={V} P={P}")
        check(res_rec(V, P, q) == rootdef, f"C0 recursion q={q}")
        cnt += 1
print(f"C0 Res = det Syl = root definition (split V) and = the recursion of Section 3: {cnt} cases")

# ---------- C1: swap, constant cases, P mod V, multiplicativity, L1
cnt = cntz = 0
for q in PR1 + PR3:
    for _ in range(N_REP):
        U, V = rand_pair(q)
        u, v = deg(U), deg(V)
        check(res(U, V, q) == (-1) ** (u * v) * res(V, U, q) % q, "C1 swap")
        cc = RNG.randrange(1, q)
        check(res([cc], U, q) == powm(cc, u, q), "C1 constant V")
        check(res(V, [cc], q) == powm(cc, v, q), "C1 constant P")
        Rr = pmod(U, V, q)
        if Rr:
            check(res(V, U, q) == powm(V[-1], u - deg(Rr), q) * res(V, Rr, q) % q, "C1 P mod V")
        else:
            check(res(V, U, q) == 0, "C1 P mod V zero")
        P2 = rand_poly(RNG.randrange(0, 4), q)
        if P2:
            check(res(V, mul(U, P2, q), q) == res(V, U, q) * res(V, P2, q) % q, "C1 multiplicative")
        lu, lv = U[-1], V[-1]
        W = sub(scal(U, lv, q), shift(scal(V, lu, q), u - v), q)
        check(deg(W) <= u - 1, "L1 degree drops")
        if W:
            e = deg(W) - u + v
            check(res(V, W, q) == powm(lv, e, q) * res(V, U, q) % q, f"L1 q={q} U={U} V={V}")
            cnt += 1
        else:
            check(res(V, U, q) == 0, "L1 zero W")
            cntz += 1
print(f"C1 swap / constant / P mod V / multiplicativity / L1: {cnt} steps with W != 0, {cntz} with W = 0")

# ---------- C2: discriminant: roots definition vs (-1)^(n(n-1)/2) Res(f,f'); reduction mod q
cnt = 0
for _ in range(N_REP):
    n = RNG.randrange(1, 9)
    roots = [RNG.randrange(-9, 10) for _ in range(n)]
    f = from_roots_int(roots)
    D = 1
    for i in range(n):
        for j in range(i + 1, n):
            D *= (roots[i] - roots[j]) ** 2
    check(D == (-1) ** (n * (n - 1) // 2) * res_int(f, deriv_int(f)), f"C2 disc formula roots={roots}")
    check(D == (-1) ** (n * (n - 1) // 2) * res_int(deriv_int(f), f), "C2 disc formula with Res(f',f)")
    cnt += 1
disc_samples = []
for _ in range(N_REP):
    n = RNG.randrange(2, 11)
    f = [RNG.randrange(-30, 31) for _ in range(n)] + [1]
    D = (-1) ** (n * (n - 1) // 2) * res_int(f, deriv_int(f))
    disc_samples.append((f, D))
    for q in PR1 + PR3:
        if n % q == 0: continue
        check(D % q == (-1) ** (n * (n - 1) // 2) * res(red(f, q), red(deriv_int(f), q), q) % q, "C2 reduction mod q")
print(f"C2 disc = prod (r_i-r_j)^2 = (-1)^(n(n-1)/2) Res(f,f') on {cnt} split integer polynomials; reduction mod q on {len(disc_samples)} x {len(PR1+PR3)}")
if GP_FILE:
    with open(GP_FILE, "w") as fh:
        for f, D in disc_samples[:300]:
            fh.write("if(poldisc(Pol(%s)) != %d, print(\"MISMATCH %s\"));\n" % (str(list(reversed(f))), D, str(f)))
        fh.write("print(\"gp: disc samples checked\");\n")
    print("  wrote", GP_FILE)

# ---------- C3: L2, one chain_block step = lc^2 * (A mod B)
def block_step(A, B, m, q):
    """One step of chain_block on coefficient arrays (index k = coefficient of x^k)."""
    lc = B[m]; top = A[m + 1]
    R = [0] * (m + 1)
    R[0] = lc * A[0] % q
    for k in range(1, m + 1):
        R[k] = (lc * A[k] - top * B[k - 1]) % q
    t = R[m]
    R2 = [(lc * R[k] - t * B[k]) % q for k in range(m)]
    return R, R2
def arr(p, length):
    return [(p[k] if k < len(p) else 0) for k in range(length)]
cnt = 0
for q in PR1 + PR3:
    for _ in range(N_REP):
        m = RNG.randrange(1, 10)
        B = rand_poly(m, q)                       # degree exactly m
        A = red([RNG.randrange(q) for _ in range(m + 2)], q)   # degree <= m+1 (any)
        if RNG.random() < 0.3:                    # degenerate: A close to a multiple of B
            A = add(mul(rand_poly(1, q), B, q), rand_poly(RNG.randrange(0, m), q), q)
            A = red(A[:m + 2], q)
        R, R2 = block_step(arr(A, m + 2), arr(B, m + 1), m, q)
        lc = B[-1]
        check(red(R2, q) == scal(pmod(A, B, q), lc * lc, q), f"L2 q={q} A={A} B={B}")
        cnt += 1
print(f"C3 L2 (R2 = lc^2 * (A mod B)): {cnt} random steps")

# ---------- C4: Theorem A
def chain_block(f, q):
    """chain_block for monic f of degree n (ascending integer coefficients), as the Fortran with 9 -> n-1."""
    n = deg(f)
    A = [c % q for c in f]
    B = [(k * A[k]) % q for k in range(1, n + 1)] + [0]
    ok = True
    P = [A[:], B[:]]
    checks = []
    for m in range(n - 1, 0, -1):
        R, R2 = block_step(A, B, m, q)
        R2 = R2 + [0] * (n + 1 - len(R2))
        if m >= 2:
            checks.append(R2[m - 1])
            ok = ok and R2[m - 1] != 0
        A, B = B, R2
        P.append(R2[:])
    base = B[0]
    acc = sym(pow(base, (q - 1) // 2, q), q)
    return ok, base, acc, (not (ok and acc == -1)), P, checks
def remainders(f, q):
    r = [red(f, q), red(deriv_int(f), q)]
    while r[-1] and deg(r[-1]) > 0:
        r.append(pmod(r[-2], r[-1], q))
    return r
def gen_f(n, kind, lo=-20, hi=20):
    f = [RNG.randint(lo, hi) for _ in range(n)] + [1]
    if kind == 1: f[n - 1] = 0
    if kind == 2: f[n - 1] = 0; f[n - 2] = 0
    if kind == 3:
        for k in range(1, n, 2): f[k] = 0
    if kind == 4:   # repeated root: disc = 0
        g = [RNG.randint(-5, 5) for _ in range(n - 2)] + [1]
        r = RNG.randint(-5, 5)
        f = mul_int(g, [r * r, -2 * r, 1])
    if kind == 5:   # shift of a polynomial with zero x^(n-1), x^(n-2)
        f[n - 1] = 0; f[n - 2] = 0
        f = translate_int(f, RNG.randint(-6, 6))
    return f
cntA = {"1mod4": 0, "3mod4": 0}; normal = abnormal = 0
for q in PR1 + PR3:
    for rep in range(N_REP):
        n = RNG.choice([2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 12])
        if n % q == 0: continue
        f = gen_f(n, rep % 6)
        D = (-1) ** (n * (n - 1) // 2) * res_int(f, deriv_int(f))
        ok, base, acc, keep, P, checks = chain_block(f, q)
        r = remainders(f, q)
        isnormal = all(s < len(r) and r[s] and deg(r[s]) == n - s for s in range(n))
        check(ok == isnormal, f"ThA ok <=> normal, q={q} f={f}")
        if not ok:
            abnormal += 1
            continue
        normal += 1
        sig = [1, 1]
        for s in range(1, n):
            lam = P[s][n - s] % q
            sig.append(lam * lam * sig[s - 1] % q)
        rn = r[n] if len(r) > n else []
        for s in range(n + 1):
            rs = r[s] if s < len(r) else []
            check(red(P[s], q) == scal(rs, sig[s], q), f"ThA P_s = sigma_s r_s, q={q} f={f} s={s}")
            check(sig[s] != 0 and chi(sig[s], q) == 1, "ThA sigma nonzero square")
        prod = 1
        for s in range(1, n): prod = prod * r[s][-1] % q
        rn0 = rn[0] if rn else 0
        check(res(r[1], r[0], q) == prod * prod * rn0 % q, "ThA Res(f',f) = (prod lc)^2 r_n")
        check(chi(base, q) == chi((-1) ** (n * (n - 1) // 2), q) * chi(D, q), f"ThA sharp form q={q} f={f}")
        if q % 4 == 1:
            check(chi(base, q) == chi(D, q), "ThA")
            cntA["1mod4"] += 1
        else:
            cntA["3mod4"] += 1
print(f"C4 Theorem A: {normal} normal chains ({cntA}), {abnormal} abnormal; ok <=> normal checked on all")

# ---------- C5: Theorem B with the invariant after every iteration
def chain_general(f, q, trace=None):
    n = deg(f)
    U = red(f, q); V = red(deriv_int(f), q)
    Rff = res(V, U, q)        # Res(f', f)
    mm, eps, c = 1, 1, 1
    it = 0
    while True:
        it += 1
        du, dv = deg(U), deg(V)
        check(du >= dv >= 1, "ThB du >= dv >= 1")
        lu, lv = U[-1], V[-1]
        W = sub(scal(U, lv, q), shift(scal(V, lu, q), du - dv), q)
        check(deg(W) <= du - 1, "ThB step lowers degree")
        if not W:
            check(Rff == 0, "ThB zero W => Res(f',f) = 0")
            if trace is not None: trace.append(dict(it=it, U=U, V=V, W=W, end="zero"))
            return dict(resz=True, mm=mm, it=it, eps=eps)
        dw = deg(W); e = dw - du + dv
        if e % 2:
            mm = mm * lv % q; c = c * powm(lv, -(e + 1) // 2, q) % q
        else:
            c = c * powm(lv, -e // 2, q) % q
        swapped = dw < dv
        if swapped:
            eps *= (-1) ** (dv * dw)
            U, V = V, W
        else:
            U = W
        rec = dict(it=it, e=e, swapped=swapped, mm=mm, eps=eps, c=c, U=U, V=V, dw=dw)
        if deg(V) == 0:
            c0 = V[0]
            if deg(U) % 2: mm = mm * c0 % q
            c = c * powm(c0, deg(U) // 2, q) % q
            check(Rff == eps * mm * c * c % q, "ThB final identity")
            rec.update(mm=mm, c=c, end="const")
            if trace is not None: trace.append(rec)
            return dict(resz=False, mm=mm, it=it, eps=eps, c=c)
        rv = res(V, U, q)
        check(Rff == eps * mm * c * c * rv % q, f"ThB invariant q={q} f={f} it={it}")
        rec["ResVU"] = rv
        if trace is not None: trace.append(rec)
cnt = {"1mod4": 0, "3mod4": 0}; zeros = 0; maxit = {}
for q in PR1 + PR3:
    for rep in range(N_REP):
        n = RNG.choice([2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 12])
        if n % q == 0: continue
        f = gen_f(n, rep % 6)
        D = (-1) ** (n * (n - 1) // 2) * res_int(f, deriv_int(f))
        g = chain_general(f, q)
        check(g["it"] <= 2 * n - 2, f"ThB termination bound n={n} it={g['it']}")
        maxit[n] = max(maxit.get(n, 0), g["it"])
        if g["resz"]:
            check(D % q == 0, "ThB resz => disc = 0 mod q"); zeros += 1
            continue
        check(g["mm"] % q != 0, "ThB mm nonzero")
        check(chi(g["mm"], q) == chi((-1) ** (n * (n - 1) // 2) * g["eps"], q) * chi(D, q), "ThB sharp form")
        if q % 4 == 1:
            check(chi(g["mm"], q) == chi(D, q), f"ThB q={q} f={f}"); cnt["1mod4"] += 1
        else:
            cnt["3mod4"] += 1
        # (A first draft also claimed "normal sequence => exactly 2n-2 iterations"; that is false,
        #  e.g. f1 = x^4+3x^2+2x+8 takes 5: a_{n-1} = 0 makes the first elementary step drop by 2.
        #  Only the bound it <= 2n-2 above is claimed.)
print(f"C5 Theorem B: {cnt} chains ending in a constant, {zeros} ending in 0 (all with disc = 0 mod q); max iterations by n: {maxit}")

# ---------- C6: L3, shift invariance, depressed form
cnt = 0
for _ in range(N_REP * 4):
    n = RNG.randrange(3, 13)
    f = gen_f(n, RNG.randrange(6), -50, 50) if RNG.random() < 0.5 else [RNG.randint(-99, 99) for _ in range(n)] + [1]
    an1, an2 = f[n - 1], f[n - 2]
    A = f[:]; B = [k * f[k] for k in range(1, n + 1)] + [0]
    BIG = 10 ** 60   # exact integer arithmetic through the same code (no reduction happens)
    R, R2 = block_step(A, B, n - 1, BIG)
    piv = R2[n - 2] if R2[n - 2] < BIG // 2 else R2[n - 2] - BIG
    check(piv == 2 * n * an2 - (n - 1) * an1 * an1, f"L3 pivot formula n={n}")
    p1 = -an1; p2 = an1 * an1 - 2 * an2
    check(piv == p1 * p1 - n * p2, "L3 power sums")
    t = RNG.randint(-30, 30)
    g = translate_int(f, t)
    check(2 * n * g[n - 2] - (n - 1) * g[n - 1] ** 2 == piv, "L3 shift invariance")
    cnt += 1
for _ in range(N_REP):
    n = RNG.randrange(2, 10)
    roots = [RNG.randrange(-20, 21) for _ in range(n)]
    f = from_roots_int(roots)
    S = sum((roots[i] - roots[j]) ** 2 for i in range(n) for j in range(i + 1, n))
    check(2 * n * f[n - 2] - (n - 1) * f[n - 1] ** 2 == -S, "L3 = -sum (r_i - r_j)^2")
cntd = 0
for q in [x for x in PR1 + PR3 if x > 12]:
    for rep in range(N_REP // 3):
        n = RNG.randrange(3, 12)
        f = gen_f(n, rep % 6)
        t = (-f[n - 1] * inv(n, q)) % q
        g = translate_mod(f, t, q)
        g = arr(g, n + 1)
        check(g[n - 1] % q == 0, "depressed")
        ks = [k for k in range(n - 1) if g[k] % q]
        expect = max(ks) if ks else -1
        r = remainders(f, q)
        check(deg(r[2]) == expect, f"Prop depressed form q={q} f={f}")
        cntd += 1
print(f"C6 L3 (pivot formula, power sums, shift invariance) {cnt} + {N_REP} root cases; depressed-form proposition {cntd} cases")

# ---------- C7: filter primes and Theorem C
def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True
def usable(q):
    if q % 4 != 1 or not is_prime(q): return False
    return all(p % q != 0 and pow(p, (q - 1) // 2, q) == 1 for p in (2, 5))
small, q = [], 2 ** 10 - 1
while len(small) < 6 and q > 256:
    if usable(q): small.append(q)
    q -= 2
large, q = [], 2 ** 25 - 1
while len(large) < 24 and q > 2 ** 24:
    if usable(q): large.append(q)
    q -= 2
FP = small + large
print("filter primes: small", small, "large", large[:3], "...", large[-1], "count", len(FP), "all 1 or 9 mod 40:", all(p % 40 in (1, 9) for p in FP))

def filter_stage(f, q, usegen):
    if usegen:
        g = chain_general(f, q)
        return not (not g["resz"] and sym(pow(g["mm"], (q - 1) // 2, q), q) == -1), "general"
    ok, base, acc, keep, P, _ = chain_block(f, q)
    if ok: return keep, "block"
    g = chain_general(f, q)
    return not (not g["resz"] and sym(pow(g["mm"], (q - 1) // 2, q), q) == -1), "pass2"
def gen_poly_test():
    """genPoly() of test_pdfilter_linux.cpp, without the irreducibility test (not needed for the discriminant)."""
    c = [RNG.randrange(100) for _ in range(5)]
    if RNG.randrange(4) == 0: c[4] = 0
    if 55 <= c[0] < 70: c[0] = -1
    elif 70 <= c[0] < 85: c[0] = 2
    elif c[0] >= 85: c[0] = 5
    else: c[0] = c[0] * c[0]
    f = [c[0], 0, c[1], 0, c[2], 0, c[3], 0, c[4], 0, 1]
    if RNG.randrange(4) != 0:
        t = 1 + RNG.randrange(30)
        if RNG.randrange(2): t = -t
        f = translate_int(f, t)
    return f, c
def passes(D):
    if D == 0: return False
    D = abs(D)
    for p in (2, 5):
        while D % p == 0: D //= p
    s = int(D ** 0.5)
    while s * s > D: s -= 1
    while (s + 1) ** 2 <= D: s += 1
    return s * s == D
stagesC = {"block": 0, "pass2": 0, "general": 0}; rejected = 0; npass = 0
NC = 60 if QUICK else 400
for rep in range(NC):
    if rep % 2 == 0:
        f, c = gen_poly_test()
    else:
        f = gen_f(10, rep % 6, -99, 99)
    D = -res_int(f, deriv_int(f))   # n = 10: (-1)^45 = -1
    ps = passes(D)
    npass += ps
    for q in (FP if rep < 40 else small + large[:4]):
        for usegen in (False, True):
            keep, how = filter_stage(f, q, usegen)
            stagesC[how] += 1
            if not keep:
                rejected += 1
                check(chi(D, q) == -1, "ThC reject => chi = -1")
                check(not ps, "ThC never rejects a passing polynomial")
print(f"C7 Theorem C: {NC} degree-10 polynomials ({npass} passing), stage runs {stagesC}, {rejected} rejections, each with chi(disc) = -1")

# genPoly discriminant formula disc(g(x^2)) = (-4)^m g(0) disc(g)^2, and the test's polynomials pass
for _ in range(N_REP // 3):
    m = RNG.randrange(1, 7)
    g = [RNG.randint(-20, 20) for _ in range(m)] + [1]
    f = [0] * (2 * m + 1)
    for k, a in enumerate(g): f[2 * k] = a
    Df = (-1) ** ((2 * m) * (2 * m - 1) // 2) * res_int(f, deriv_int(f))
    Dg = (-1) ** (m * (m - 1) // 2) * res_int(g, deriv_int(g))
    check(Df == (-4) ** m * g[0] * Dg * Dg, f"disc(g(x^2)) formula g={g}")
for _ in range(N_REP // 3):
    f, c = gen_poly_test()
    D = -res_int(f, deriv_int(f))
    check(D == 0 or passes(D), "genPoly polynomials pass")
    check(2 * 10 * f[8] - 9 * f[9] ** 2 == 20 * c[4], "genPoly pivot = 20 c4")
print("C7b disc(g(x^2)) = (-4)^m g(0) disc(g)^2 and genPoly pivot = 20*a8 checked")

# NF_PDFDEBUG polynomial x^10 + c x^8 + 5x^6 + 7x^4 + 11x^2 + 4 at the 30 filter primes (model, not the Fortran)
for cc in (3, 0):
    f = [4, 0, 11, 0, 7, 0, 5, 0, cc, 0, 1]
    oks = [chain_block(f, q)[0] for q in FP]
    keeps = [filter_stage(f, q, False)[0] for q in FP]
    D = -res_int(f, deriv_int(f))
    print(f"  c = {cc}: model chain_block ok at {sum(oks)} of 30 filter primes; kept at all: {all(keeps)}; disc passes: {passes(D)}")

# ---------- C8: the worked example (q = 41)
q = 41
print("\n===== worked example, q = 41 =====")
for name, f in (("f1", [8, 2, 3, 0, 1]), ("f2", [8, 2, 0, 0, 1]), ("f3", [8, 0, 0, 0, 1])):
    n = 4
    D = (-1) ** (n * (n - 1) // 2) * res_int(f, deriv_int(f))
    print(f"{name}: disc = {D}, disc mod 41 = {D % 41} (sym {sym(D, 41)}), chi = {chi(D, 41)}, Res(f',f) mod 41 = {sym(res(red(deriv_int(f), q), red(f, q), q), q)}, pivot 8a2-3a3^2 = {8 * f[2] - 3 * f[3] ** 2}")
    r = remainders(f, q)
    print("   true remainders r_s (descending, symmetric):", [show(x, q) for x in r])
    ok, base, acc, keep, P, checks = chain_block(f, q)
    print("   chain_block P_s:", [show(red(x, q), q) for x in P], "ok", ok, "base", sym(base, q), "acc", acc)
    if ok:
        sig = [1, 1]
        for s in range(1, n):
            lam = P[s][n - s] % q
            sig.append(lam * lam * sig[s - 1] % q)
        print("   sigma_s:", [sym(x, q) for x in sig], " lambda_s:", [sym(P[s][n - s], q) for s in range(1, n)])
        for s in range(1, n):
            A, B = red(P[s - 1], q), red(P[s], q)
            mo = pmod(A, B, q)
            print(f"   step s={s} (m={n - s}): lc={sym(B[-1], q)}, A mod B = {show(mo, q)}, lc^2*(A mod B) = {show(scal(mo, B[-1] ** 2, q), q)}, R2 = {show(red(P[s + 1], q), q)}")
        prod = 1
        for s in range(1, n): prod = prod * r[s][-1] % q
        print(f"   prod lc(r_s) = {sym(prod, q)}, prod^2 = {sym(prod * prod, q)}, r_4 = {show(r[4], q)}, (prod)^2 r_4 = {sym(prod * prod * r[4][0], q)}; disc/base = {sym(D * inv(base, q), q)}")
    tr = []
    gres = chain_general(f, q, tr)
    print("   chain_general:", gres)
    for rec in tr:
        print("     ", {k: (show(v, q) if isinstance(v, list) else (sym(v, q) if k in ("mm", "c", "ResVU") else v)) for k, v in rec.items()})
# Euler by repeated squaring, exactly as lines 309-316
for b0 in (8, -14, 5):
    base, acc, e = b0 % q, 1, (q - 1) // 2
    steps = []
    while e > 0:
        bit = e & 1
        if bit: acc = acc * base % q
        steps.append((e, bit, sym(acc, q), sym(base * base % q, q)))
        base = base * base % q
        e >>= 1
    print(f"Euler {b0}^20 mod 41: rows (e, bit, acc after, base after) = {steps}, result {sym(acc, q)}")
print("squares mod 41:", sorted({x * x % 41 for x in range(1, 41)}))
print("-disc(f2) = 27*2^4 - 256*8^3 =", 27 * 2 ** 4 - 256 * 8 ** 3)

# ---------- C9: Section 9, zeroing the x^8 coefficient
# Its own random generator, so the counts of C0-C8 do not change; the same cases in quick and full runs.
# The draws follow the order used for the numbers quoted in Section 9 of pdfilter_chains.tex.
print("\n===== C9: zeroing the x^8 coefficient (Section 9) =====")
RNG = random.Random(8)        # rand_poly() and rand_pair() use this name
def x8_deriv(f, q): return red([k * f[k] for k in range(1, len(f))], q)
def x8_defect(f):              # n f - x f' = sum (n-k) a_k x^k  (Lemma "the Euler defect")
    n = deg(f); return trim([(n - k) * f[k] for k in range(n + 1)])
def x8_remseq(f, q):
    r = [red(f, q), x8_deriv(f, q)]
    while r[-1] and deg(r[-1]) > 0:
        r.append(pmod(r[-2], r[-1], q))
    return r
def x8_degseq(f, q): return tuple(deg(x) for x in x8_remseq(f, q))
def x8_neg(p): return [c if k % 2 == 0 else -c for k, c in enumerate(p)]     # p(-x)
def x8_gen(c4zero=None):
    """f = G(x^2) as genPoly() builds it, c4 = 0 if c4zero is True; no irreducibility test, no shift."""
    c = [RNG.randrange(100) for _ in range(5)]
    if c4zero is True: c[4] = 0
    if 55 <= c[0] < 70: c[0] = -1
    elif 70 <= c[0] < 85: c[0] = 2
    elif c[0] >= 85: c[0] = 5
    else: c[0] = c[0] * c[0]
    return [c[0], 0, c[1], 0, c[2], 0, c[3], 0, c[4], 0, 1], c

# C9.1 Euler defect; f mod f' for a centred f; deg r_2 and its x^(n-2) coefficient pi/n^2 for any f
cnt = 0
for q in (41, 761, 33554249):
    for _ in range(500):
        n = RNG.choice([10, 10, 10, 3, 4, 5, 6, 7, 8, 9, 12])
        if q <= n: continue
        f = [RNG.randint(-99, 99) for _ in range(n)] + [1]
        lhs = trim([n * a for a in f]); rhs = shift(deriv_int(f), 1)
        diff = trim([(lhs[k] if k < len(lhs) else 0) - (rhs[k] if k < len(rhs) else 0) for k in range(max(len(lhs), len(rhs)))])
        check(diff == x8_defect(f), "C9.1 Euler defect identity over Z")
        g = f[:]
        if RNG.random() < 0.5: g[n - 1] = 0          # centred
        r2 = pmod(red(g, q), x8_deriv(g, q), q)
        if g[n - 1] % q == 0:
            check(r2 == scal(x8_defect(g), inv(n, q), q), "C9.1 f mod f' = (1/n) Euler defect")
        t = (-g[n - 1] * inv(n, q)) % q
        dep = arr(translate_mod(g, t, q), n + 1)
        ks = [k for k in range(n - 1) if dep[k] % q]
        check(deg(r2) == (max(ks) if ks else -1), "C9.1 deg r_2 = max k with g_k != 0")
        piv = 2 * n * g[n - 2] - (n - 1) * g[n - 1] ** 2
        c = r2[n - 2] if len(r2) > n - 2 else 0
        check(c % q == piv * inv(n * n, q) % q, "C9.1 x^(n-2) coefficient of r_2 = pi/n^2")
        check(dep[n - 2] % q == piv * inv(2 * n, q) % q, "C9.1 g_(n-2) = pi/(2n)")
        cnt += 1
print("C9.1 Euler defect, f mod f', deg r_2, coefficient pi/n^2:", cnt, "cases")

# C9.2 parity: f = G(x^2) => r_k(-x) = (-1)^k r_k(x); every W of chain_general has the parity of U, so dw <= du - 2
cnt = 0
for q in (761, 41):
    for _ in range(500):
        f, c = x8_gen(c4zero=RNG.random() < 0.5)
        for k, rk in enumerate(x8_remseq(f, q)):
            check(red(x8_neg(rk), q) == scal(rk, (-1) ** k, q), "C9.2 parity of r_k")
            if rk: check(deg(rk) % 2 == k % 2, "C9.2 deg r_k = k mod 2")
        U = red(f, q); V = x8_deriv(f, q)
        while True:
            du, dv = deg(U), deg(V)
            W = sub(scal(U, V[-1], q), shift(scal(V, U[-1], q), du - dv), q)
            if not W: break
            check(red(x8_neg(W), q) == scal(W, (-1) ** du, q), "C9.2 parity of W")
            check(deg(W) <= du - 2, "C9.2 even f: dw <= du - 2")
            if deg(W) < dv: U, V = V, W
            else: U = W
            if deg(V) == 0: break
        cnt += 1
print("C9.2 parity of the remainders and of the elementary steps:", cnt, "even polynomials")

# C9.3 degree sequences of the generator's polynomials, and accidental failures of later checks (observations, not theorems)
for q in (761, 41):
    RNG.random()               # one draw, kept so that the counts quoted in Section 9 stay reproducible
    for label, c4z, shiftit in (("c4=0 unshifted", True, False), ("c4=0 shifted", True, True), ("c4 random, shifted", None, True)):
        rows = []; later = 0; firstok = 0
        for _ in range(3000):
            f, c = x8_gen(c4zero=c4z)
            if shiftit:
                t = 1 + RNG.randrange(30); t = -t if RNG.randrange(2) else t
                f = translate_int(f, t)
            ds = x8_degseq(f, q); rows.append(ds)
            if len(ds) > 2 and ds[2] == 8:
                firstok += 1
                if ds != tuple(range(10, -1, -1)): later += 1
        freq = {}
        for d in rows: freq[d] = freq.get(d, 0) + 1
        top = sorted(freq.items(), key=lambda kv: -kv[1])[:4]
        print(f"C9.3 q={q} {label}: " + "; ".join(f"{100*v/len(rows):.1f}% {','.join(map(str, k))}" for k, v in top)
              + (f"; first check passed {firstok}, later degree drop or zero remainder in {later} ({100*later/firstok:.2f}%; 7/q = {700/q:.2f}%)" if firstok else ""))

# C9.4 exponent of a k-fold pseudo-step: Res(V, l_v^k (U mod V)) = l_v^(k v - (u - w)) Res(V, U); first step E_2 = v + w - 1
cnt = 0
for q in (41, 43, 761):
    for _ in range(1500):
        U, V = rand_pair(q)
        R = pmod(U, V, q)
        if not R: continue
        k = RNG.randrange(0, 4); lv = V[-1]
        E = k * deg(V) - (deg(U) - deg(R))
        check(res(V, scal(R, pow(lv, k, q), q), q) == powm(lv, E, q) * res(V, U, q) % q, "C9.4 E_k formula")
        cnt += 1
cases = {}
for q in (41, 761):
    for _ in range(600):
        f = [RNG.randint(-30, 30) for _ in range(10)] + [1]
        f[9] = 0; f[8] = 0
        if RNG.random() < 0.5: f[7] = 0            # 10,9,6 if a6 != 0, otherwise 10,9,7
        A = red(f, q); B = x8_deriv(f, q)
        R2 = scal(pmod(A, B, q), 100, q)           # chain_block's first R_2 = lambda^2 (A mod B), lambda = 10
        if not R2: continue
        w = deg(R2); E2 = 2 * 9 - (10 - w)
        check(E2 == 9 + w - 1, "C9.4 E_2 = v + w - 1")
        check(res(B, R2, q) == powm(10, E2, q) * res(B, A, q) % q, "C9.4 first-step exponent")
        cases[w] = cases.get(w, 0) + 1
print("C9.4 E_k formula:", cnt, "random steps; first steps by w =", dict(sorted(cases.items())), "(E_2 = 8 + w: 14 for w = 6, 15 for w = 7)")

# C9.5 two chain_general iterations with the same V make one chain_block step
cnt2 = cnt1 = 0
for q in (41, 761):
    for _ in range(1500):
        m = RNG.randrange(1, 10)
        B = rand_poly(m, q); A = rand_poly(m + 1, q)
        lam, tau = B[-1], A[-1]
        R = sub(scal(A, lam, q), shift(scal(B, tau, q), 1), q)       # round 1 = chain_general iteration 1 (u = m+1, v = m)
        t = R[m] if len(R) > m else 0
        R2 = sub(scal(R, lam, q), scal(B, t, q), q)                  # round 2
        if t % q:
            check(deg(R) == m, "C9.5 t != 0: no swap")
            W2 = sub(scal(R, lam, q), scal(B, R[-1], q), q)            # iteration 2, u = v = m, l_u = t
            check(W2 == R2, "C9.5 iteration 2 = R_2"); cnt2 += 1
        else:
            check(R2 == scal(R, lam, q), "C9.5 t = 0: R_2 = lambda R"); cnt1 += 1
print(f"C9.5 two iterations = one chain_block step: {cnt2} cases with t != 0, {cnt1} with t = 0")

# C9.6 the elementary step is not translation-covariant; the characters of mm agree anyway (Theorem B)
def x8_mm(f, q):
    U = red(f, q); V = x8_deriv(f, q); mm = 1; it = 0
    while True:
        it += 1
        du, dv = deg(U), deg(V); lu, lv = U[-1], V[-1]
        W = sub(scal(U, lv, q), shift(scal(V, lu, q), du - dv), q)
        if not W: return None, it
        if (deg(W) - du + dv) % 2: mm = mm * lv % q
        if deg(W) < dv: U, V = V, W
        else: U = W
        if deg(V) == 0:
            if deg(U) % 2: mm = mm * V[0] % q
            return mm, it
same_mm = diff_it = n_ok = noncov = 0
for q in (41, 761):
    for _ in range(400):
        f, c = x8_gen(c4zero=True)
        tt = 1 + RNG.randrange(30)
        g = translate_int(f, tt)
        mf, itf = x8_mm(f, q); mg, itg = x8_mm(g, q)
        if mf is None or mg is None:
            check(mf is None and mg is None, "C9.6 zero case agrees"); continue
        check(chi(mf, q) == chi(mg, q), "C9.6 same character of mm")
        same_mm += (mf == mg); diff_it += (itf != itg); n_ok += 1
        W1f = sub(scal(red(f, q), 10, q), shift(x8_deriv(f, q), 1), q)
        W1g = sub(scal(red(g, q), 10, q), shift(x8_deriv(g, q), 1), q)
        check(translate_mod(W1f, tt, q) == sub(W1g, scal(x8_deriv(g, q), tt, q), q), "C9.6 tau_t(W_f) = W_g - t g'")
        noncov += (W1g != translate_mod(W1f, tt, q))
print(f"C9.6 even f with c4 = 0 against f(x+t): {n_ok} pairs; same chi(mm) in all, same mm in {same_mm}, "
      f"different iteration counts in {diff_it}, first W not the translate in {noncov}")

# C9.7 probability that a call of s polynomials sets prefgen: P(Bin(s, p) > s/2), p = 1/4 + (3/4)(1/100)
from fractions import Fraction
p = Fraction(1, 4) + Fraction(3, 4) * Fraction(1, 100)
for s in (1, 7, 63, 64, 65, 127, 255, 1000):
    P = sum(Fraction(comb(s, k)) * p**k * (1 - p)**(s - k) for k in range(s + 1) if 2 * k > s)
    print(f"C9.7 s = {s:5d}: P(Bin(s, {float(p):.4f}) > s/2) = {float(P):.3g}")

# C9.8 Tschirnhaus transforms of h = H(x^r) by g of degree d: pi = 0 when 2d < r (Proposition "Tschirnhaus transforms").
# Power sums of h by Newton's identities (exact integers); then p1(f) = sum g_j p_j, p2(f) = sum g_j g_k p_(j+k), pi = p1^2 - n p2.
def x8_powersums(h, K):
    n = deg(h); p = [n]
    for k in range(1, K + 1):
        s = -sum(h[n - i] * p[k - i] for i in range(1, min(k - 1, n) + 1))
        if k <= n: s -= k * h[n - k]
        p.append(s)
    return p
def x8_tsch_pi(h, g):
    n = deg(h); p = x8_powersums(h, 2 * deg(g))
    p1 = sum(gj * p[j] for j, gj in enumerate(g))
    p2 = sum(gj * gk * p[j + k] for j, gj in enumerate(g) for k, gk in enumerate(g))
    return p1 * p1 - n * p2
seed = [-4, 0, 0, 0, 0, -12, 0, 0, 0, 0, 1]                        # x^10 - 12x^5 - 4
check(x8_powersums(seed, 10) == [10, 0, 0, 0, 0, 60, 0, 0, 0, 0, 760], "C9.8 power sums of the seed")
check(all(x8_tsch_pi(seed, [c, b, a]) == 0 for a in (1, 2, 3) for b in range(-4, 4) for c in range(-4, 4)), "C9.8 the 192 substitutions")
check(x8_tsch_pi(seed, [1, 0, 1, 1]) == -1200, "C9.8 x^3 + x^2 + 1 breaks it")
zero = nonzero = 0
for _ in range(2000):
    r = RNG.randrange(2, 7); m = RNG.randrange(1, 4)
    H = [RNG.randint(-9, 9) for _ in range(m)] + [1]
    h = [0] * (r * m + 1)
    for k, a in enumerate(H): h[r * k] = a
    d = RNG.randrange(1, 5)
    g = [RNG.randint(-5, 5) for _ in range(d)] + [RNG.randint(1, 5)]
    pv = x8_tsch_pi(h, g)
    if 2 * d < r:
        check(pv == 0, f"C9.8 pi = 0 for 2d < r, h={h} g={g}"); zero += 1
    else:
        nonzero += (pv != 0)
print(f"C9.8 Tschirnhaus: seed power sums, the 192 substitutions and x^3+x^2+1 checked; {zero} random (h = H(x^r), 2d < r) with pi = 0, "
      f"and pi != 0 in {nonzero} of {2000 - zero} cases with 2d >= r")

print("\nFAILURES:", len(FAIL))
sys.exit(1 if FAIL else 0)
