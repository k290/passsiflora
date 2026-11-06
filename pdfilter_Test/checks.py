#!/usr/bin/env python3
# checks.py -- exact-integer checks of the statements about resultants, discriminants and the two
# Euclidean chains that docs/whitepaper.tex proves (Section 2 "The rejection rule", Section 3 "Two Euclidean chains for
# one Legendre symbol", Section 4 "Exact integer arithmetic").
# Independent of chain_trace.cpp: own polynomial arithmetic, own models of chain_block and
# chain_general, resultants as Sylvester determinants (a quoted standard result), checked against
# the root definition for split polynomials.  Python 3 standard library only, fixed seed.
#
# Usage:   python3 checks.py [quick] [--gp FILE]
#   quick      a tenth of the cases (a few seconds; the full run takes about a minute)
#   --gp FILE  also write a PARI/GP script that compares 300 of the discriminants with poldisc():
#              gp -q FILE < /dev/null   prints "gp: disc samples checked" and nothing else if all agree
# Exit status 1 if any check fails.  Section C10 checks the size limits of the integer arithmetic (whitepaper, Section 4.4
# "The size limits, and the guard") with its own random generator.  Two checks tie the documents to the code: C11 that the
# limits in src/pdfilter.f90 and src/nf_gpu_common.cpp are the proved ones, C12 that every citation of the whitepaper in the
# repository resolves, with the title of every section or table it cites by number.
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
    """The resultant by recursion on degrees: the rules (R1)-(R3) of the whitepaper, Section 2 "The rejection rule"."""
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
print(f"C0 Res = det Syl = root definition (split V) and = the recursion (R1)-(R3): {cnt} cases")

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
small, q = [], 2 ** 13 - 1          # pdf_init's rule: the 30 largest usable primes below 2^13
while len(small) < 30 and q > 256:
    if usable(q): small.append(q)
    q -= 2
FP = small
print("filter primes:", FP[:3], "...", FP[-1], "count", len(FP), "all 1 or 9 mod 40:", all(p % 40 in (1, 9) for p in FP))
# The sets S of the project's searches (decic fields with S = {2,5}; finished ones with {2,11} and {5,11}; batch status pages
# of NumberFields@home, read 2026-10-30): the range has at least the 30 primes the filter takes (whitepaper, Section 9.2, "The prime split").
for S in ((2, 5), (2, 11), (5, 11)):
    use = lambda q: q % 4 == 1 and is_prime(q) and all(pow(p, (q - 1) // 2, q) == 1 for p in S)
    u13 = [q for q in range(2 ** 13 - 1, 256, -2) if use(q)]
    u10 = [q for q in range(2 ** 10 - 1, 256, -2) if use(q)]
    check(len(u13) >= 30, f"C7 S={S}: only {len(u13)} usable primes below 2^13")
    print(f"C7 S = {S}: {len(u13)} usable primes in (256, 2^13) ({u13[0]} ... the 30th {u13[29]}), {len(u10)} in (256, 2^10)")
# The GPU table (buildPrimeTable of src/nf_gpu_common.cpp): the 24 smallest usable primes in (256, 8192); GPU mode needs at least 8.
# For every S of one to three primes below 60: how many sets leave fewer than 8, or fewer than 24, in (256, 1024) (the earlier range)
# and in (256, 8192).
import itertools
Q4 = [q for q in range(257, 8192, 4) if is_prime(q)]
gpu_short = {1024: [0, 0], 8192: [0, 0]}
small_p = [p for p in range(2, 60) if is_prime(p)]
for k in (1, 2, 3):
    for S in itertools.combinations(small_p, k):
        u = [q for q in Q4 if all(p % q and pow(p, (q - 1) // 2, q) == 1 for p in S)]
        for hi in (1024, 8192):
            n = sum(1 for q in u if q < hi)
            gpu_short[hi][0] += n < 8; gpu_short[hi][1] += n < 24
check(gpu_short[8192] == [0, 0], f"C7 GPU table: sets with fewer than 8 / 24 usable primes below 8192: {gpu_short[8192]}")
gpu25 = [q for q in Q4 if pow(2, (q - 1) // 2, q) == 1 and pow(5, (q - 1) // 2, q) == 1][:24]
print(f"C7 GPU table for S = (2, 5): {gpu25[0]} ... {gpu25[-1]} (24 primes, none of the CPU's: {not set(gpu25) & set(small)});"
      f" sets of 1-3 primes below 60 with fewer than 8 / 24 usable primes: {gpu_short[1024]} below 1024, {gpu_short[8192]} below 8192")

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
    for q in FP:
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

# ---------- C10: the size limits of the integer arithmetic (whitepaper, Section 4.4 "The size limits, and the guard"; own
# random generator, so the output of C0-C7 does not depend on it)
print("\n===== C10: size limits of the integer arithmetic =====")
RNG10 = random.Random(20261026)

def s_k(v, k):
    v %= 1 << k
    return v - (1 << k) if v >= 1 << (k - 1) else v

def mred(x, q, k):
    """signed Montgomery reduction with R = 2^k as in the whitepaper, Lemma "Montgomery reduction": m = s_k(x q'), u = (x - m q) / R"""
    R = 1 << k
    m = s_k(x * pow(q, -1, R), k)
    assert (x - m * q) % R == 0
    return (x - m * q) // R, m

def mont16(a, b, q):
    """the Fortran mont, step by step: lo, m, h1, h2 as 16/32-bit integers"""
    qinv = s_k(pow(q, -1, 1 << 16), 16)
    x = a * b
    lo = s_k(x, 16); m = s_k(lo * qinv, 16)
    h1 = s_k(x >> 16, 16); h2 = s_k((m * q) >> 16, 16)
    return s_k(h1 - h2, 16)

# C10.1 when the result leaves (-q,q): u >= q iff x = (R+m)q + Rt (t >= 0); every (m', t) gives u = q + t; same for u <= -q
cnt = 0
for k in (16, 32):
    R = 1 << k
    for _ in range(1000 if not QUICK else 100):
        q = RNG10.randrange(3, (1 << (k - 1)) - 1) | 1
        mp = RNG10.randrange(-R // 2, R // 2); t = RNG10.randrange(0, 50)
        x = (R + mp) * q + R * t
        u, m = mred(x, q, k)
        check(m == mp and u == q + t, "C10.1 x = (R+m')q + Rt gives m = m', u = q + t")
        x = (mp - R) * q - R * t
        u, m = mred(x, q, k)
        check(m == mp and u == -(q + t), "C10.1 x = (m'-R)q - Rt gives u = -(q + t)")
        x = RNG10.randrange(-3 * R * q, 3 * R * q)
        u, m = mred(x, q, k)
        check((u >= q) == (x - (R + m) * q >= 0), "C10.1 u >= q iff x >= (R+m) q")
        check((u <= -q) == (x - (m - R) * q <= 0), "C10.1 u <= -q iff x <= (m-R) q")
        cnt += 1
    # the smallest |x| that fails: (R/2) q on the positive side, (R/2 + 1) q on the negative one
    for q in (3, 41, 8191, 8193, (1 << 30) + 1, (1 << 31) - 1):
        if q >= R // 2: continue
        check(mred(q * R // 2, q, k)[0] == q, "C10.1 x = qR/2 gives u = q")
        check(abs(mred(q * R // 2 - 1, q, k)[0]) < q, "C10.1 x = qR/2 - 1 gives |u| < q")
        check(mred(-(R // 2 + 1) * q, q, k)[0] == -q, "C10.1 x = -(R/2+1)q gives u = -q")
print(f"C10.1 characterisation of |u| >= q for R = 2^16 and 2^32: {cnt} random q, each with a constructed x on both sides and a random x")

# C10.2 mont: the replayed counterexamples, and the scan of [8193, 8300] by the characterisation (positive and negative side)
check(mont16(-16481, -16473, 8241) == 8262, "C10.2 mont q=8241 a=-16481 b=-16473 gives 8262 = q + 21")
check(mont16(16481, -16473, 8241) == -8262, "C10.2 mont q=8241 negative side")
for (a, b) in ((-16481, -16473), (16473, 16481), (-16481, 16473)):
    u, _ = mred(a * b, 8241, 16)
    check(u == mont16(a, b, 8241), "C10.2 the Fortran steps equal the definition")

def mont_witness(q, neg):
    R = 1 << 16; B = 2 * q - 1
    ms = range(-R // 2, R // 2) if not neg else range(R // 2 - 1, -R // 2 - 1, -1)
    for m in ms:
        x0 = (R + m) * q if not neg else (R - m) * q
        if x0 > B * B: break
        for x in range(x0, B * B + 1, R):
            lo = -(-x // B)
            a = lo
            while a * a <= x:
                if x % a == 0 and x // a <= B: return (a, x // a if not neg else -(x // a))
                a += 1
    return None

clean = []
for q in range(8193, 8301, 2):
    wp, wn = mont_witness(q, False), mont_witness(q, True)
    if wp: check(mont16(wp[0], wp[1], q) >= q, f"C10.2 positive witness q={q}")
    if wn: check(mont16(wn[0], wn[1], q) <= -q, f"C10.2 negative witness q={q}")
    if not wp and not wn: clean.append(q)
check(clean == list(range(8193, 8240, 2)) + [8243, 8247, 8251, 8253, 8255, 8257, 8259], "C10.2 the odd q in [8193, 8300] without a counterexample")
check(min(q for q in range(8193, 8301, 2) if q not in clean) == 8241, "C10.2 first failing odd q = 8241")
print(f"C10.2 mont: counterexample at q = 8241 replayed; odd q in [8193, 8300] without a counterexample: {len(clean)} (8193..8239, 8243, 8247, 8251..8259)")

# C10.3 the bound propagation (whitepaper, Lemma "Bound propagation" and Theorem "The limits of the whole filter") for every odd q <= 9001
def red_bound(X, q, R): return (X + q * (R // 2)) // R
def load_bound(q): return red_bound(2 * ((1 << 21) - 1) * (q - 1) + (1 << 21) * (q - 1), q, 1 << 32)
def bounds_ok(q):
    R = 1 << 16
    def red(X):
        if X >= q * (R // 2): raise ValueError
        return red_bound(X, q, R)
    def mo(x, y): return red(x * y)
    def df(x, y):
        if x + y > 32767: raise ValueError
        return x + y
    def pr(x1, y1, x2, y2): return df(mo(x1, y1), mo(x2, y2))
    try:
        L = load_bound(q)
        A = [L] * 10 + [1]; B = [mo(A[k], k) for k in range(1, 11)] + [0]
        for mm in range(9, 0, -1):
            lc, top = B[mm], A[mm + 1]
            Rr = [0] * 11; P = [0] * 11
            Rr[0] = mo(lc, A[0])
            for k in range(1, mm + 1): Rr[k] = pr(lc, A[k], top, B[k - 1])
            t = Rr[mm]
            for k in range(mm): P[k] = pr(lc, Rr[k], t, B[k])
            A, B = B, P
        mo(B[0], B[0]); mo(q - 1, B[0]); mo(q - 1, q - 1)
        bu = max(L, 1); bv = mo(10, L); mmb = 1
        for _ in range(18):
            bw = pr(bv, bu, bu, bv)
            mmb = mo(max(mmb, q - 1), max(bu, bv, bw))
            bu = bv = max(bu, bv, bw)
        mo(q - 1, bu); mo(q - 1, q - 1)
    except ValueError:
        return False
    return True
fails = [q for q in range(3, 9002, 2) if not bounds_ok(q)]
check(fails[0] == 8387 and fails == list(range(8387, 9002, 2)), "C10.3 passes exactly for odd q <= 8385 (up to 9001)")
print(f"C10.3 bound propagation: passes for every odd q <= 8385 and fails for every odd q in [8387, 9001]")

# C10.4 overflow limit: a difference of two mont results (inputs |a|,|b| <= 2q-1) leaves 16 bits from q = 12801
def mont_extremes(q):
    B = 2 * q - 1; best = bestn = 0; R = 1 << 16
    for a in range(B, 0, -1):
        if a * B < R * min(best, bestn) - q * R // 2: break
        for b in range(B, 0, -1):
            if a * b < R * min(best, bestn) - q * R // 2: break
            best = max(best, mont16(a, b, q)); bestn = max(bestn, -mont16(a, -b, q))
    return best, bestn
e1, e2 = mont_extremes(12799), mont_extremes(12801)
check(sum(e1) < 32768 <= sum(e2), "C10.4 i2 overflow of a difference of two mont results first possible at q = 12801")
print(f"C10.4 max - min of mont over |a|,|b| <= 2q-1: {sum(e1)} at q = 12799, {sum(e2)} at q = 12801 (16-bit overflow of a "
      f"difference from 32768)")

# ---------- C11: the limits in the source are the proved ones.  pdf_init's guard (QMIN, QMAX in src/pdfilter.f90) must be the
# limits C10.3 derives, its prime search (QSTART, QSTOP) must lie inside the guard, and the GPU
# table's range (buildPrimeTable in src/nf_gpu_common.cpp) inside the 11 <= q < 8192 of the kernel's proof.  A changed constant
# fails here until the proof is redone.
import os, re
print("\n===== C11: the limits in the source =====")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def src_text(rel): return open(os.path.join(ROOT, rel), encoding="utf-8", errors="replace").read()
F90 = src_text("src/pdfilter.f90")
def f90_const(name):   # the value of an integer parameter, e.g. "QSTART = 2_i8**13 - 1"
    m = re.search(r"parameter\s*::.*?\b" + name + r"\s*=\s*([0-9_i8*+\- ]+?)\s*(?:,|!|$)", F90, re.M)
    check(m is not None, f"C11 parameter {name} not found in src/pdfilter.f90")
    expr = m.group(1).replace("_i8", "") if m else "0"
    return eval(expr, {"__builtins__": {}}) if re.fullmatch(r"[0-9*+\- ]+", expr) else 0
K = {n: f90_const(n) for n in ("QMIN", "QMAX", "QSTART", "QSTOP", "MAXS")}
limit = fails[0] - 2                            # C10.3: every odd q up to here passes, the next fails
check(K["QMAX"] == limit, f"C11 QMAX = {K['QMAX']} is not the proved limit {limit}")
check(K["QMIN"] >= 11, f"C11 QMIN = {K['QMIN']}: the proofs need q >= 11 (q does not divide 10, the constants 1..10 lie below 2q)")
check(K["QMIN"] <= K["QSTOP"] < K["QSTART"] <= K["QMAX"], f"C11 the prime search ({K['QSTOP']}, {K['QSTART']}] leaves the guard")
gm = re.search(r"for \(long long q = (\d+); q < (\d+)", src_text("src/nf_gpu_common.cpp"))
check(gm is not None, "C11 the prime loop of buildPrimeTable was not found in src/nf_gpu_common.cpp")
glo, ghi = (int(gm.group(1)), int(gm.group(2))) if gm else (0, 0)
check(glo >= 11 and ghi <= 8192, f"C11 the GPU table's range [{glo}, {ghi}) leaves 11 <= q < 8192 (Lemma redc, Theorem E)")
print(f"C11 guard QMIN = {K['QMIN']}, QMAX = {K['QMAX']} (the proved limits); search ({K['QSTOP']}, {K['QSTART']}] "
      f"for {K['MAXS']} primes; GPU table [{glo}, {ghi})")

# ---------- C12: every citation of docs/whitepaper.tex in the comments and documents of the repository names something that
# exists there: a section or appendix by its number and the quoted title of it or of a subsection or paragraph in it, a section
# title, a lettered theorem, the quoted title of a theorem, proposition, lemma, definition or corollary, a paragraph title, or a
# table by its number and the quoted opening words of its caption.  A number without its title fails, since a renumbering can
# leave it valid but pointing elsewhere.  A citation is looked for in the text around each mention of the whitepaper, with
# comment markers removed and lines joined, so that citations which wrap are found too; a quoted theorem, proposition, lemma,
# definition or corollary is looked for anywhere, since the repository cites only the whitepaper's.  This file's own C11 and
# C12 are not scanned.
print("\n===== C12: citations of the whitepaper =====")
TEX = src_text("docs/whitepaper.tex")
def brace_arg(s, i):   # s[i] is "{": the balanced argument and the index after it
    depth = 0
    for j in range(i, len(s)):
        depth += (s[j] == "{") - (s[j] == "}")
        if depth == 0: return s[i + 1:j], j + 1
    return s[i + 1:], len(s)
def plain(t):          # TeX or cited text -> lower case words: \code{x} -> x, \_ -> _, no $, ~, `, \, or other macros
    t = re.sub(r"\\(?:code|emph|texttt|textsf|mathrm)\{([^{}]*)\}", r"\1", t)
    t = t.replace("\\_", "_").replace("~", " ").replace("\\,", " ").replace("$", "").replace("`", "").replace("''", "")
    t = re.sub(r"\\[a-zA-Z]+\*?", " ", t).replace("{", "").replace("}", "")
    return " ".join(t.split()).lower()
heads = {}             # number -> (title, start, end) for every section and subsection
appendix_at = TEX.find("\\appendix")
marks = [(m.start(), m.group(1), m.end() - 1) for m in re.finditer(r"\\(section|subsection)\{", TEX)]
sec = sub = 0; in_app = False
for k, (pos, kind, brace) in enumerate(marks):
    title, _ = brace_arg(TEX, brace)
    if kind == "section":
        if appendix_at >= 0 and pos > appendix_at and not in_app: in_app, sec = True, 0
        sec += 1; sub = 0
        num = chr(ord("A") + sec - 1) if in_app else str(sec)
        end = next((p for p, kd, _ in marks[k + 1:] if kd == "section"), len(TEX))
    else:
        sub += 1
        num = f"{chr(ord('A') + sec - 1) if in_app else sec}.{sub}"
        end = marks[k + 1][0] if k + 1 < len(marks) else len(TEX)
    heads[num] = (plain(title), pos, end)
envs = []              # (kind, letter or "", plain title)
for m in re.finditer(r"\\begin\{(proposition|lemma|theorem|corollary|definition|remark|thm([A-E]))\}(?:\[([^\]]*)\])?", TEX):
    envs.append(("theorem" if m.group(2) else m.group(1), m.group(2) or "", plain(m.group(3) or "")))
tables = [plain(brace_arg(TEX, m.end() - 1)[0]) for m in re.finditer(r"\\caption\{", TEX)]
PARA = re.compile(r"\\paragraph\{")
def titles_in(lo, hi):           # the titles of the sections, subsections and paragraphs that start in TEX[lo:hi], without a final "."
    return ({h[0] for h in heads.values() if lo <= h[1] < hi}
            | {plain(brace_arg(TEX, m.end() - 1)[0]).rstrip(".") for m in PARA.finditer(TEX, lo, hi)})
TITLES = titles_in(0, len(TEX))
def resolves(kind, a, b=None):   # True if the citation exists
    if kind in ("section", "appendix"):   # a: number or letter; b: quoted title (optional) of it or of a part of it
        return a in heads and (b is None or plain(b).rstrip(".") in titles_in(heads[a][1], heads[a][2]))
    if kind == "sectitle": return any(plain(a) == h[0] for h in heads.values())
    if kind == "letter": return any(e[1] == a for e in envs)
    if kind == "named": return any(e[0] == a.lower() and plain(b) in e[2] for e in envs)
    if kind == "para": return plain(b).rstrip(".") in TITLES
    if kind == "table":          # a: number; b: quoted opening words of the caption (optional)
        return 1 <= int(a) <= len(tables) and (b is None or tables[int(a) - 1].startswith(plain(b).rstrip(".")))
NUM = r"([0-9]+(?:\.[0-9]+)?|[A-Z](?:\.[0-9]+)?)\b"
DASH = "(?:-|--|–|—|and|to)"
TITLE = r"(?:[\s,]*\(?\s*\"([^\"]+)\")?"   # the quoted title after a number, if there is one
PATTERNS = [
    ("section", re.compile(r"\bSections?\s+" + NUM + r"(?:\s*" + DASH + r"\s*" + NUM + r")?"
                           r"(?:\s+of\s+`?(?:docs/)?whitepaper(?:\.tex)?`?)?" + TITLE)),
    ("sectitle", re.compile(r"\bSection\s+\"([^\"]+)\"")),
    ("letter", re.compile(r"\bTheorems?\s+([A-E])(?:\s*" + DASH + r"\s*([A-E]))?\b")),
    ("named", re.compile(r"\b(Theorem|Proposition|Lemma|Definition|Corollary)\s+\"([^\"]+)\"")),
    ("para", re.compile(r"whitepaper(?:\.tex)?`?'s\s+\"([^\"]+)\"")),
    ("table", re.compile(r"\bTables?\s+([0-9]+)\b" + TITLE)),
    ("appendix", re.compile(r"\bAppendix\s+([A-Z])\b" + TITLE)),
]
SKIP_DIRS = {"dist", "pdfilter_Test/build", "pdfilter_Test/fault_build", "pdfilter_Test/arith_build", "pdfilter_Test/tiny_build",
             "pdfilter_Test/poldisc_build", "pdfilter_Test/repro/expected"}    # outputs of the build and test scripts
EXTS = (".f90", ".inc", ".cpp", ".h", ".c", ".cu", ".sh", ".py", ".md", ".ps1")
cites = []             # (file, kind, cited text, ok, number only)
for dirpath, dirnames, filenames in os.walk(ROOT):
    reldir = os.path.relpath(dirpath, ROOT).replace(os.sep, "/")
    dirnames[:] = [d for d in dirnames if not d.startswith(".") and (d if reldir == "." else f"{reldir}/{d}") not in SKIP_DIRS]
    for fn in filenames:
        if not fn.endswith(EXTS): continue
        path = os.path.join(dirpath, fn); rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
        text = open(path, encoding="utf-8", errors="replace").read()
        if rel == "pdfilter_Test/checks.py": text = text[:text.find("# ---------- C11:")]
        text = " ".join(re.sub(r"^\s*(?://+|/\*+|\*+/?|!+|#+)?\s*", "", ln) for ln in text.splitlines())
        seen = set()
        near = [w.start() for w in re.finditer(r"whitepaper", text)]
        for kind, rx in PATTERNS:
            for m in rx.finditer(text):
                if kind != "named" and not any(w - 60 <= m.start() < w + 220 for w in near): continue
                if (kind, m.group(0)) in seen: continue
                seen.add((kind, m.group(0)))
                t = m.group(rx.groups) if kind in ("section", "table", "appendix") else None   # the quoted title
                if kind == "section":
                    ok = resolves("section", m.group(1), t) and (m.group(2) is None or resolves("section", m.group(2)))
                elif kind == "letter":
                    ok = resolves("letter", m.group(1)) and (m.group(2) is None or resolves("letter", m.group(2)))
                elif kind == "named": ok = resolves("named", m.group(1), m.group(2))
                elif kind == "sectitle": ok = resolves(kind, m.group(1))
                elif kind == "para": ok = resolves("para", None, m.group(1))
                else: ok = resolves(kind, m.group(1), t)
                by_number = kind in ("section", "table", "appendix") and (t is None or kind == "section" and m.group(2) is not None)
                cites.append((rel, kind, m.group(0), ok, by_number))
                check(ok, f"C12 {rel}: '{m.group(0)}' does not resolve in docs/whitepaper.tex")
                check(not by_number, f"C12 {rel}: '{m.group(0)}' cites by number; quote the title of each section or table it cites")
files = sorted({c[0] for c in cites})
print(f"C12 {len(cites)} citations of the whitepaper in {len(files)} files: "
      + ", ".join(f"{k} {sum(1 for c in cites if c[1] == k)}" for k, _ in PATTERNS if any(c[1] == k for c in cites))
      + f"; {sum(c[3] for c in cites)} resolve, {sum(c[4] for c in cites)} give a number without its title; "
      + f"the whitepaper has {len(heads)} numbered sections, {len(envs)} theorem-like statements, {len(tables)} tables")
check(len(cites) >= 30, f"C12 only {len(cites)} citations found: the scan is broken")

print("\nFAILURES:", len(FAIL))
sys.exit(1 if FAIL else 0)
