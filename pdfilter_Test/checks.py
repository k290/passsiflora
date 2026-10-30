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
# section C9 checks Section 10 of the document (zeroing the x^8 coefficient) with its own random generator;
# section C10 checks Section 8.8 (the size limits of the integer arithmetic) with another one; C11 checks the source-line citations
# and C12 the form of the ledger (python3 checks.py --cites / --ledger run one of them alone; see below).
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

# ---------- C11, defined here so that it can also run alone:
#   python3 checks.py --cites [--upstream DIR]          only the citation check (C11)
#   python3 checks.py --fix-cites [REV] [--upstream DIR] renumber the citations after a source edit: REV (default HEAD) is the
#                                                        version the citations were last right for; then run --cites
# pdfilter_chains.tex cites source lines as \Lc{file}{n}{anchor} and \Lcr{file}{n}{m}{anchor n}{anchor m}; an anchor is a text
# that must occur on the cited line (BLANK: an empty line, RBRACE: a line "}").  The pasted listings (lstlisting with firstnumber)
# must equal the cited range of the source.  DIR is a patched upstream GetDecicsSrc (for the citations of TgtMartinet.cpp).
import os, re, subprocess, difflib
CITE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CITE_TEX = os.path.join(CITE_ROOT, "pdfilter_Test", "pdfilter_chains.tex")
CITE_FILES = {"F": "src/pdfilter.f90", "I": "src/chain_general.inc", "FP": "327b04f:src/pdfilter.f90",
              "FPI": "327b04f:src/chain_general.inc", "TL": "pdfilter_Test/test_pdfilter_linux.cpp",
              "TT": "pdfilter_Test/test_pdfilter_tschirnhaus.cpp", "CT": "pdfilter_Test/chain_trace.cpp",
              "BL": "build/build-linux.sh", "BLFP": "327b04f:build/build-linux.sh", "OCL": "src/nf_gpu_opencl.cpp", "TGT": "upstream:TgtMartinet.cpp",
              "PD": "upstream:polDiscTest_cpu.cpp", "K": "src/nf_kernel.h", "GC": "src/nf_gpu_common.cpp"}
CITE_UPSTREAM = sys.argv[sys.argv.index("--upstream") + 1] if "--upstream" in sys.argv[1:-1] else os.environ.get("NF_UPSTREAM_SRC")
RE_LC = re.compile(r"\\Lc\{(\w+)\}\{(\d+)\}\{([^{}]*)\}")
RE_LCR = re.compile(r"\\Lcr\{(\w+)\}\{(\d+)\}\{(\d+)\}\{([^{}]*)\}\{([^{}]*)\}")
RE_LST = re.compile(r"\\begin\{lstlisting\}\[style=fortran,firstnumber=(\d+),caption=\{\\code\{src/([^}]*)\}, lines (\d+)--(\d+):")

def cite_source(key, rev=None):
    spec = CITE_FILES[key]
    try:
        if spec.startswith("upstream:"):
            if not CITE_UPSTREAM: return None
            return open(os.path.join(CITE_UPSTREAM, spec[9:])).read().split("\n")
        if ":" in spec or rev:
            r, p = spec.split(":") if ":" in spec else (rev, spec)
            out = subprocess.run(["git", "-C", CITE_ROOT, "show", f"{r}:{p}"], capture_output=True, text=True)
            return out.stdout.split("\n") if out.returncode == 0 else None
        return open(os.path.join(CITE_ROOT, spec)).read().split("\n")
    except OSError:
        return None

def anchor_on(lines, n, a):
    if not (1 <= n <= len(lines)): return False
    t = lines[n - 1]
    return t.strip() == "" if a == "BLANK" else t.strip() == "}" if a == "RBRACE" else a in t

def listings(tex):
    """(firstnumber, source key, caption first, caption last, body lines, start offset) of every pasted Fortran listing"""
    out = []
    for m in RE_LST.finditer(tex):
        key = {"pdfilter.f90": "F", "chain\\_general.inc": "I"}[m.group(2)]
        b0 = tex.index("\n", m.end()) + 1
        body = tex[b0:tex.index("\\end{lstlisting}", b0)].rstrip("\n").split("\n")
        out.append((int(m.group(1)), key, int(m.group(3)), int(m.group(4)), body, m.start()))
    return out

def check_cites():
    tex = open(CITE_TEX, encoding="utf-8").read()
    src = {k: cite_source(k) for k in CITE_FILES}
    n1 = n2 = nl = 0; skipped = set()
    for m in RE_LC.finditer(tex):
        k, n, a = m.group(1), int(m.group(2)), m.group(3)
        if src[k] is None: skipped.add(k); continue
        n1 += 1; check(anchor_on(src[k], n, a), f"C11 citation {CITE_FILES[k]} line {n}: '{a}' is not on it")
    for m in RE_LCR.finditer(tex):
        k, n, e, a, b = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), m.group(5)
        if src[k] is None: skipped.add(k); continue
        n2 += 1
        check(n < e and anchor_on(src[k], n, a) and anchor_on(src[k], e, b), f"C11 citation {CITE_FILES[k]} lines {n}--{e}: anchors '{a}' / '{b}'")
    for (first, k, c0, c1, body, _) in listings(tex):
        nl += 1
        check(c0 == first and c1 == first + len(body) - 1, f"C11 listing at line {first} of {CITE_FILES[k]}: caption {c0}--{c1}")
        check(src[k][first - 1:first - 1 + len(body)] == body, f"C11 listing {CITE_FILES[k]} {c0}--{c1} differs from the source")
    return n1, n2, nl, skipped

def auto_anchor(line):
    """the longest piece of the line without TeX-special characters, at most 48 characters"""
    t = line.strip()
    if not t: return "BLANK"
    if t == "}": return "RBRACE"
    best = max((s.strip() for s in re.split(r"[%#\\{}&~$^]", t)), key=len)
    if len(best) > 48:
        cut = best[:48]
        best = (cut[:cut.rfind(" ")] if " " in cut[20:] else cut).strip()
    return best

def fix_cites(base):
    tex = open(CITE_TEX, encoding="utf-8").read()
    maps, new = {}, {}
    for k, spec in CITE_FILES.items():
        if ":" in spec: continue                       # a fixed commit, or the upstream source: never renumbered
        old = open(os.path.join(base, spec)).read().split("\n") if os.path.isdir(base) else cite_source(k, base)
        cur = cite_source(k)
        if old is None or cur is None: continue
        mp = {}
        for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, old, cur, autojunk=False).get_opcodes():
            for d in range(i2 - i1):                   # equal lines map one to one; a rewritten block maps by position
                mp[i1 + d + 1] = (j1 + min(d, max(j2 - j1 - 1, 0)) + 1, tag == "equal")
        maps[k], new[k] = mp, cur
    changed, review = [], []
    def place(k, n, a):
        """(new line number, anchor) for the old line n with anchor a"""
        if k not in maps: return n, a
        cand, same = maps[k].get(n, (n, False))
        if anchor_on(new[k], cand, a): return cand, a
        if not same and cand <= len(new[k]):        # the cited line itself was rewritten: new anchor, to be reviewed
            a2 = auto_anchor(new[k][cand - 1])
            review.append(f"{CITE_FILES[k]} line {n} -> {cand}: '{a}' -> '{a2}'")
            return cand, a2
        hits = [i + 1 for i in range(len(new[k])) if anchor_on(new[k], i + 1, a)]
        if not hits: raise SystemExit(f"fix-cites: anchor '{a}' ({CITE_FILES[k]} line {n}) no longer occurs")
        return min(hits, key=lambda h: abs(h - cand)), a
    def one(m):
        k, n, a = m.group(1), int(m.group(2)), m.group(3)
        n2, a2 = place(k, n, a)
        if n2 != n: changed.append((k, n, n2))
        return f"\\Lc{{{k}}}{{{n2}}}{{{a2}}}"
    def rng(m):
        k, n, e, a, b = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4), m.group(5)
        (n2, a2), (e2, b2) = place(k, n, a), place(k, e, b)
        if (n2, e2) != (n, e): changed.append((k, f"{n}--{e}", f"{n2}--{e2}"))
        return f"\\Lcr{{{k}}}{{{n2}}}{{{e2}}}{{{a2}}}{{{b2}}}"
    tex = RE_LCR.sub(rng, RE_LC.sub(one, tex))
    for (first, k, c0, c1, body, pos) in reversed(listings(tex)):
        cur = new.get(k)
        if cur is None: continue
        hits = [i + 1 for i in range(len(cur) - len(body) + 1) if cur[i:i + len(body)] == body]
        if len(hits) != 1: raise SystemExit(f"fix-cites: the listing {CITE_FILES[k]} {c0}--{c1} occurs {len(hits)} times in the source")
        f2 = hits[0]
        if f2 != first:
            head = tex[pos:tex.index("\n", pos)]
            head2 = head.replace(f"firstnumber={first},", f"firstnumber={f2},").replace(f"lines {c0}--{c1}:", f"lines {f2}--{f2 + len(body) - 1}:")
            tex = tex[:pos] + head2 + tex[pos + len(head):]
            changed.append((k, f"listing {c0}--{c1}", f"{f2}--{f2 + len(body) - 1}"))
    open(CITE_TEX, "w", encoding="utf-8").write(tex)
    print(f"fix-cites (base {base}): {len(changed)} citations renumbered")
    for c in changed: print("  ", c)
    for r in review: print("   REVIEW (the cited line was rewritten; does the text still mean it?):", r)

# ---------- C12, the ledger (appendix "The ledger"):  python3 checks.py --ledger   prints it as text; every run checks its form:
# four cells per row, a status macro of the legend in the second, no unfilled placeholder anywhere in the document.
LEDGER_STATUS = {"\\stproved": "proved", "\\stquoted": "quoted", "\\stnum": "checked numerically", "\\streport": "reported",
                 "\\stassume": "named assumption"}
PLACEHOLDERS = ("POLDISCLEDGER", "TINYLEDGER", "WULEDGER", "GPULEDGER", "POLDISCROWS", "POLDISCTEXT", "GPUSHIM", "TINYROWS",
                "TINYTEXT", "GPURESULTS", "CPUSPEED")   # the markers used while the document was written, before the runs had finished

def ledger_rows():
    tex = open(CITE_TEX, encoding="utf-8").read()
    body = tex[tex.index("% LEDGER-BEGIN") + len("% LEDGER-BEGIN"):tex.index("% LEDGER-END")]
    return [[c.strip() for c in re.split(r"(?<!\\)&", r)] for r in body.split("\\\\\n") if r.strip()]

def detex(s):
    for m, w in LEDGER_STATUS.items(): s = s.replace(m + "{}", "[" + w + "]").replace(m, "[" + w + "]")
    s = re.sub(r"\\(?:ref|label)\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\cite(?:\[([^]]*)\])?\{([^}]*)\}", lambda m: "[" + m.group(2) + (", " + m.group(1) if m.group(1) else "") + "]", s)
    for _ in range(3): s = re.sub(r"\\(?:code|textbf|emph|texttt)\{([^{}]*)\}", r"\1", s)
    s = s.replace("~", " ").replace("\\_", "_").replace("{,}", ",").replace("\\ ", " ").replace("\\,", " ")
    s = re.sub(r"\\(le|ge|equiv|cdot|nmid|pmod|bmod|dots|chi|mid)(?![a-zA-Z])", lambda m: {"le": "<=", "ge": ">=", "equiv": "=", "cdot": "*", "dots": "..."}.get(m.group(1), m.group(1)), s)
    s = re.sub(r"\\[a-zA-Z]+", "", s)
    return re.sub(r"\s+", " ", s.replace("$", "").replace("{", "").replace("}", "")).strip()

def check_ledger():
    rows = ledger_rows()
    for i, r in enumerate(rows):
        check(len(r) == 4, f"C12 ledger row {i + 1}: {len(r)} cells")
        check(any(m in r[1] for m in LEDGER_STATUS), f"C12 ledger row {i + 1}: no status macro in '{r[1]}'")
    tex = open(CITE_TEX, encoding="utf-8").read()
    for p in PLACEHOLDERS:
        check(re.search(r"\b" + p + r"\b", tex) is None, f"C12 placeholder {p} left in the document")
    count = {}
    for r in rows:
        for m, w in LEDGER_STATUS.items():
            if m in r[1]: count[w] = count.get(w, 0) + 1
    return rows, count

if "--ledger" in sys.argv[1:]:
    rows, count = check_ledger()
    for i, r in enumerate(rows):
        print(f"[{i + 1}] {detex(r[0])}\n      status:   {detex(r[1])}\n      checked:  {detex(r[2])}\n      rests on: {detex(r[3])}")
    print(f"C12 ledger: {len(rows)} rows; " + ", ".join(f"{w} {n}" for w, n in sorted(count.items())))
    print("FAILURES:", len(FAIL))
    sys.exit(1 if FAIL else 0)

if "--fix-cites" in sys.argv[1:]:
    i = sys.argv.index("--fix-cites")
    fix_cites(sys.argv[i + 1] if i + 1 < len(sys.argv) and not sys.argv[i + 1].startswith("-") else "HEAD")
if "--cites" in sys.argv[1:] or "--fix-cites" in sys.argv[1:]:
    c1, c2, cl, sk = check_cites()
    print(f"C11 citations: {c1} lines and {c2} ranges with their anchors, {cl} listings; files not available: {sorted(sk) or 'none'}")
    print("FAILURES:", len(FAIL))
    sys.exit(1 if FAIL else 0)

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
small, q = [], 2 ** 13 - 1
while len(small) < 24 and q > 256:
    if usable(q): small.append(q)
    q -= 2
large, q = [], 2 ** 25 - 1
while len(large) < 6 and q > 2 ** 24:
    if usable(q): large.append(q)
    q -= 2
FP = small + large
print("filter primes: small", small, "large", large[:3], "...", large[-1], "count", len(FP), "all 1 or 9 mod 40:", all(p % 40 in (1, 9) for p in FP))
# The sets S of the project's searches (decic fields with S = {2,5}; finished ones with {2,11} and {5,11}; batch status pages
# of NumberFields@home, read 2026-10-30): the small range has at least the 24 primes the filter takes (Section "The prime split").
for S in ((2, 5), (2, 11), (5, 11)):
    use = lambda q: q % 4 == 1 and is_prime(q) and all(pow(p, (q - 1) // 2, q) == 1 for p in S)
    u13 = [q for q in range(2 ** 13 - 1, 256, -2) if use(q)]
    u10 = [q for q in range(2 ** 10 - 1, 256, -2) if use(q)]
    check(len(u13) >= 24, f"C7 S={S}: only {len(u13)} usable small primes")
    print(f"C7 S = {S}: {len(u13)} usable primes in (256, 2^13) ({u13[0]} ... the 24th {u13[23]}, the 30th {u13[29]}), {len(u10)} in (256, 2^10)")
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
print(f"C7 GPU table for S = (2, 5): {gpu25[0]} ... {gpu25[-1]} (24 primes, none of the CPU's small ones: {not set(gpu25) & set(small)});"
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

# ---------- C9: Section 10, zeroing the x^8 coefficient
# Its own random generator, so the counts of C0-C8 do not change; the same cases in quick and full runs.
# The draws follow the order used for the numbers quoted in Section 10 of pdfilter_chains.tex.
print("\n===== C9: zeroing the x^8 coefficient (Section 10) =====")
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
    RNG.random()               # one draw, kept so that the counts quoted in Section 10 stay reproducible
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

# ---------- C10: Section 8.8, the size limits of the integer arithmetic (own random generator; C0-C9 output unchanged)
print("\n===== C10: size limits of the integer arithmetic (Section 8.8) =====")
RNG10 = random.Random(20261026)

def s_k(v, k):
    v %= 1 << k
    return v - (1 << k) if v >= 1 << (k - 1) else v

def mred(x, q, k):
    """signed Montgomery reduction with R = 2^k as defined in Lemma 8.1: m = s_k(x q'), u = (x - m q) / R"""
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

# C10.1 the characterisation (Lemma "When the result leaves (-q,q)", Section 8.2): u >= q iff x = (R+m)q + Rt (t >= 0); every (m', t) gives u = q + t; same for u <= -q
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

# C10.3 sredc on differences (Lemma "Where sredc can fail on a difference", Section 8.8): the reduction, the bound 2(delta-1)^2 >= 2^30, and the complete search of the remaining range
def sredc_fail_search(delta, neg=False):
    """all solutions of U + V + 2 + 4 t delta = L q, mu = 2 delta - 4 - U + L - 4t >= (1 if neg else 0), V = ij + i'j', i+j+i'+j' = U"""
    q = (1 << 30) + delta
    for L in (1, 2, 3):
        for U in range(0, 2 * delta - 3 + L):
            tmax = (2 * delta - 4 + L - U - (1 if neg else 0)) // 4
            if tmax < 0: break
            top = L * q - U - 2; vmax = U * U // 4
            if top < 0: continue
            tlo = max(0, -(-(top - vmax) // (4 * delta))); thi = min(tmax, top // (4 * delta))
            for t in range(tlo, thi + 1):
                V = top - 4 * t * delta
                if not (0 <= V <= vmax): continue
                for s in range(U + 1):
                    s2 = U - s
                    for i in range(s // 2 + 1):
                        r = V - i * (s - i)
                        if r < 0: break
                        disc = s2 * s2 - 4 * r
                        if disc < 0: continue
                        sq = int(disc ** 0.5)
                        while sq * sq > disc: sq -= 1
                        while (sq + 1) ** 2 <= disc: sq += 1
                        if sq * sq == disc and (s2 - sq) % 2 == 0:
                            return (L, U, t, V, i, s - i, (s2 - sq) // 2, s2 - (s2 - sq) // 2)
    return None

for d in range(1, 23172, 2):
    if 2 * (d - 1) ** 2 >= 1 << 30: check(False, "C10.3 the bound 2(delta-1)^2 < 2^30 below 23173"); break
first = None; failing = []
for d in range(23173, 23379, 2):
    w = sredc_fail_search(d)
    if w:
        failing.append(d)
        if first is None: first = (d, w)
check(first is not None and first[0] == 23279, "C10.3 first failing odd q = 2^30 + 23279")
fprime = [d for d in failing if is_prime((1 << 30) + d)]
check(fprime[0] == 23293 and ((1 << 30) + 23293) % 4 == 1, "C10.3 first failing prime 2^30 + 23293, and it is 1 mod 4")
check([d for d in failing if usable((1 << 30) + d)][0] == 23377, "C10.3 first failing prime usable for S = {2,5}: 2^30 + 23377")
d, (L, U, t, V, i, j, i2_, j2) = first
q = (1 << 30) + d; n = q - 1
a, b, c, dd = n - i, n - j, -(n - i2_), n - j2
x = a * b - c * dd
u, m = mred(x, q, 32)
check(u >= q and x == 2 * n * n - (U * n - V), "C10.3 the witness fails and satisfies the reduction")
check(sredc_fail_search(d, neg=True) is not None, "C10.3 the negative side fails at the same q")
reps = [(1073765103, 1073765102, 1073765101, -1073764896, 1073764888), (1073765117, 1073765111, 1073765111, -1073764914, 1073764866),
        (1073765201, 1073765200, 1073765198, -1073765108, 1073765063), (1073774593, 1073741825, 1073741824, -1073774592, 1073774592)]
for (qq, a1, b1, c1, d1) in reps:
    uu, _ = mred(a1 * b1 - c1 * d1, qq, 32)
    check(uu >= qq and max(abs(a1), abs(b1), abs(c1), abs(d1)) <= qq - 1, f"C10.3 counterexample at q = {qq}")
check(mred((1 << 31) * 1073774593, 1073774593, 32)[0] == 1073774593, "C10.3 x = 2^31 q returns +q (a true zero)")
print(f"C10.3 sredc: no failure for odd q = 2^30 + delta, delta <= 23171 (bound), none for delta in [23173, 23277] (search), "
      f"first at delta = 23279: a,b,c,d = {a}, {b}, {c}, {dd}, sredc = q + {u - q}; first prime delta = 23293 (= 1 mod 4), "
      f"first usable 23377; {len(failing)} of the {len(range(23279, 23379, 2))} odd delta in [23279, 23377] fail; 4 counterexamples replayed")

# C10.4 the bound propagation of Section 8.8 (Theorem "The limits of the whole filter"): small stage for every odd q <= 9001, large stage at its limit
def red_bound(X, q, R): return (X + q * (R // 2)) // R
def load_bound(q): return red_bound(2 * ((1 << 21) - 1) * (q - 1) + (1 << 21) * (q - 1), q, 1 << 32)
def bounds_ok(q, small):
    R = 1 << 16 if small else 1 << 32
    def red(X):
        if X >= q * (R // 2): raise ValueError
        return red_bound(X, q, R)
    def mo(x, y): return red(x * y)
    def df(x, y):
        if small and x + y > 32767: raise ValueError
        return x + y
    def pr(x1, y1, x2, y2): return df(mo(x1, y1), mo(x2, y2)) if small else red(x1 * y1 + x2 * y2)
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
fails = [q for q in range(3, 9002, 2) if not bounds_ok(q, True)]
check(fails[0] == 8387 and fails == list(range(8387, 9002, 2)), "C10.4 small stage: passes exactly for odd q <= 8385 (up to 9001)")
check(bounds_ok(1097010265, False) and not bounds_ok(1097010267, False), "C10.4 large stage: limit 1097010265")
samp = [RNG10.randrange(3, 1097010265) | 1 for _ in range(300 if not QUICK else 30)]
check(all(bounds_ok(qq, False) for qq in samp), "C10.4 large stage: random odd q below the limit pass")
print(f"C10.4 bound propagation: small stage passes for every odd q <= 8385 and fails for every odd q in [8387, 9001]; "
      f"large stage passes at 1097010265 and {len(samp)} random odd q below, fails at 1097010267")

# C10.5 overflow limits: a difference of two mont results (inputs |a|,|b| <= 2q-1) leaves 16 bits from q = 12801; sredc's x - m q
# can leave 64 bits only if 2(q-1)^2 + 2^31 q >= 2^63
def mont_extremes(q):
    B = 2 * q - 1; best = bestn = 0; R = 1 << 16
    for a in range(B, 0, -1):
        if a * B < R * min(best, bestn) - q * R // 2: break
        for b in range(B, 0, -1):
            if a * b < R * min(best, bestn) - q * R // 2: break
            best = max(best, mont16(a, b, q)); bestn = max(bestn, -mont16(a, -b, q))
    return best, bestn
e1, e2 = mont_extremes(12799), mont_extremes(12801)
check(sum(e1) < 32768 <= sum(e2), "C10.5 i2 overflow of a difference of two mont results first possible at q = 12801")
lo = max(qq for qq in range(1676704560, 1676704570) if 2 * (qq - 1) ** 2 + (1 << 31) * qq < 1 << 63)
check(lo == 1676704566, "C10.5 sredc: no 64-bit overflow for q <= 1676704566")
print(f"C10.5 max - min of mont over |a|,|b| <= 2q-1: {sum(e1)} at q = 12799, {sum(e2)} at q = 12801 (16-bit overflow of a "
      f"difference from 32768); sredc stays within 64 bits for q <= {lo}")

# ---------- C11: the source-line citations and the pasted listings of pdfilter_chains.tex (the functions are at the top)
print("\n===== C11: source-line citations of the document =====")
c1, c2, cl, sk = check_cites()
print(f"C11 citations: {c1} lines and {c2} ranges with their anchors, {cl} listings; files not available: {sorted(sk) or 'none'}")

# ---------- C12: the ledger (the functions are at the top; --ledger prints it)
print("\n===== C12: the ledger of the document =====")
rows, count = check_ledger()
print(f"C12 ledger: {len(rows)} rows; " + ", ".join(f"{w} {n}" for w, n in sorted(count.items())))

print("\nFAILURES:", len(FAIL))
sys.exit(1 if FAIL else 0)
