#!/usr/bin/env python3
"""
nb_gate_count.py -- referee check of N_b for 'Kardashev's Conundrum' (Section 2.5).
Sources for the algorithm: NIST FIPS 180-4 (2015), Sec. 3.2, 4.1.2, 4.2.2,
5.1.1, 5.3.3, 6.2.2; RFC 6234 (2011). Header layout and test vector:
Bitcoin wiki 'Block hashing algorithm' (block 125552). Adders: Mano,
Digital Logic and Computer Design (Prentice-Hall), Sec. 4-3 and 5-2.
Usage: python3 nb_gate_count.py [samples, multiple of 64]
Requires Python 3 and NumPy (1.x or 2.x).

Gate-level double SHA-256 for Bitcoin mining, built only from two-input
XOR, AND and OR gates, exactly as counted in the N_b derivation:
  * ripple-carry adders (half adder at bit 0, 30 full adders, sum-only top bit)
  * Ch = g ^ (e & (f ^ g)), Maj = (a & b) | (c & (a ^ b))
  * Sigma/sigma as XORs of rotated/shifted copies (rotations are wiring)
The circuit is evaluated bit-sliced on S random nonces for a fixed header.

It reports
  1. correctness: digests against hashlib and against two real blocks;
  2. the plain gate count (every gate of the two compressions per nonce);
  3. the folded count (gates whose inputs vary with the nonce only);
  4. the early-exit count (only gates that feed H7 of the final digest);
  5. per-gate information lost, H(inputs) - H(output), and H(output),
     estimated from the samples (bits).
"""
import hashlib
import math
import sys
import numpy as np

S = int(sys.argv[1]) if len(sys.argv) > 1 else 8192     # samples (multiple of 64)
L = S // 64
rng = np.random.default_rng(20260926)
FULL = np.uint64(0xFFFFFFFFFFFFFFFF)


class Wire:
    __slots__ = ("c", "v", "id", "s")

    def __init__(self, c=None, v=None, nid=None, s=False):
        self.c, self.v, self.id, self.s = c, v, nid, s


def const(bit, structural=False):
    return Wire(c=int(bit), s=structural)


ZS = const(0, structural=True)          # structural zero (SHR fill)

# ------------------------------------------------------------ bookkeeping
node_type, node_in1, node_in2, node_tag = [], [], [], []
node_hin, node_hout = [], []
plain = {}                              # tag -> gate count before folding
nots = {}
TAG = ["none"]


def h2(p):
    return 0.0 if p <= 0.0 or p >= 1.0 else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


def hjoint(counts):
    return -sum((c / S) * math.log2(c / S) for c in counts if c > 0)


if hasattr(np, "bitwise_count"):          # NumPy 2.0 or later
    def pc(x):
        return int(np.bitwise_count(x).sum())
else:                                     # NumPy 1.x: 8-bit lookup table
    _POP8 = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)

    def pc(x):
        b = np.ascontiguousarray(x, dtype=np.uint64).view(np.uint8)
        return int(_POP8[b].sum(dtype=np.int64))


def new_node(kind, a, b, val, hin, hout):
    nid = len(node_type)
    node_type.append(kind)
    node_in1.append(a.id)
    node_in2.append(b.id if b is not None else -1)
    node_tag.append(TAG[0])
    node_hin.append(hin)
    node_hout.append(hout)
    return Wire(v=val, nid=nid)


def count_plain(a, b):
    if not (a.s or b.s):
        plain[TAG[0]] = plain.get(TAG[0], 0) + 1


def stats2(a, b, out_ones):
    n11 = pc(a.v & b.v)
    n10 = pc(a.v & ~b.v)
    n01 = pc(~a.v & b.v)
    n00 = S - n11 - n10 - n01
    return hjoint((n00, n01, n10, n11)), h2(out_ones / S)


def NOT(a):
    if a.c is not None:
        return const(1 - a.c)
    nots[TAG[0]] = nots.get(TAG[0], 0) + 1
    h = h2(pc(a.v) / S)
    return new_node("NOT", a, None, ~a.v, h, h)


def XOR(a, b):
    count_plain(a, b)
    if a.c is not None and b.c is not None:
        return const(a.c ^ b.c)
    if a.c is not None:
        a, b = b, a
    if b.c is not None:                  # one constant input
        return a if b.c == 0 else NOT(a)
    if a.id == b.id:
        return const(0)
    v = a.v ^ b.v
    hin, hout = stats2(a, b, pc(v))
    return new_node("XOR", a, b, v, hin, hout)


def AND(a, b):
    count_plain(a, b)
    if a.c is not None and b.c is not None:
        return const(a.c & b.c)
    if a.c is not None:
        a, b = b, a
    if b.c is not None:
        return const(0) if b.c == 0 else a
    if a.id == b.id:
        return a
    v = a.v & b.v
    hin, hout = stats2(a, b, pc(v))
    return new_node("AND", a, b, v, hin, hout)


def OR(a, b):
    count_plain(a, b)
    if a.c is not None and b.c is not None:
        return const(a.c | b.c)
    if a.c is not None:
        a, b = b, a
    if b.c is not None:
        return const(1) if b.c == 1 else a
    if a.id == b.id:
        return a
    v = a.v | b.v
    hin, hout = stats2(a, b, pc(v))
    return new_node("OR", a, b, v, hin, hout)



# ------------------------------------------------------------ cell-level accounting
cell_cost = {}      # cell kind -> [count, sum of H(in)-H(out)]
CELL_TAGS = {}      # tag prefix -> sum of cell costs


def hwires(ws):
    ws = [w for w in ws if w.c is None]
    if not ws:
        return 0.0
    k = len(ws)
    counts = []
    for pat in range(1 << k):
        m = np.full(L, FULL, dtype=np.uint64)
        for j, w in enumerate(ws):
            m &= w.v if (pat >> j) & 1 else ~w.v
        counts.append(pc(m))
    return hjoint(counts)


def record_cell(kind, ins, outs):
    if all(w.c is not None for w in ins):
        return
    cost = hwires(ins) - hwires(outs)
    rec = cell_cost.setdefault(kind, [0, 0.0])
    rec[0] += 1
    rec[1] += cost
    key = TAG[0].split(":")[0]
    CELL_TAGS[key] = CELL_TAGS.get(key, 0.0) + cost
    CELL_TAGS[key + ":n"] = CELL_TAGS.get(key + ":n", 0) + 1

# ------------------------------------------------------------ words (bit 0 = LSB)
def word_const(x):
    return [const((x >> i) & 1) for i in range(32)]


def is_const(w):
    return all(b.c is not None for b in w)


def rotr(w, n):
    return [w[(i + n) % 32] for i in range(32)]


def shr(w, n):
    return [w[i + n] if i + n < 32 else ZS for i in range(32)]


def xor3(x, y, z):
    out = []
    for i in range(32):
        w = XOR(XOR(x[i], y[i]), z[i])
        record_cell("XOR3", [x[i], y[i], z[i]], [w])
        out.append(w)
    return out


def Sig0(x): return xor3(rotr(x, 2), rotr(x, 13), rotr(x, 22))
def Sig1(x): return xor3(rotr(x, 6), rotr(x, 11), rotr(x, 25))
def sig0(x): return xor3(rotr(x, 7), rotr(x, 18), shr(x, 3))
def sig1(x): return xor3(rotr(x, 17), rotr(x, 19), shr(x, 10))


def Ch(e, f, g):
    out = []
    for i in range(32):
        w = XOR(g[i], AND(e[i], XOR(f[i], g[i])))
        record_cell("MUX", [e[i], f[i], g[i]], [w])
        out.append(w)
    return out


def Maj(a, b, c):
    out = []
    for i in range(32):
        w = OR(AND(a[i], b[i]), AND(c[i], XOR(a[i], b[i])))
        record_cell("MAJ", [a[i], b[i], c[i]], [w])
        out.append(w)
    return out


def add(x, y):
    s = [None] * 32
    s[0] = XOR(x[0], y[0])
    c = AND(x[0], y[0])
    record_cell("HA", [x[0], y[0]], [s[0], c])
    for i in range(1, 31):
        cin = c
        p = XOR(x[i], y[i])
        s[i] = XOR(p, c)
        c = OR(AND(x[i], y[i]), AND(c, p))
        record_cell("FA", [x[i], y[i], cin], [s[i], c])
    s[31] = XOR(XOR(x[31], y[31]), c)
    record_cell("XOR3", [x[31], y[31], c], [s[31]])
    return s


def add_many(terms):
    """Sum words mod 2^32, constant words first (maximises folding;
    the number of adders, len(terms) - 1, is the same in any order)."""
    terms = sorted(terms, key=lambda w: 0 if is_const(w) else 1)
    acc = terms[0]
    for t in terms[1:]:
        acc = add(acc, t)
    return acc


K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]
H0 = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]


def compress(state, block, tag):
    W = list(block)
    TAG[0] = tag + ":schedule"
    for t in range(16, 64):
        W.append(add_many([sig1(W[t - 2]), W[t - 7], sig0(W[t - 15]), W[t - 16]]))
    a, b, c, d, e, f, g, h = state
    for t in range(64):
        TAG[0] = tag + ":round:%02d" % t
        T1 = add_many([word_const(K[t]), W[t], h, Sig1(e), Ch(e, f, g)])
        T2 = add(Sig0(a), Maj(a, b, c))
        h, g, f = g, f, e
        e = add(d, T1)
        d, c, b = c, b, a
        a = add(T1, T2)
    TAG[0] = tag + ":final"
    return [add(state[j], v) for j, v in enumerate((a, b, c, d, e, f, g, h))]


# ------------------------------------------------------------ the header
# Bitcoin wiki 'Block hashing algorithm', block 125552 (serialized header)
HDR_HEX = ("01000000"
           "81cd02ab7e569e8bcd9317e2fe99f2de44d49ab2b8851ba4a308000000000000"
           "e320b6c2fffc8d750423db8b1eb942ae710e951ed797f7affc8892b0f1fc122b"
           "c7f5d74d" "f2b9441a" "42a14695")
WIKI_HASH = "00000000000000001e8d6829a8a21adc5d38d0a473b144b6765798e61f98bd1d"
hdr = bytes.fromhex(HDR_HEX)
true_nonce = int.from_bytes(hdr[76:80], "little")

# nonces: sample 0 is the real nonce, the rest uniform random 32-bit values
nonces = rng.integers(0, 2**32, size=S, dtype=np.uint64)
nonces[0] = true_nonce


def pack_bits(values, bit):
    bits = ((values >> np.uint64(bit)) & np.uint64(1)).astype(np.uint8)
    return np.packbits(bits.reshape(L, 64)[:, ::-1], axis=1).view(">u8").astype(np.uint64).ravel()


def unpack_word(w):
    out = np.zeros(S, dtype=np.uint64)
    for i, b in enumerate(w):
        if b.c is not None:
            if b.c:
                out |= np.uint64(1 << i)
        else:
            bits = np.unpackbits(b.v.astype(">u8").view(np.uint8).reshape(L, 8), axis=1)[:, ::-1]
            out |= bits.reshape(S).astype(np.uint64) << np.uint64(i)
    return out


# nonce input wires (primary inputs)
TAG[0] = "input"
nonce_bits = []
for i in range(32):
    v = pack_bits(nonces, i)
    nid = len(node_type)
    node_type.append("IN"); node_in1.append(-1); node_in2.append(-1); node_tag.append("input")
    node_hin.append(0.0); node_hout.append(0.0)
    nonce_bits.append(Wire(v=v, nid=nid))

# block 1: header bytes 0-63, constant for every nonce -> midstate
blk1 = [word_const(int.from_bytes(hdr[4 * j:4 * j + 4], "big")) for j in range(16)]
mid = compress([word_const(x) for x in H0], blk1, "midstate")

# block 2: header bytes 64-79 + padding (0x80, zeros, length 640)
w3 = [None] * 32                         # big-endian read of the little-endian nonce
for k in range(4):
    for j in range(8):
        w3[8 * (3 - k) + j] = nonce_bits[8 * k + j]
blk2 = ([word_const(int.from_bytes(hdr[64 + 4 * j:68 + 4 * j], "big")) for j in range(3)]
        + [w3, word_const(0x80000000)] + [word_const(0)] * 10 + [word_const(640)])
inner = compress(mid, blk2, "inner")

# outer SHA-256 of the 32-byte inner digest: one block
blk3 = inner + [word_const(0x80000000)] + [word_const(0)] * 6 + [word_const(256)]
outer = compress([word_const(x) for x in H0], blk3, "outer")

# ------------------------------------------------------------ 1. correctness
words = [unpack_word(w) for w in outer]
ok = 0
for s in range(S):
    h = hdr[:76] + int(nonces[s]).to_bytes(4, "little")
    ref = hashlib.sha256(hashlib.sha256(h).digest()).digest()
    got = b"".join(int(words[j][s]).to_bytes(4, "big") for j in range(8))
    ok += (ref == got)
real = b"".join(int(words[j][0]).to_bytes(4, "big") for j in range(8))
print(f"1. digests equal to hashlib: {ok} of {S}")
print(f"   block 125552, real nonce {true_nonce}: display hash {real[::-1].hex()}")
print(f"   matches the Bitcoin wiki value: {real[::-1].hex() == WIKI_HASH};"
      f" H7 of that digest = {int(words[7][0]):#010x}")

# ------------------------------------------------------------ 2. plain counts
def tot(d, key):
    return sum(v for k, v in d.items() if k.startswith(key))

print("\n2. Plain gate count (every two-input gate instantiated)")
for key in ("midstate", "inner", "outer"):
    print(f"   {key:9s}: {tot(plain, key):7d}   (schedule {tot(plain, key + ':schedule')},"
          f" rounds {tot(plain, key + ':round')}, final {tot(plain, key + ':final')})")
print(f"   per nonce (inner + outer) = {tot(plain, 'inner') + tot(plain, 'outer')}")

# ------------------------------------------------------------ 3. folded counts
typ = np.array(node_type)
tags = np.array(node_tag)
two = np.isin(typ, ["XOR", "AND", "OR"])
print("\n3. After constant folding (header fixed, only the nonce varies)")
for key in ("midstate", "inner", "outer"):
    m = two & np.char.startswith(tags.astype(str), key)
    print(f"   {key:9s}: {int(m.sum()):7d} two-input gates, {tot(nots, key)} inverters")
fold_total = int((two & ~np.char.startswith(tags.astype(str), "input")).sum())
print(f"   per nonce = {fold_total}")

# ------------------------------------------------------------ 4. early exit (cone of H7)
need = np.zeros(len(node_type), dtype=bool)
stack = [b.id for b in outer[7] if b.c is None]
while stack:
    n = stack.pop()
    if n < 0 or need[n]:
        continue
    need[n] = True
    stack.extend((node_in1[n], node_in2[n]))
cone = need & two
print("\n4. Early exit: gates that feed H7 of the final digest")
for key in ("inner:schedule", "inner:round", "inner:final", "outer:schedule", "outer:round", "outer:final"):
    m = cone & np.char.startswith(tags.astype(str), key)
    print(f"   {key:15s}: {int(m.sum()):7d}")
print(f"   per nonce = {int(cone.sum())}   (= {cone.sum() / 120448:.3f} compressions of 120,448)")

# ------------------------------------------------------------ 5. information per gate
hin = np.array(node_hin)
hout = np.array(node_hout)
loss = hin - hout
for label, m in (("all folded gates", two & ~np.char.startswith(tags.astype(str), "input")),
                 ("early-exit cone", cone)):
    n = int(m.sum())
    print(f"\n5. {label}: {n} gates")
    print(f"   information lost, sum of H(in)-H(out) = {loss[m].sum():9.0f} bits"
          f" = {loss[m].mean():.3f} bits per gate")
    print(f"   entropy of outputs, sum of H(out)     = {hout[m].sum():9.0f} bits"
          f" = {hout[m].mean():.3f} bits per gate")
    for k in ("XOR", "AND", "OR"):
        mk = m & (typ == k)
        print(f"   {k:3s}: {int(mk.sum()):6d} gates, lost {loss[mk].mean():.3f},"
              f" H(out) {hout[mk].mean():.3f} bits per gate")
print(f"\n(S = {S} samples; nodes created = {len(node_type)})")
print("\n6. Per round: folded gates / gates in the H7 cone (plain count is 1,430 per round)")
for comp in ("inner", "outer"):
    rows = []
    for t in list(range(0, 6)) + list(range(58, 64)):
        key = "%s:round:%02d" % (comp, t)
        m2 = two & (tags == key)
        rows.append("r%02d %5d/%5d" % (t, int(m2.sum()), int((cone & (tags == key)).sum())))
    print("   " + comp + ": " + ", ".join(rows))


print("\n7. Cell-level Landauer cost, sum of H(in)-H(out) over cells (bits)")
for k, (n, cst) in sorted(cell_cost.items()):
    print(f"   {k:4s}: {n:6d} cells, {cst:9.0f} bits, {cst / n:.3f} bits per cell")
for key in ("midstate", "inner", "outer"):
    print(f"   {key:9s}: {CELL_TAGS.get(key + ':n', 0):6d} cells, {CELL_TAGS.get(key, 0.0):9.0f} bits")
print(f"   per nonce (inner + outer): {CELL_TAGS.get('inner', 0) + CELL_TAGS.get('outer', 0):.0f} bits")
