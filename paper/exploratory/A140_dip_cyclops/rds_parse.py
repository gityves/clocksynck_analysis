"""Minimal pure-python reader of R serialization (XDR, v2/v3) for untrusted .rds.gz files.
Nothing is executed (no R). S4 objects and lists become dicts / Obj; big numeric vectors are
summarised (BigVec) instead of kept. Run with python -I."""
import gzip, struct, sys
import numpy as np

BIG = 3_000_000


class Obj(dict):
    """R object: payload in .data, attributes in dict items."""
    data = None
    kind = None


class BigVec:
    def __init__(self, kind, n, head, integer_valued=None, vmin=None, vmax=None):
        self.kind, self.n, self.head, self.integer_valued, self.vmin, self.vmax = kind, n, head, integer_valued, vmin, vmax

    def __repr__(self):
        return f"BigVec({self.kind},n={self.n},int={self.integer_valued},min={self.vmin},max={self.vmax})"


class Reader:
    def __init__(self, f):
        self.f = f
        self.refs = []

    def rd(self, n):
        b = self.f.read(n)
        while len(b) < n:
            c = self.f.read(n - len(b))
            if not c:
                raise EOFError
            b += c
        return b

    def i(self):
        return struct.unpack(">i", self.rd(4))[0]

    def length(self):
        n = self.i()
        if n == -1:
            a, b = self.i(), self.i()
            n = ((a & 0xFFFFFFFF) << 32) + (b & 0xFFFFFFFF)
        return n

    def header(self):
        assert self.rd(2) == b"X\n"
        v = self.i(); self.i(); self.i()
        if v == 3:
            self.rd(self.i())
        return v

    def attrs(self):
        return self.item()

    def pairlist(self, flags):
        out = []
        while True:
            t = flags & 0xFF
            if t == 254:
                return out
            attr = self.item() if flags & 0x200 else None
            tag = self.item() if flags & 0x400 else None
            car = self.item()
            out.append((tag, car))
            flags = self.i()
            if flags & 0xFF not in (2, 6, 3, 17, 254):
                # tail is not a pairlist cell; read as item and stop
                out.append((None, self.item_flags(flags)))
                return out

    def item(self):
        return self.item_flags(self.i())

    def num(self, t, n):
        dt, w = {10: (">i4", 4), 13: (">i4", 4), 14: (">f8", 8), 15: (">f8", 8)}[t]
        if t == 15:
            n *= 2
        if n > BIG:
            first = np.frombuffer(self.rd(w * min(n, 1000)), dtype=dt)
            head = first.copy()
            done = len(first)
            iv, mn, mx = True, float(first.min()), float(first.max())
            if t == 14:
                iv = bool(np.all(first == np.round(first)))
            while done < n:
                k = min(n - done, 4_000_000)
                a = np.frombuffer(self.rd(w * k), dtype=dt)
                done += k
                mn, mx = min(mn, float(np.nanmin(a))), max(mx, float(np.nanmax(a)))
                if t == 14 and iv:
                    iv = bool(np.all(a == np.round(a)))
            return BigVec({10: "lgl", 13: "int", 14: "real", 15: "cplx"}[t], n, head, iv, mn, mx)
        return np.frombuffer(self.rd(w * n), dtype=dt).copy()

    def item_flags(self, flags):
        t = flags & 0xFF
        hasattr_, hastag, isobj = flags & 0x200, flags & 0x400, flags & 0x100
        if t in (254,):
            return None
        if t in (242, 241, 253, 251, 250, 252):
            return {242: "<emptyenv>", 241: "<baseenv>", 253: "<globalenv>", 251: "<missing>",
                    250: "<unbound>", 252: "<basens>"}[t]
        if t == 255:
            idx = flags >> 8
            if idx == 0:
                idx = self.i()
            return self.refs[idx - 1]
        if t in (247, 248, 249):
            self.i()
            n = self.i()
            v = [self.item() for _ in range(n)]
            self.refs.append(("pkg", v))
            return ("pkg", v)
        if t == 1:
            s = self.item()
            self.refs.append(s)
            return s
        if t == 9:
            n = self.i()
            return None if n == -1 else self.rd(n).decode("utf8", "replace")
        if t == 4:
            self.i()
            k = len(self.refs)
            self.refs.append("<env>")
            for _ in range(4):
                self.item()
            return "<env>"
        if t in (2, 6, 3, 17, 5):
            pl = self.pairlist(flags) if t in (2, 17) else self._lang(flags)
            if t == 2:
                o = Obj(); o.kind = "pairlist"
                for tag, car in pl:
                    o[tag if isinstance(tag, str) else len(o)] = car
                return o
            return "<lang/closure>"
        if t in (10, 13, 14, 15):
            n = self.length()
            data = self.num(t, n)
        elif t == 16:
            n = self.length()
            data = [self.item() for _ in range(n)]
        elif t in (19, 20):
            n = self.length()
            data = [self.item() for _ in range(n)]
        elif t == 24:
            n = self.length()
            data = self.rd(n)
        elif t == 25:
            data = None
        elif t == 238:
            info = self.item()
            state = self.item()
            attr = self.item()
            cls = info[0][1] if isinstance(info, Obj) is False and info else None
            return self.altrep(info, state, attr)
        else:
            raise ValueError(f"unsupported SEXP type {t}")
        o = Obj(); o.data = data; o.kind = t
        if hasattr_:
            a = self.item()
            if isinstance(a, Obj):
                for k, v in a.items():
                    o[k] = v
        return o

    def _lang(self, flags):
        # closure / lang / promise: consume attr, tag(s), car, cdr structurally
        if flags & 0x200:
            self.item()
        if flags & 0x400:
            self.item()
        self.item(); self.item()
        return []

    def altrep(self, info, state, attr):
        cls = None
        if isinstance(info, Obj):
            vals = list(info.values())
            cls = vals[0] if vals else None
        o = Obj(); o.kind = "altrep:" + str(cls)
        if cls == "compact_intseq" or cls == "compact_realseq":
            n, start, step = state.data[:3]
            o.data = (start + step * np.arange(int(n))).astype(np.int64 if cls == "compact_intseq" else float)
        elif cls and cls.startswith("wrap_"):
            o = state.data[0] if isinstance(state, Obj) else state
        elif cls == "deferred_string":
            arg = state[0] if isinstance(state, Obj) and 0 in state else state
            o.data = ("deferred_string", getattr(state, "data", None))
        else:
            o.data = state
        if isinstance(attr, Obj) and isinstance(o, Obj):
            for k, v in attr.items():
                o[k] = v
        return o


def read_rds(path):
    with gzip.open(path, "rb") as f0:
        f = f0
        if f0.peek(2)[:2] == b"\x1f\x8b":  # saveRDS gzip inside the GEO .gz
            f = gzip.GzipFile(fileobj=f0, mode="rb")
        r = Reader(f)
        r.header()
        return r.item()


def strs(o):
    return None if o is None else o.data


if __name__ == "__main__":
    obj = read_rds(sys.argv[1])
    print(type(obj), list(obj.keys())[:30])
