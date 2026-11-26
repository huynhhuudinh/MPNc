"""Minimal, dependency-free reader for R's RDX2/RDX3 XDR save files.

It supports the object types that occur in AHMbook's ``cswa.rda`` (pairlists,
symbols, character/integer/double/logical/string vectors, generic vectors and
attributes). It is used in two ways:

* as a fallback in ``run_application.load_data`` when the ``rdata`` package is
  not installed; and
* by ``test_application.py`` to check the normalized data independently of
  ``rdata``/``xarray``.

Unsupported types raise ``ValueError`` instead of being skipped silently.
"""
from __future__ import annotations

import bz2
import gzip
import lzma
import struct
from pathlib import Path

import numpy as np

NILVALUE, EMPTYENV, BASEENV, GLOBALENV, UNBOUND, MISSINGARG, BASENAMESPACE = (
    254, 242, 241, 253, 252, 251, 247)
REFSXP = 255
NA_INTEGER = -2147483648


class RObject:
    """An R value with optional attributes."""

    def __init__(self, rtype, value, attributes=None):
        self.rtype = rtype
        self.value = value
        self.attributes = attributes or {}

    def attr(self, name, default=None):
        return self.attributes.get(name, default)

    def __repr__(self):
        return f"RObject({self.rtype}, attrs={list(self.attributes)})"


class _Reader:
    def __init__(self, data: bytes):
        self.data = data
        self.pos = 0
        self.refs = []

    def _int(self):
        (v,) = struct.unpack_from(">i", self.data, self.pos)
        self.pos += 4
        return v

    def _length(self):
        n = self._int()
        if n == -1:
            hi, lo = self._int(), self._int()
            n = (hi << 32) + lo
        return n

    def _attributes(self, has_attr):
        if not has_attr:
            return {}
        pairs = self.item()
        if not isinstance(pairs, list):
            raise ValueError("Attributes must be a pairlist")
        return {name: value for name, value in pairs}

    def item(self):
        flags = self._int()
        rtype = flags & 0xFF
        has_attr = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))
        if rtype == NILVALUE:
            return None
        if rtype in (EMPTYENV, BASEENV, GLOBALENV, UNBOUND, MISSINGARG, BASENAMESPACE):
            return RObject("special", rtype)
        if rtype == REFSXP:
            index = flags >> 8
            if index == 0:
                index = self._int()
            return self.refs[index - 1]
        if rtype == 1:  # SYMSXP
            name = self.item()
            self.refs.append(name)
            return name
        if rtype == 2:  # LISTSXP (pairlist): returned as [(tag, value), ...]
            out = []
            while True:
                attributes = self._attributes(has_attr)
                tag = self.item() if has_tag else None
                out.append((tag, self.item()))
                del attributes  # pairlist-node attributes are not used here
                flags = self._int()
                rtype = flags & 0xFF
                has_attr = bool(flags & (1 << 9))
                has_tag = bool(flags & (1 << 10))
                if rtype == NILVALUE:
                    return out
                if rtype != 2:
                    raise ValueError("Unsupported pairlist continuation")
        if rtype == 9:  # CHARSXP
            n = self._int()
            if n == -1:
                return None
            s = self.data[self.pos:self.pos + n].decode("utf-8")
            self.pos += n
            return s
        if rtype in (10, 13):  # LGLSXP, INTSXP
            n = self._length()
            v = np.frombuffer(self.data, ">i4", n, self.pos).astype(np.int64)
            self.pos += 4 * n
            return RObject("logical" if rtype == 10 else "integer", v,
                           self._attributes(has_attr))
        if rtype == 14:  # REALSXP
            n = self._length()
            v = np.frombuffer(self.data, ">f8", n, self.pos).astype(float)
            self.pos += 8 * n
            return RObject("double", v, self._attributes(has_attr))
        if rtype == 16:  # STRSXP
            n = self._length()
            v = [self.item() for _ in range(n)]
            return RObject("character", v, self._attributes(has_attr))
        if rtype == 19:  # VECSXP
            n = self._length()
            v = [self.item() for _ in range(n)]
            return RObject("list", v, self._attributes(has_attr))
        raise ValueError(f"Unsupported R SEXP type {rtype}")


def read_rda(path):
    """Return ``{name: RObject}`` for every object saved in an .rda file."""
    raw = Path(path).read_bytes()
    for decompress in (gzip.decompress, bz2.decompress, lzma.decompress):
        try:
            raw = decompress(raw)
            break
        except (OSError, ValueError, lzma.LZMAError):
            continue
    if raw[:5] not in (b"RDX2\n", b"RDX3\n") or raw[5:7] != b"X\n":
        raise ValueError("Only XDR RDX2/RDX3 files are supported")
    reader = _Reader(raw)
    reader.pos = 7
    version = reader._int()
    reader._int()  # writer R version
    reader._int()  # minimal reader R version
    if version == 3:
        n = reader._int()
        reader.pos += n  # native encoding name
    objects = reader.item()
    return {name: value for name, value in objects}


def _vector(obj):
    """Numeric vector with R missing values converted to NaN."""
    if obj.rtype in ("integer", "logical"):
        values = obj.value.astype(float)
        values[obj.value == NA_INTEGER] = np.nan
        levels = obj.attr("levels")
        if levels is not None:  # factor -> labels
            return np.array([levels.value[int(i) - 1] if np.isfinite(i) else None
                             for i in values], dtype=object)
        return values
    if obj.rtype == "double":
        return obj.value.copy()
    if obj.rtype == "character":
        return np.array(obj.value, dtype=object)
    raise ValueError(f"Not a vector: {obj.rtype}")


def data_frame(obj):
    """Convert an R data.frame to (columns dict, row names list)."""
    names = obj.attr("names").value
    columns = {name: _vector(col) for name, col in zip(names, obj.value)}
    rows = obj.attr("row.names")
    if rows is not None and rows.rtype == "character":
        row_names = list(rows.value)
    else:
        row_names = [str(i + 1) for i in range(len(next(iter(columns.values()))))]
    return columns, row_names


def read_cswa(path):
    """Decode AHMbook ``cswa`` into plain NumPy/pandas structures.

    Returns ``counts`` in R's dimension order (point x interval x occasion),
    the three dimnames lists, the covariate data frame indexed by point id and
    the spot-mapping data frame indexed by patch id.
    """
    import pandas as pd

    cswa = read_rda(path)["cswa"]
    parts = dict(zip(cswa.attr("names").value, cswa.value))
    counts = parts["counts"]
    dims = tuple(int(v) for v in counts.attr("dim").value)
    dimnames = [list(d.value) for d in counts.attr("dimnames").value]
    array = _vector(counts).reshape(dims, order="F")  # R arrays are column-major
    cov_cols, cov_rows = data_frame(parts["covs"])
    spot_cols, spot_rows = data_frame(parts["spotMaps"])
    covs = pd.DataFrame(cov_cols, index=cov_rows)
    spot = pd.DataFrame(spot_cols, index=spot_rows)
    return array, dimnames, covs, spot
