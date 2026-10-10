"""Locate (or build) the hypersample shared library and declare its C ABI for ctypes.

Search order: $HYPERSAMPLE_LIB; the repository's build/ directory (from `make`); otherwise
the C sources are compiled once with the system C compiler into a cache directory
($HYPERSAMPLE_CACHE, else ~/.cache/hypersample), keyed by a hash of the sources.
"""
from __future__ import annotations

import ctypes
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ABI_VERSION = 1
_REPO = Path(__file__).resolve().parents[2]
_SUFFIX = ".dylib" if sys.platform == "darwin" else (".dll" if os.name == "nt" else ".so")

i32, i64, u64, dbl = ctypes.c_int32, ctypes.c_int64, ctypes.c_uint64, ctypes.c_double
P32, P64, PD = ctypes.POINTER(i32), ctypes.POINTER(i64), ctypes.POINTER(dbl)

_SIGNATURES = {
    "hs_abi_version": ([], ctypes.c_int),
    "hs_version_string": ([], ctypes.c_char_p),
    "hs_status_string": ([ctypes.c_int], ctypes.c_char_p),
    "hs_config_degrees": ([i32, i32, P32, u64, i64, P32, P64], ctypes.c_int),
    "hs_config_regular": ([i32, i32, i32, u64, i64, P32, P64], ctypes.c_int),
    "hs_gnp": ([i32, i32, dbl, u64, P32, i64, P64], ctypes.c_int),
    "hs_gnm": ([i32, i32, i64, u64, P32], ctypes.c_int),
    "hs_sequential_degrees": ([i32, i32, P32, u64, P32, PD], ctypes.c_int),
    "hs_switch_chain": ([i32, i32, i64, P32, i64, u64, P64], ctypes.c_int),
    "hs_is_simple": ([i32, i32, i64, P32], ctypes.c_int),
    "hs_degrees": ([i32, i32, i64, P32, P32], ctypes.c_int),
    "hs_binomial": ([i32, i32], i64),
}


def _sources():
    return sorted((_REPO / "src").glob("*.c")) + [_REPO / "src" / "hs_internal.h", _REPO / "include" / "hypersample.h"]


def _compile() -> Path:
    srcs = _sources()
    if not all(p.exists() for p in srcs):
        raise OSError("hypersample C sources not found; set HYPERSAMPLE_LIB or run `make`")
    tag = hashlib.sha256(b"".join(p.read_bytes() for p in srcs)).hexdigest()[:16]
    cache = Path(os.environ.get("HYPERSAMPLE_CACHE") or Path.home() / ".cache" / "hypersample")
    target = cache / f"libhypersample_{tag}{_SUFFIX}"
    if target.exists():
        return target
    cc = os.environ.get("CC") or shutil.which("cc") or shutil.which("clang") or shutil.which("gcc")
    if cc is None:
        raise OSError("no C compiler found; install one or set HYPERSAMPLE_LIB")
    cache.mkdir(parents=True, exist_ok=True)
    flags = ["-dynamiclib"] if sys.platform == "darwin" else ["-shared"]
    with tempfile.TemporaryDirectory(dir=cache) as tmp:            # same filesystem for os.replace
        out = Path(tmp) / target.name
        cmd = [cc, "-O2", "-std=c11", "-fPIC", *flags, "-I", str(_REPO / "include"),
               *[str(p) for p in srcs if p.suffix == ".c"], "-o", str(out), "-lm"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode:
            raise OSError("compiling hypersample failed:\n" + res.stderr)
        os.replace(out, target)
    return target


def load() -> ctypes.CDLL:
    path = os.environ.get("HYPERSAMPLE_LIB")
    if not path:
        built = _REPO / "build" / f"libhypersample{_SUFFIX}"
        path = str(built) if built.exists() else str(_compile())
    lib = ctypes.CDLL(path)
    for name, (args, res) in _SIGNATURES.items():
        fn = getattr(lib, name)
        fn.argtypes, fn.restype = args, res
    if lib.hs_abi_version() != ABI_VERSION:
        raise OSError(f"hypersample ABI {lib.hs_abi_version()} at {path}, binding expects {ABI_VERSION}")
    return lib
