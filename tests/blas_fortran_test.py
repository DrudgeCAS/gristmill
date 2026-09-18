"""Tests for the BLAS-based Fortran printer."""

import os
import re
import shutil
import subprocess

import numpy as np
import pytest
from drudge import PartHoleDrudge, Range, Perm, NEG, Drudge
from sympy import IndexedBase, Rational, Symbol, symbols

from gristmill import (
    optimize,
    verify_eval_seq,
    FortranPrinter,
    BlasFortranPrinter,
    EinsumPrinter,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_REF_DGEMM = os.path.join(_HERE, "ref_dgemm.f90")

needs_gfortran = pytest.mark.skipif(
    shutil.which("gfortran") is None, reason="gfortran is not available"
)


@pytest.fixture(scope="module")
def offset_parthole(spark_ctx):
    """The particle-hole drudge with the virtual range after the occupied."""

    no = Symbol("no")
    nv = Symbol("nv")
    dr = PartHoleDrudge(
        spark_ctx,
        part_orb=(Range("V", no, no + nv), PartHoleDrudge.DEFAULT_PART_DUMMS),
    )
    return dr


@pytest.fixture(scope="module")
def ccd_doubles(offset_parthole):
    """The CCD doubles amplitude equation with symmetric left-hand side."""

    dr = offset_parthole
    p = dr.names
    a, b = p.V_dumms[:2]
    i, j = p.O_dumms[:2]
    c_ = p.c_
    c_dag = p.c_dag

    t = IndexedBase("t")
    dr.set_dbbar_base(t, 2)
    t2 = dr.einst(
        Rational(1, 4) * t[a, b, i, j] * c_dag[a] * c_dag[b] * c_[j] * c_[i]
    )

    curr = dr.ham
    h_bar = dr.ham
    for order in range(2):
        curr = (curr | t2).simplify() / (order + 1)
        h_bar += curr
        continue
    h_bar = h_bar.simplify()

    proj = c_dag[i] * c_dag[j] * c_[b] * c_[a]
    r = IndexedBase("r")
    dr.set_symm(r, Perm([1, 0, 2, 3], NEG), Perm([0, 1, 3, 2], NEG), valence=4)
    return dr.define(r[a, b, i, j], (proj * h_bar).eval_fermi_vev().simplify())


_PROGRAM = """
program main
  implicit none
{params}
{inputs}
{outputs}
  integer :: unit_
{decls}

{reads}
{evals}
{writes}
end program main
"""


def run_fortran(
    tmp_path, printer, eval_seq, params, inputs, outputs, tag="run", flags=()
):
    """Compile and run the generated Fortran code, returning the outputs.

    Parameters
    ----------

    params
        Mapping from integer parameter names to their values.

    inputs
        Mapping from input array names to pairs of the NumPy array and the
        Fortran bounds declaration.

    outputs
        Mapping from output array names to pairs of the shape and the Fortran
        bounds declaration.

    """

    decls, evals = printer.doprint(eval_seq, separate_decls=True)

    # All the indices need to be declared as loop variables.
    index_names = set()
    for def_ in eval_seq:
        index_names.update(str(i) for i, _ in def_.exts)
        for term in def_.rhs_terms:
            index_names.update(str(i) for i, _ in term.sums)
            continue
        continue
    decls = "  integer :: {}\n".format(", ".join(sorted(index_names))) + decls

    param_lines = []
    for k, v in params.items():
        param_lines.append("  integer, parameter :: {} = {}".format(k, v))
        continue

    input_lines = []
    read_lines = []
    for name, (arr, bounds) in inputs.items():
        input_lines.append(
            "  real(kind=8), dimension({}) :: {}".format(bounds, name)
        )
        fname = os.path.join(tmp_path, "{}_{}.bin".format(tag, name))
        np.asarray(arr, dtype=np.float64).flatten(order="F").tofile(fname)
        read_lines.append(
            "  open(newunit=unit_, file='{}', access='stream', "
            "form='unformatted', status='old')\n  read(unit_) {}\n"
            "  close(unit_)".format(fname, name)
        )
        continue

    output_lines = []
    write_lines = []
    out_files = {}
    for name, (shape, bounds) in outputs.items():
        output_lines.append(
            "  real(kind=8), dimension({}) :: {}".format(bounds, name)
        )
        fname = os.path.join(tmp_path, "{}_{}_out.bin".format(tag, name))
        out_files[name] = (fname, shape)
        write_lines.append(
            "  open(newunit=unit_, file='{}', access='stream', "
            "form='unformatted', status='replace')\n  write(unit_) {}\n"
            "  close(unit_)".format(fname, name)
        )
        continue

    code = _PROGRAM.format(
        params="\n".join(param_lines),
        inputs="\n".join(input_lines),
        outputs="\n".join(output_lines),
        decls=decls,
        reads="\n".join(read_lines),
        evals=evals,
        writes="\n".join(write_lines),
    )
    src = os.path.join(tmp_path, "{}.f90".format(tag))
    with open(src, "w") as fp:
        fp.write(code)

    exe = os.path.join(tmp_path, "{}.exe".format(tag))
    comp = subprocess.run(
        [
            "gfortran",
            "-O1",
            "-Wall",
            "-ffree-line-length-none",
            *flags,
            "-o",
            exe,
            src,
            _REF_DGEMM,
        ],
        capture_output=True,
        text=True,
    )
    assert comp.returncode == 0, comp.stderr + "\n" + code
    run = subprocess.run([exe], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr

    res = {}
    for name, (fname, shape) in out_files.items():
        res[name] = np.fromfile(fname, dtype=np.float64).reshape(
            shape, order="F"
        )
        continue

    return res, code


def _dbbar_random(rng, shape):
    """Form a random array with the double-bar symmetry."""

    arr = rng.random(shape)
    arr = arr - arr.transpose(1, 0, 2, 3)
    arr = arr - arr.transpose(0, 1, 3, 2)
    if shape[0] == shape[2]:
        arr = arr + arr.transpose(2, 3, 0, 1)
    return arr


@pytest.fixture(scope="module")
def matrix_drudge(spark_ctx):
    """A drudge with two ranges starting from zero for GEMM pattern tests."""

    dr = Drudge(spark_ctx)
    m = Symbol("m")
    n = Symbol("n")
    r = Range("R", 0, m)
    s = Range("S", 0, n)
    dr.set_dumms(r, symbols("i j k l"))
    dr.set_dumms(s, symbols("a b c d"))
    dr.add_resolver_for_dumms()
    return dr


@needs_gfortran
def test_gemm_patterns(matrix_drudge, tmp_path):
    """Test the GEMM planning on representative contraction patterns.

    The patterns cover direct GEMM without copies, copies of interleaved
    operands, mismatched order of the summed indices, permuted accumulation of
    the result, contraction to a scalar, outer products, and the fallback to
    loops for terms that are not binary contractions.
    """

    dr = matrix_drudge
    p = dr.names
    i, j, k, l = p.R_dumms[:4]
    a, b, c, d = p.S_dumms[:4]
    x = IndexedBase("x")  # Over R, R, R, R.
    y = IndexedBase("y")  # Over R, R, S, S.
    w = IndexedBase("w")  # Over S, S, R, R.
    z = IndexedBase("z")  # Over R.
    v = IndexedBase("v")  # Over S.

    defs = [
        # Direct: summed indices trailing in x and leading in y.
        dr.define_einst(
            IndexedBase("d1")[i, j, a, b], x[i, j, k, l] * y[k, l, a, b]
        ),
        # Both operands transposed, still direct.
        dr.define_einst(
            IndexedBase("d2")[i, j, a, b], x[k, l, i, j] * w[a, b, k, l]
        ),
        # Interleaved operand needs a copy.
        dr.define_einst(
            IndexedBase("d3")[i, j, a, b], x[i, k, j, l] * y[k, l, a, b]
        ),
        # Summed indices in different orders.
        dr.define_einst(
            IndexedBase("d4")[i, j, a, b], x[i, j, k, l] * y[l, k, a, b]
        ),
        # External indices interleaved in the target.
        dr.define_einst(
            IndexedBase("d5")[i, a, j, b], x[i, j, k, l] * y[k, l, a, b]
        ),
        # Contraction to a scalar and an outer product.
        dr.define_einst(Symbol("e"), y[i, j, a, b] * w[a, b, i, j]),
        dr.define_einst(IndexedBase("d6")[i, a], z[i] * v[a]),
        # Fallbacks: repeated index, and a power from the unary extraction.
        dr.define_einst(IndexedBase("d7")[i, j], 2 * y[j, i, a, a]),
        dr.define_einst(IndexedBase("d8")[i, a], z[i] * v[a] * v[a]),
    ]

    printer = BlasFortranPrinter(
        openmp=False, default_type="real(kind=8)", heap_interm=False
    )
    code = printer.doprint(defs)
    assert code.count("call dgemm(") == 7
    # Copies are made for the interleaved operand and for the mismatched
    # orders of the summed indices in d4 and e.
    assert len(re.findall(r"(?<!de)allocate\(gm_1\(", code)) == 1
    assert len(re.findall(r"(?<!de)allocate\(gm_2\(", code)) == 2
    # Temporaries for the permuted target and the scalar target.
    assert len(re.findall(r"(?<!de)allocate\(gm_c\(", code)) == 1
    assert "gm_c(1)" in code
    assert code.count("end block") == 7
    # Each term is preceded by its Einstein-notation comment.
    assert "! d1(i, j, a, b) += x(i, j, k, l) * y(k, l, a, b)" in code
    assert "! d7(i, j) += 2 * y(j, i, a, a)" in code
    plain = BlasFortranPrinter(
        openmp=False, default_type="real(kind=8)", comments=False
    ).doprint(defs)
    assert "! d1(i, j, a, b)" not in plain

    m, n = 3, 2
    rng = np.random.default_rng(3)
    arrays = {
        "x": rng.random((m, m, m, m)),
        "y": rng.random((m, m, n, n)),
        "w": rng.random((n, n, m, m)),
        "z": rng.random(m),
        "v": rng.random(n),
    }
    inputs = {
        "x": (arrays["x"], "m, m, m, m"),
        "y": (arrays["y"], "m, m, n, n"),
        "w": (arrays["w"], "n, n, m, m"),
        "z": (arrays["z"], "m"),
        "v": (arrays["v"], "n"),
    }
    outputs = {
        "d1": ((m, m, n, n), "m, m, n, n"),
        "d2": ((m, m, n, n), "m, m, n, n"),
        "d3": ((m, m, n, n), "m, m, n, n"),
        "d4": ((m, m, n, n), "m, m, n, n"),
        "d5": ((m, n, m, n), "m, n, m, n"),
        "e": ((1,), "1"),
        "d6": ((m, n), "m, n"),
        "d7": ((m, m), "m, m"),
        "d8": ((m, n), "m, n"),
    }
    res, _ = run_fortran(
        tmp_path, printer, defs, {"m": m, "n": n}, inputs, outputs, tag="pat"
    )

    xa, ya, wa = arrays["x"], arrays["y"], arrays["w"]
    za, va = arrays["z"], arrays["v"]
    assert np.allclose(res["d1"], np.einsum("ijkl,klab->ijab", xa, ya))
    assert np.allclose(res["d2"], np.einsum("klij,abkl->ijab", xa, wa))
    assert np.allclose(res["d3"], np.einsum("ikjl,klab->ijab", xa, ya))
    assert np.allclose(res["d4"], np.einsum("ijkl,lkab->ijab", xa, ya))
    assert np.allclose(res["d5"], np.einsum("ijkl,klab->iajb", xa, ya))
    assert np.allclose(res["e"][0], np.einsum("ijab,abij->", ya, wa))
    assert np.allclose(res["d6"], np.outer(za, va))
    assert np.allclose(res["d7"], 2 * np.einsum("jiaa->ij", ya))
    assert np.allclose(res["d8"], np.outer(za, va ** 2))


@needs_gfortran
@pytest.mark.parametrize("m,n", [(4, 2), (1, 2), (4, 0)])
def test_composite_and_empty_dimensions(spark_ctx, tmp_path, m, n):
    """Products of offset range sizes and empty contractions are safe."""
    dr = Drudge(spark_ctx)
    ms, ns, hs = symbols("m n h")
    i, j, k, l = symbols("i j k l")
    a, b = symbols("a b")
    dr.set_dumms(Range("R", hs, ms), (i, j, k, l))
    dr.set_dumms(Range("S", 0, ns), (a, b))
    dr.add_resolver_for_dumms()
    x, y, z = map(IndexedBase, ("x", "y", "z"))
    defs = [dr.define_einst(z[i, a], x[i, k, l] * y[l, k, a])]
    rng = np.random.default_rng(12)
    xa = rng.random((m - 1, m - 1, m - 1))
    ya = rng.random((m - 1, m - 1, n))
    printer = BlasFortranPrinter(
        openmp=True, explicit_bounds=True, default_type="real(kind=8)"
    )
    res, _ = run_fortran(
        tmp_path,
        printer,
        defs,
        {"m": m, "n": n, "h": 1},
        {"x": (xa, "h+1:m,h+1:m,h+1:m"), "y": (ya, "h+1:m,h+1:m,n")},
        {"z": ((m - 1, n), "h+1:m,n")},
        flags=("-fcheck=all", "-fopenmp"),
    )
    assert np.allclose(res["z"], np.einsum("ikl,lka->ia", xa, ya))


@needs_gfortran
def test_gemm_packs_smaller_operand(matrix_drudge, tmp_path):
    """Prefer packing amplitudes over a much larger four-virtual ERI."""
    dr = matrix_drudge
    i, j = dr.names.R_dumms[:2]
    a, b, c, d = dr.names.S_dumms[:4]
    aa, zz, out = map(IndexedBase, ("aa", "zz", "out"))
    defs = [dr.define_einst(out[i, j, a, b], aa[d, c, i, j] * zz[a, b, c, d])]
    printer = BlasFortranPrinter(
        openmp=False,
        default_type="real(kind=8)",
        size_substs={Symbol("m"): 2, Symbol("n"): 6},
    )
    code = printer.doprint(defs)
    assert re.search(r"gm_[12]\([^\n]+ = aa\(", code)
    assert not re.search(r"gm_[12]\([^\n]+ = zz\(", code)
    rng = np.random.default_rng(14)
    small = rng.random((6, 6, 2, 2))
    large = rng.random((6, 6, 6, 6))
    res, _ = run_fortran(
        tmp_path,
        printer,
        defs,
        {"m": 2, "n": 6},
        {"aa": (small, "n,n,m,m"), "zz": (large, "n,n,n,n")},
        {"out": ((2, 2, 6, 6), "m,m,n,n")},
        flags=("-fcheck=all",),
    )
    assert np.allclose(res["out"], np.einsum("dcij,abcd->ijab", small, large))


@needs_gfortran
def test_ccd_doubles_blas(ccd_doubles, tmp_path):
    """Test the BLAS Fortran code for the CCD doubles equation.

    The two-body interaction and the Fock matrix are stored over the whole
    orbital space, so they are always copied into block temporaries.  The
    result is compared with the NumPy code with input blocks and with the
    naive Fortran loops.
    """

    tdef = ccd_doubles
    dr = tdef.drudge
    p = dr.names
    u = dr.two_body
    f = dr.fock
    no_, nv_ = p.no, p.nv

    eval_seq = optimize(
        [tdef], substs={nv_: 10 * no_}, interm_fmt="tau{}", lhs_symm=True
    )
    assert verify_eval_seq(eval_seq, [tdef], simplify=True)

    no, nv = 3, 4
    n = no + nv
    rng = np.random.default_rng(11)
    u_arr = _dbbar_random(rng, (n, n, n, n))
    f_arr = rng.random((n, n))
    t_arr = _dbbar_random(rng, (nv, nv, no, no))

    # Reference from the NumPy block code.
    env = {
        "einsum": np.einsum,
        "no": no,
        "nv": nv,
        "u": u_arr,
        "f": f_arr,
        "t": t_arr,
    }
    exec(EinsumPrinter(blocks=[u, f], base_indent=0).doprint(eval_seq), env)
    r_ref = env["r"]

    params = {"no": no, "nv": nv, "n": n}
    inputs = {
        "u": (u_arr, "n, n, n, n"),
        "f": (f_arr, "n, n"),
        "t": (t_arr, "no + 1:n, no + 1:n, no, no"),
    }
    outputs = {"r": ((nv, nv, no, no), "no + 1:n, no + 1:n, no, no")}

    blas_printer = BlasFortranPrinter(
        openmp=False,
        default_type="real(kind=8)",
        explicit_bounds=True,
        copy_inputs=[u, f],
    )
    res, code = run_fortran(
        tmp_path, blas_printer, eval_seq, params, inputs, outputs, tag="blas"
    )
    assert np.allclose(res["r"], r_ref)

    # All contractions go through GEMM; only the single-factor terms and the
    # assembly are loops.
    n_contr = sum(
        1 for def_ in eval_seq for term in def_.rhs_terms if len(term.sums) > 0
    )
    assert code.count("call dgemm(") == n_contr
    # Temporaries are all released.
    n_alloc = len(re.findall(r"(?<!de)allocate\(gm_", code))
    assert n_alloc == code.count("deallocate(gm_")

    # With OpenMP: loop nests get their own parallel do, GEMM calls do not.
    omp_printer = BlasFortranPrinter(
        openmp=True,
        default_type="real(kind=8)",
        explicit_bounds=True,
        copy_inputs=[u, f],
    )
    res_omp, omp_code = run_fortran(
        tmp_path,
        omp_printer,
        eval_seq,
        params,
        inputs,
        outputs,
        tag="omp",
        flags=("-fopenmp",),
    )
    assert np.allclose(res_omp["r"], r_ref)
    assert omp_code.count("!$omp parallel do") == omp_code.count(
        "!$omp end parallel do"
    )
    assert omp_code.count("!$omp parallel do") > 0
    assert "!$omp parallel default" not in omp_code
    assert "!$omp single" not in omp_code
    # The GEMM calls are never inside a parallel construct.
    depth = 0
    for line in omp_code.splitlines():
        stripped = line.strip()
        if stripped.startswith("!$omp parallel do"):
            depth += 1
        elif stripped.startswith("!$omp end parallel do"):
            depth -= 1
        elif "call dgemm(" in stripped:
            assert depth == 0
        continue
    assert depth == 0

    naive_printer = FortranPrinter(
        openmp=False, default_type="real(kind=8)", explicit_bounds=True
    )
    res_naive, _ = run_fortran(
        tmp_path, naive_printer, eval_seq, params, inputs, outputs, tag="naive"
    )
    assert np.allclose(res_naive["r"], r_ref)
