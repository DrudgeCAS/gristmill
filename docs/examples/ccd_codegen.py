"""Generate block-einsum and BLAS-Fortran code for a CCD residual.

This example deliberately keeps the generated source in memory.  Applications
can embed either returned string in their own Python or Fortran driver.
"""

from drudge import NEG, PartHoleDrudge, Perm
from pyspark import SparkConf, SparkContext
from sympy import IndexedBase, Rational

from gristmill import BlasFortranPrinter, EinsumPrinter, optimize


def form_ccd_sequence():
    """Derive and optimize a CCD doubles residual using its LHS symmetry."""
    conf = SparkConf().setMaster("local[2]").setAppName("ccd-codegen")
    ctx = SparkContext.getOrCreate(conf)
    dr = PartHoleDrudge(ctx)
    p = dr.names
    a, b = p.V_dumms[:2]
    i, j = p.O_dumms[:2]

    t = IndexedBase("t")
    dr.set_dbbar_base(t, 2)
    t2 = dr.einst(
        Rational(1, 4)
        * t[a, b, i, j]
        * p.c_dag[a]
        * p.c_dag[b]
        * p.c_[j]
        * p.c_[i]
    )

    curr = dr.ham
    h_bar = dr.ham
    for order in range(2):
        curr = (curr | t2).simplify() / (order + 1)
        h_bar += curr
    h_bar = h_bar.simplify()

    r = IndexedBase("r")
    dr.set_symm(
        r,
        Perm([1, 0, 2, 3], NEG),
        Perm([0, 1, 3, 2], NEG),
        valence=4,
    )
    proj = p.c_dag[i] * p.c_dag[j] * p.c_[b] * p.c_[a]
    target = dr.define(
        r[a, b, i, j], (proj * h_bar).eval_fermi_vev().simplify()
    )
    sequence = optimize(
        [target],
        substs={p.nv: 10 * p.no},
        interm_fmt="tau{}",
        lhs_symm=True,
    )
    return sequence, dr


def generate_code(sequence, dr):
    """Return two alternative source strings for the same CCD sequence."""
    # The Python driver should provide ``einsum`` (for example,
    # ``from pyscf.lib import einsum``), together with NumPy ``zeros`` and
    # ``dtype``.  Full-space Fock and ERI arrays are sliced into occupied and
    # virtual blocks before their first use and released after their last use.
    einsum_code = EinsumPrinter(
        einsum="einsum",
        blocks=[dr.fock, dr.two_body],
        base_indent=0,
    ).doprint(sequence)

    # Full-space Fock and ERI inputs have nonzero virtual lower bounds.  They
    # are packed only when their storage is not directly usable by DGEMM.
    fortran_printer = BlasFortranPrinter(
        gemm="dgemm",
        default_type="real(kind=8)",
        explicit_bounds=True,
        copy_inputs=[dr.fock, dr.two_body],
    )
    declarations, evaluations = fortran_printer.print_decl_eval(sequence)
    fortran_code = "\n".join((*declarations, *evaluations))
    return einsum_code, fortran_code


if __name__ == "__main__":
    seq, ccd_drudge = form_ccd_sequence()
    python_source, fortran_source = generate_code(seq, ccd_drudge)
    print("# Blockwise einsum\n", python_source)
    print("! BLAS Fortran\n", fortran_source)
