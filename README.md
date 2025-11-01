![CI](https://github.com/DrudgeCAS/gristmill/workflows/CI/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/DrudgeCAS/gristmill/badge.svg?branch=master)](https://coveralls.io/github/DrudgeCAS/gristmill?branch=master)
[![Cite this repo](https://img.shields.io/badge/cite-CITATION.cff-blue)](./CITATION.cff)

<h1 align="center">Gristmill</h1>

Gristmill is a package built upon the
[drudge](https://github.com/DrudgeCAS/drudge) computer algebra system for
automatic optimization and code generation of tensor computations. While
designed for quantum chemistry and many-body theory, it is suitable for any
scientific computing problem involving tensors.

The optimizer utilizes advanced algorithms to efficiently parenthesize and
factorize tensor computations, reducing the floating-point operation (FLOP)
count. For example, a matrix chain product

$$
\mathbf{R} = \mathbf{A} \mathbf{B} \mathbf{C}
$$

can be parenthesized as

$$
\mathbf{R} = ( \mathbf{A} \mathbf{B} ) \mathbf{C}
$$

or

$$
\mathbf{R} = \mathbf{A} ( \mathbf{B} \mathbf{C} ),
$$

depending on which has fewer FLOPs given the shapes of the matrices. General
tensor contractions are supported, with minimal overhead relative to
specialized dynamic programming algorithms for matrix chain products. As an
example, when evaluating the ladder term in coupled cluster doubles (CCD)
residuals

$$
r_{abij} = \sum_{c,d=1}^v \sum_{k,l=1}^o v_{klcd} t_{cdij} t_{abkl},
$$

Gristmill can automatically generate a two-step contraction scheme:

$$
\begin{aligned}
    p_{klij} &= \sum_{c,d=1}^v v_{klcd} t_{cdij}\\
    r_{abij} &= \sum_{k,l=1}^o p_{klij} t_{abkl}
\end{aligned}
$$

Gristmill's algorithm can efficiently handle complicated contractions with,
e.g., twenty factors.

When evaluating sums of multiple contractions, Gristmill factors each term to
reduce computational cost. For example, the coupled cluster singles and doubles
(CCSD) correlation energy

$$
E = \frac{1}{4} \sum_{i,j=1}^o \sum_{a,b=1}^{v} u_{ijab} t^{(2)}_{abij} + \frac{1}{2} \sum_{i,j=1}^o \sum_{a,b=1}^v u_{ijab} t^{(1)}_{ai} t^{(1)}_{bj}
$$

can be automatically simplified into

$$
E = \frac{1}{4} \sum_{i,j=1}^o \sum_{a,b=1}^v u_{ijab} (t^{(2)}_{abij} + 2 t^{(1)}_{ai} t^{(1)}_{bj}),
$$

which has lower FLOP count.

Additionally, Gristmill includes optimization heuristics such as common
symmetrization, which ensures intermediates that are equivalent by symmetry are
computed only once using the canonicalization capability in
[drudge](https://github.com/DrudgeCAS/drudge).

The code generator in Gristmill is an orthogonal component to the optimizer.
Both optimized and unoptimized computations can be fed into the code generator
to yield naive Fortran and C code (with optional OpenMP parallelization) as well
as Python code using NumPy and Julia code using
[OMEinsum.jl](https://github.com/under-Peter/OMEinsum.jl).


## Installation

Gristmill can be installed directly from the GitHub repository using
[uv](https://github.com/astral-sh/uv) (recommended)
```bash
uv pip install git+https://github.com/DrudgeCAS/gristmill.git
```
or [pip](https://pypi.org/project/pip/)
```bash
pip install git+https://github.com/DrudgeCAS/gristmill.git
```

> **Note:** Native Windows builds are currently not supported. Please use WSL
> (Windows Subsystem for Linux) to install and run Gristmill on Windows for now.


## Documentation

Please refer to the documentation at
[https://drudgecas.github.io/drudge/](https://drudgecas.github.io/drudge/).
Additional examples can be found in the `./docs/examples` directory.


## Citation

If you use Drudge and Gristmill in your work, please cite their GitHub
repositories and Jinmo Zhao's Ph.D. thesis:

**1. The Drudge GitHub repository**  
```bibtex
@misc{DrudgeCAS,
  author       = {Jinmo Zhao and Guo P. Chen and Gaurav Harsha and Matthew Wholey and Thomas M. Henderson and Gustavo E. Scuseria},
  title        = {Drudge: A symbolic algebra system for tensorial and noncommutative algebras},
  publisher    = {GitHub},
  year         = {2016--2025},
  url          = {https://github.com/DrudgeCAS/drudge},
  note         = {GitHub repository}
}
```

**2. The Gristmill GitHub repository**  
```bibtex
@misc{Gristmill,
  author       = {Jinmo Zhao and Guo P. Chen and Gaurav Harsha and Thomas M. Henderson and Gustavo E. Scuseria},
  title        = {Gristmill: A tensor contraction optimizer and code generator based on Drudge},
  publisher    = {GitHub},
  year         = {2016--2025},
  url          = {https://github.com/DrudgeCAS/gristmill},
  note         = {GitHub repository}
}
```

**3. Jinmo Zhao’s Ph.D. thesis**  
```bibtex
@phdthesis{Zhao2018Drudge,
  author       = {Jinmo Zhao},
  title        = {Symbolic Solution for Computational Quantum Many-Body Theory Development},
  school       = {Rice University},
  year         = {2018},
  month        = {April},
  address      = {Houston, Texas, USA},
  type         = {PhD thesis},
  url          = {https://www.proquest.com/openview/61a9a86c07dbb6e5270bdeb1c84384db/1?pq-origsite=gscholar&cbl=18750&diss=y}
}
```
Link: [Symbolic Solution for Computational Quantum Many-Body Theory Development — Jinmo Zhao (2018)](https://www.proquest.com/openview/61a9a86c07dbb6e5270bdeb1c84384db/1?pq-origsite=gscholar&cbl=18750&diss=y)

---

You may also use the [`CITATION.cff`](./CITATION.cff) file provided in this
repository, which is compatible with citation managers such as Zotero and
Mendeley.


## Acknowledgments

Gristmill was originally developed by Jinmo Zhao during his Ph.D. at Rice
University, under the supervision of Prof. Gustavo E. Scuseria. The project was
supported as part of the Center for the Computational Design of Functional
Layered Materials, an Energy Frontier Research Center funded by the U.S.
Department of Energy, Office of Science, Basic Energy Sciences under Award
DE-SC0012575. The package is currently maintained by Guo P. Chen, Gaurav
Harsha, and members of the Scueria group.

