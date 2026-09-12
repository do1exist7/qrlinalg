# Third-party notices

The BSD-3-Clause license in [`LICENSE`](LICENSE) applies to the original
qrlinalg code and documentation. The following bundled components retain
their own copyright notices and licenses.

## qrupdate-ng

The files named `src/qrupdate/qrupdate*.f90` are modified, generic-precision
versions of code derived from Martin Köhler's
[qrupdate-ng](https://gitlab.mpi-magdeburg.mpg.de/koehlerm/qrupdate-ng)
project. qrupdate-ng is itself a modernization of the original
[QRUPDATE](https://sourceforge.net/projects/qrupdate/) project.

The qrupdate-ng-derived files are licensed under the GNU General Public
License, version 3 or, at the recipient's option, any later version
(`GPL-3.0-or-later`). Their copyright and license headers have been retained.
The complete license text is in [`src/qrupdate/COPYING`](src/qrupdate/COPYING).

The copies in this repository were adapted to use qrlinalg's compile-time
working kind and bundled BLAS/LAPACK interfaces. These are local modifications
and are not endorsed by the qrupdate-ng authors.

Because qrlinalg links the qrupdate-ng-derived code into the same library, a
distributed combined binary is subject to the GPL-3.0-or-later terms even
though the original qrlinalg portions are available under BSD-3-Clause.

## Reference BLAS and LAPACK

`src/qrupdate/BLAS.f` and `src/qrupdate/LAPACK.f` aggregate modified routines
from the Netlib reference BLAS and LAPACK distributions. They were adapted for
the compile-time working kind used by qrlinalg; `BLAS.f` also contains local
DGEMM/ZGEMM optimizations and optional OpenMP worksharing.

The upstream copyright and BSD-style license text is reproduced at the top of
each aggregate file. See the [LAPACK license](https://www.netlib.org/lapack/LICENSE.txt)
and the [Netlib LAPACK project](https://www.netlib.org/lapack/) for the
upstream material. These modified copies are not endorsed by the LAPACK or
BLAS authors.
