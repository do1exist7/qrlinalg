# Generic-precision QR updates

This directory vendors a modified, generic-precision implementation derived
from Martin Köhler's official
[`qrupdate-ng`](https://gitlab.mpi-magdeburg.mpg.de/koehlerm/qrupdate-ng)
repository. The adaptations select real and complex precision through
`wp_def` and provide explicit interfaces to qrlinalg's bundled BLAS/LAPACK.

The public `qrupdate` module provides generic real/complex interfaces for
`qr1up`, `qrinc`, `qrdec`, `qrinr`, `qrder`, `qrshc`, and `gqvec`. Its
`real(wp)` and `complex(wp)` procedures use the same build-selected `wp_def`
module as qrlinalg.

`BLAS.f` and `LAPACK.f` are the precision-independent bundled Netlib
aggregates required by both the QR updates and qrlinalg's factorization and
solve paths.
They intentionally use D/Z symbol names for build-selected `real(wp)` and
`complex(wp)` routines. Do not link a conventional fixed-double BLAS/LAPACK
into `wp=10` or `wp=16` builds.

The qrupdate-ng-derived implementation is GPL-3.0-or-later; see `COPYING`.
The imported Netlib routines retain their upstream notices and licensing
documentation in the aggregate source files. The repository-wide attribution
and modification summary is in [`../../THIRD_PARTY_NOTICES.md`](../../THIRD_PARTY_NOTICES.md).
