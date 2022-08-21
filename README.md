# mpLBFGS

**mpLBFGS** is a header-only C++ library that implements the Polak-Ribière
conjugate gradient descent and the limited-memory BFGS algorithm (L-BFGS) for
unconstrained optimization. The library uses multiple precision scalars instead
of the default double machine precision. This allows to obtain much more
accurate approximations of the minimizers, effectively to arbitrary precision.

# Build instructions and dependencies

If you have the **Eigen** and **mpfr** libraries installed on your system,
`make` builds the wtpframe target by default. Alternatively, to download the
dependencies first and then build wtpframe, issue

`make all`

To build with cmake:

`mkdir build && cd build && cmake ..`

The library requires **Eigen** for matrix computations and **mpfr** for multiple
precision floating point scalars, respectively. You will need **gmp** and
**mpfr** installed on your system in order to use the multiple precision
scalars. The optimizers are templatized in the scalar type, so modifying that
should be straightforward. 

# Credits

The code is derived from [libLBFGS](https://github.com/chokkan/liblbfgs) and
[LBFGS++](http://yixuan.cos.name/LBFGSpp/doc/) libraries; the original algorithm
and implementation goes back to Jorge Nocedal and his collaborators.  The
multiple precision numbers are provided by the [GNU
MPFR](https://www.mpfr.org/), wrapped in `mpreal` class by [Pavel
Holoborodko](http://www.holoborodko.com/pavel/mpfr/).  

**mpLBFGS** is licensed under the GPL license.
