# mpLBFGS

**mpLBFGS** is a header-only C++ library that implements the Polak-Ribière
conjugate gradient descent and the limited-memory BFGS algorithm (L-BFGS) for unconstrained optimization. The library uses multiple precision scalars instead of the default double machine precision. This allows to obtain much more accurate approximations of the minimizers, effectively to arbitrary precision.

# Build instructions and dependencies

You can either take the provided Makefile as an example, or use cmake like so:

`mkdir build && cd build && cmake ..`

The library itself is header only; however, it requires **Eigen** for matrix
computations and **mpfr** for multiple precision floating point scalars. You
will need **gmp** and **mpfr** installed on your system to link the executable.  

## Credits

The code is derived from [libLBFGS](https://github.com/chokkan/liblbfgs) and
[LBFGS++](http://yixuan.cos.name/LBFGSpp/doc/) libraries; the original algorithm
and implementation goes back to Jorge Nocedal.  The multiple precision class
mpreal uses [MPFR C++](http://www.holoborodko.com/pavel/mpfr/).  

**mpLBFGS** is licensed under the GPL license.

Copyright (c) 1990, Jorge Nocedal<br>
Copyright (c) 2007–2010, Naoaki Okazaki<br>
Copyright (c) 2016, Yixuan Qiu<br>
Copyright (c) 2018–2019, Alex Vlasiuk<br>
