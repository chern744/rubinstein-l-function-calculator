# Random L-function test

This project evaluates a sample Selberg-class L-function using the supplied
coefficients and prints the absolute value along the line `Re(s) = w / 2`.

# Limitations of this algorithm: 

1. L-functions with poles or non-trivial zeros that has a higher order of vanishing won't return right values. 
2. This algorithm is precision sensitive. L-functions with irrational coefficients will hardly converge.
3. Improve the convergence behavior by changing the precision.
4. Author has only tested availability inside critical line. 

## Files

- `main.cpp` - sample input and evaluation loop
- `L.h` - L-function implementation
- `math.h` - high-precision numeric types and gamma approximation

## Dependencies

Building from source requires a C++20 compiler, Boost headers with
Boost.Multiprecision, and the MPFR, MPC, and GMP libraries used by the
Multiprecision backends in `math.h`. Install/configure these dependencies for
your compiler before building. The `.exe` is prebuilt and does not replace
those build dependencies.
