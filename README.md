# Random L-function test

This project evaluates a sample Selberg-class L-function using the supplied
coefficients and prints the absolute value along the line `Re(s) = w / 2`.

## Files

- `main.cpp` - sample input and evaluation loop
- `L.h` - L-function implementation
- `math.h` - high-precision numeric types and gamma approximation
- `Random-L-function-test.exe` - supplied prebuilt Windows executable

## Dependencies

Building from source requires a C++20 compiler, Boost headers with
Boost.Multiprecision, and the MPFR, MPC, and GMP libraries used by the
Multiprecision backends in `math.h`. Install/configure these dependencies for
your compiler before building. The `.exe` is prebuilt and does not replace
those build dependencies.
