#pragma once
#include<cmath>
#include<vector>

#pragma warning(push)
#pragma warning(disable: 4146)

#include<boost/multiprecision/mpc.hpp>
#include<boost/multiprecision/mpfr.hpp>

using namespace std;
using namespace boost::multiprecision;

using f = number<mpfr_float_backend<2000>>;
using c = number<mpc_complex_backend<2000>>;
using v = vector<c>;

//constants
const f pi = f("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679");
const f e = f("2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274");

c gamma(c s) {
    if (s.real() < 10) return gamma(s + 1) / s;
    c result = sqrt(2 * pi / s) * pow(s / e, s);
    v numerator = { 1, -1, 1, -1, 1, -691, 1, -3617, 43867, -174611, 77683, -236364091 };
    v denominator = { 12, 360, 1260, 1680, 1188, 360360, 156, 122400, 244188, 125400, 5796, 1506960 };
    c err_term = 0;
    c s_ = s;
    c s_sqr = s * s;
    for (int i = 0; i < 12; ++i) {
        err_term += numerator[i] / (s_ * denominator[i]);
        s_ *= s_sqr;
    }
	result *= exp(err_term);
    return result;
}