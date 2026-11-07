#pragma once
#include "math.h"
#include <map>

class Lfunction {
private:
	v a;
	f N;
	v rfactors;
	c eps;
	f w;
	f Q;
	int terms;
	int max_deriv;
	map<pair<int, int>, c> phi_cache{};

	vector<v> phi_coeffs;

	//Build coeff table for phi function (inv mellin transform of gamma factors)
	//coeffs are residues of gamma factors at poles
	void init_phi() {
		for (int i = 0; i < rfactors.size(); i++) {
			c res = 2;
			c r = rfactors[i];
			for (c r_ : rfactors) {
				if (r != r_) res *= gamma((r_ - r) / 2);
			}
			phi_coeffs.push_back({ res });
		}
		for (int i = 0; i < rfactors.size(); i++) {
			c res = phi_coeffs[i][0];
			c r = rfactors[i];
			for (f k = 1; k < terms; k++) {
				for (c r_ : rfactors) {
					res *= 2 / (r_ - r - 2 * k);
				}
				phi_coeffs[i].push_back(res);
			}
		}
		return;
	}

	//Compute the phi function at a given point t
	c phi_value(f t, int deriv) {
		c result = 0;
		v term;
		f tsqr = t * t;
		for (c r_ : rfactors) {
			term.push_back(pow(t, r_ - deriv));
		}
		c contribution = 1;
		for (int k = 0; k < terms && (contribution == 0 || abs(contribution) > 1e-200); k++) {
			for (int i = 0; i < rfactors.size(); i++) {
				contribution = phi_coeffs[i][k] * term[i];
				for (int j = 0; j < deriv; j++) {
					contribution *= (rfactors[i] + 2 * k - j);
				}
				result += contribution;
				term[i] *= tsqr;
			}
		}
		return result;
	}

	//Compute the phi function at a given point t, caching to avoid recomputation
	c phi(f t, int deriv) {
		int index = round((double)(t * Q)) - 1;
		pair<int, int> cache_key = make_pair(deriv, index);
		if (!phi_cache.contains(cache_key)) {
			phi_cache[cache_key] = phi_value(t, deriv);
		}
		return phi_cache[cache_key];
	}

	//Compute the G function at a given point s and t, with error tolerance err
	c G(c s, f t, f err) {
		c term = 1 / s;
		c result = 0;
		c contribution = 1;
		for (int k = 1; abs(contribution) > err && k < max_deriv; k++) {
			contribution = term * phi(t, k - 1);
			term *= -t / (s + k);
			result += contribution;
		}
		return result;
	}

public:
	// Constructor of Selberg class L-functions
	Lfunction(v coeff, f condutor, v r, c epsilon, f weight)
		: a(coeff), N(condutor), rfactors(r), eps(epsilon), w(weight) {
		cout << "Initializing..." << endl;
		terms = 4000;
		max_deriv = 1000;
		init_phi();
		Q = N;
		for (c r : rfactors) {
			Q /= pi;
		}
		Q = sqrt(Q);
		cout << "L function initialized with " << a.size() <<
			" coefficients, conductor " << N << ", " << 
			rfactors.size() << " gamma factors, epsilon " << eps << ", weight " << w << endl;
	}

	c evaluate(c s) {
		c result = 0;
		c g = 1;
		c dual_g = 1;
		c term = 1;
		for (c r_ : rfactors) {
			g *= gamma((s + r_) / 2);
			dual_g *= gamma((w - s + r_) / 2);
		}
		f err = abs(g) * 1e-15;
		for (int n = 1; n <= a.size() && abs(term) > 1e-3; n++) {
			c a_n = a[n - 1];
			if (abs(a_n) < 1e-200) continue;
			f t = f(n) / Q;
			term = a_n * (pow(t, -s) - G(s, t, err) / g) +
				eps * a_n * (dual_g * pow(t, s - w) - G(w - s, t, err)) / g;
			result += term;
		}
		result *= pow(Q, -s);
		return result;
	}
};