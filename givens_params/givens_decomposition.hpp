#ifndef GIVENS_DECOMPOSITION_HPP
#define GIVENS_DECOMPOSITION_HPP

#include <complex>
#include <itpp/itcomm.h>

typedef struct {
        itpp::vec t;
        itpp::vec p;
} GIVENSPARAMS;

itpp::cmat given_matrix(double theta, int a, int b, int s);
itpp::cmat givens_reconstruction(itpp::vec phis, itpp::vec thetas,int t,int n);
GIVENSPARAMS givens_decomposition(itpp::cmat V);

#endif
