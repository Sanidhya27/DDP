#ifndef UTILS_HPP
#define UTILS_HPP

namespace itpp {

        void logm(const cmat &A, cmat &eA)
        {
                it_error_if(A.rows() != A.cols() || A.rows() < 1, "expm - need square matrix");
                if (A.cols() == 1) {
                        // Square matrix
                        cvec v(1);
                        v(1) = A(0, 0);
                        eA = diag(log(v));
                }
                cmat V;
                cvec d;
                it_error_if(!eig(A, d, V), "expm - error in eig");
                eA = V * diag(log(d)) * hermitian_transpose(V);
        }

        void logm(const mat &A, mat &eA)
        {
                cmat ceA;
                logm(to_cmat(A), ceA);;
                eA = real(ceA);
        }

        void expm(const cmat &A, cmat &eA)
        {
                it_error_if(A.rows() != A.cols() || A.rows() < 1, "expm - need square matrix");
                if (A.cols() == 1) {
                        // Square matrix
                        cvec v(1);
                        v(1) = A(0, 0);
                        eA = diag(exp(v));
                }
                cmat V;
                cvec d;
                it_error_if(!eig(A, d, V), "expm - error in eig");
                eA = V * diag(exp(d)) * hermitian_transpose(V);
        }

        void expm(const mat &A, mat &eA)
        {
                cmat ceA;
                expm(to_cmat(A), ceA);
                eA = real(ceA);
        }
}


#endif
