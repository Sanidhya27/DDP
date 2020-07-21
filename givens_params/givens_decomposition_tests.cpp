#include <iomanip>
#include "givens_decomposition.hpp"
#include "gtest/gtest.h"
#include <cmath>

TEST(Givens, All) {
        int N = 3;
        itpp::cmat V(N, N);
        for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                        V(i, j) = (i == j) ? 1 : 0;
        auto params = givens_decomposition(V);
        for (int i = 0; i < params.t.size(); ++i)
                ASSERT_TRUE(fabs(params.t[i]) < 1e-5);
        for (int i = 0; i < params.p.size(); ++i)
                ASSERT_TRUE(fabs(params.t[i]) < 1e-5);
}
