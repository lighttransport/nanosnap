#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <iostream>

using namespace doctest;

TEST_CASE("random_uniform") {
#include "testvector/random_uniform.inc"

  std::vector<float> result = nanosnap::random_uniform(k_n, k_seed);

  for (size_t i = 0; i < k_n; i++) {
    CHECK(g_reference[i] == Approx(result[i]));
  }
}
