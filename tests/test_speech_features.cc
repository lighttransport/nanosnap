#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <iostream>

using namespace doctest;

TEST_CASE("lifter") {

#include "testvector/lifter.inc"

  size_t n = k_nframes * k_ncoeffs;

  std::vector<float> result(n);

  bool ret = nanosnap::lifter(g_input, k_nframes, k_ncoeffs, result.data());
  CHECK(ret == true);

  for (size_t i = 0; i < n; i++) {
    CHECK(g_reference[i] == Approx(result[i]));
  }


}

