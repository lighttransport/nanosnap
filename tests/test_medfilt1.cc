#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <iostream>

using namespace doctest;

TEST_CASE("medfilt1") {

#include "testvector/medfilt1.inc"

  float result[k_input_n];
  bool ret = nanosnap::medfilt1(k_input_n, g_input, k_window_size, result);

  CHECK(ret == true);

  for (size_t i = 0; i < k_input_n; i++) {
    CHECK(g_reference[i] == Approx(result[i]));
  }


}

