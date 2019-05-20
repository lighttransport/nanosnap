#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <complex>
#include <iostream>

using namespace doctest;

TEST_CASE("convolve1d_full") {
#include "testvector/convolve_full.inc"

  std::vector<float> result(k_outnum);
  bool ret =
      nanosnap::convolve(g_a, k_n, g_v, k_m, &result, /* mode 0 = full */0);

  CHECK(ret == true);

  std::cout << "convolve(full): len = " << sizeof(g_reference) / sizeof(g_reference[0]) << std::endl;

  for (size_t i = 0; i < k_outnum; i++) {
    CHECK(g_reference[i] == Approx(result[i]));
  }
}

TEST_CASE("convolve1d_same") {
#include "testvector/convolve_same.inc"

  std::vector<float> result(k_outnum);
  bool ret =
      nanosnap::convolve(g_a, k_n, g_v, k_m, &result, /* mode 1 = same */1);

  CHECK(ret == true);

  std::cout << "convolve(same): len = " << sizeof(g_reference) / sizeof(g_reference[0]) << std::endl;

  for (size_t i = 0; i < k_outnum; i++) {
    CHECK(g_reference[i] == Approx(result[i]));
  }
}

TEST_CASE("convolve1d_valid") {
#include "testvector/convolve_valid.inc"

  std::vector<float> result(k_outnum);
  bool ret =
      nanosnap::convolve(g_a, k_n, g_v, k_m, &result, /* mode 2 = valid */2);

  CHECK(ret == true);

  std::cout << "convolve(valid): len = " << sizeof(g_reference) / sizeof(g_reference[0]) << std::endl;

  for (size_t i = 0; i < k_outnum; i++) {
    CHECK(g_reference[i] == Approx(result[i]));
  }
}
