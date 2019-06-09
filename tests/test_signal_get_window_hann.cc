#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <complex>
#include <iostream>

using namespace doctest;

TEST_CASE("signal_get_window_hann") {
#include "testvector/signal_get_window_hann.inc"

  std::vector<float> window;
  bool ret = nanosnap::get_window("hann", k_win_length, &window, /* periodic */true);
  CHECK(ret == true);

  for (size_t i = 0; i < window.size(); i++) {
    CHECK(g_reference[i] == Approx(window[i]));
  }
}
