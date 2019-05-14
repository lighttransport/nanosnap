#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <complex>
#include <iostream>

using namespace doctest;

TEST_CASE("rfft") {
#include "testvector/rfft.inc"

  size_t out_npoints = (k_fft_len / 2) + 1;
  std::vector<std::complex<float>> result(out_npoints * k_nrows);
  bool ret =
      nanosnap::rfft(g_input, k_fft_len, k_nframes, k_nrows, result.data());

  CHECK(ret == true);

  std::cout << "result = \n";
  for (size_t j = 0; j < k_nrows; j++) {
    for (size_t i = 0; i < out_npoints; i++) {
      std::cout << "[" << j << "][" << i
                << "] = " << result[j * out_npoints + i].real() << ", "
                << result[j * out_npoints + i].imag() << "\n";
    }
  }

  std::cout << "\nreference = \n";
  for (size_t j = 0; j < k_nrows; j++) {
    for (size_t i = 0; i < out_npoints; i++) {
      std::cout << "[" << j << "][" << i
                << "] = " << g_reference[2 * (j * out_npoints + i) + 0] << ", "
                << g_reference[2 * (j * out_npoints + i) + 1] << "\n";
    }
  }

  for (size_t i = 0; i < k_nrows * out_npoints; i++) {
    CHECK(g_reference[2 * i + 0] == Approx(result[i].real()));
    CHECK(g_reference[2 * i + 1] == Approx(result[i].imag()));
  }
}
