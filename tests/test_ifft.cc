#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <complex>
#include <iostream>

using namespace doctest;

TEST_CASE("ifft") {
#include "testvector/ifft.inc"


  size_t input_len = sizeof(g_input) / sizeof(g_input[0]);
  std::cout << "len = " << input_len << std::endl;

  std::vector<std::complex<float>> input;
  for (size_t i = 0; i < input_len / 2; i++) {
    input.push_back(std::complex<float>(g_input[2 * i + 0], g_input[2 * i + 1]));
  }

  std::vector<std::complex<float>> result;

  bool ret =
      nanosnap::ifft(input.data(), k_ncolumns, k_nrows, k_n, &result);

  CHECK(ret == true);
  CHECK(result.size() == (k_n * k_nrows));

  //std::cout << "result = \n";
  //for (size_t j = 0; j < k_nrows; j++) {
  //  for (size_t i = 0; i < k_n; i++) {
  //    std::cout << "[" << j << "][" << i
  //              << "] = " << result[j * k_n + i].real() << ", "
  //              << result[j * k_n + i].imag() << "\n";
  //  }
  //}

  //std::cout << "\nreference = \n";
  //for (size_t j = 0; j < k_nrows; j++) {
  //  for (size_t i = 0; i < out_npoints; i++) {
  //    std::cout << "[" << j << "][" << i
  //              << "] = " << g_reference[2 * (j * out_npoints + i) + 0] << ", "
  //              << g_reference[2 * (j * out_npoints + i) + 1] << "\n";
  //  }
  //}

  for (size_t i = 0; i < k_n * k_nrows; i++) {
    CHECK(g_reference[2 * i + 0] == Approx(result[i].real()));
    CHECK(g_reference[2 * i + 1] == Approx(result[i].imag()));
  }
}
