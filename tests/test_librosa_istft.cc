#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <complex>
#include <iostream>

using namespace doctest;

TEST_CASE("librosa_istft") {
#include "testvector/librosa_istft.inc"

  // Build an array of complex values;
  const size_t input_len = sizeof(g_input) / sizeof(g_input[0]);

  std::vector<std::complex<float>> input;

  for (size_t i = 0; i < input_len / 2; i++) {
    input.push_back(std::complex<float>(g_input[2 * i + 0], g_input[2 * i + 1]));
  }

  std::vector<float> result;

  std::cout << "k_ncolumns = " << k_ncolumns << std::endl;
  std::cout << "k_nrows = " << k_nrows << std::endl;
  bool ret = nanosnap::istft(input.data(), k_ncolumns, k_nrows, k_hop_length, k_win_length, &result);

  CHECK(ret == true);

  const size_t out_len = sizeof(g_reference) / sizeof(g_reference[0]);
  CHECK(result.size() == out_len);

  std::cout << "result = \n";
  for (size_t i = 0; i < result.size(); i++) {
      std::cout << "[" << i << "] " << result[i] << "\n";
  }

  std::cout << "\nreference = \n";
  for (size_t i = 0; i < result.size(); i++) {
    std::cout << "[" << i << "] = " << g_reference[i] << "\n";
  }

  for (size_t i = 0; i < result.size(); i++) {
    CHECK(g_reference[i] == Approx(result[i]));
  }
}
