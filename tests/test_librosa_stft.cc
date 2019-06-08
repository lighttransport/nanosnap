#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <complex>
#include <iostream>

using namespace doctest;

TEST_CASE("librosa_stft") {
#include "testvector/librosa_stft.inc"

  size_t out_npoints = k_n_fft / 2 + 1;
  std::vector<std::complex<float>> result;
  std::cout << "input len = " << sizeof(g_input) / sizeof(g_input[0]) << "\n";
  std::cout << "nframes = " << k_nframes << "\n";

  bool ret =
      nanosnap::stft(g_input, k_nframes, k_n_fft, k_hop_length, k_win_length, &result);

  CHECK(ret == true);

  //std::cout << "result = \n";
  //for (size_t j = 0; j < k_nrows; j++) {
  //  for (size_t i = 0; i < out_npoints; i++) {
  //    std::cout << "[" << j << "][" << i
  //              << "] = " << result[j * out_npoints + i].real() << ", "
  //              << result[j * out_npoints + i].imag() << "\n";
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

  std::cout << "reference len = " << sizeof(g_reference) / sizeof(g_reference[0]) << std::endl;
  std::cout << "n = " << out_npoints << std::endl;
  std::cout << "result len = " << result.size() << std::endl;

  size_t nrows = result.size() / out_npoints;

  CHECK(out_npoints == k_out_dim0);
  CHECK(nrows == k_out_dim1);

  CHECK((2 * result.size()) == (sizeof(g_reference) / sizeof(g_reference[0])));

  //for (size_t i = 0; i < nrows * out_npoints; i++) {
  //  CHECK(g_reference[2 * i + 0] == Approx(result[i].real()));
  //  CHECK(g_reference[2 * i + 1] == Approx(result[i].imag()));
  //}
}
