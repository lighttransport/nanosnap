#include "doctest/doctest.h"

#include "nanosnap/nanosnap.h"

#include <complex>
#include <iostream>

using namespace doctest;

TEST_CASE("librosa_filters_mel") {
#include "testvector/librosa_filters_mel.inc"

  std::vector<float> result;

  std::cout << "k_sr = " << k_sr << std::endl;
  std::cout << "k_nfft = " << k_nfft << std::endl;
  bool ret = nanosnap::mel_filter(k_sr, k_nfft, &result);

  CHECK(ret == true);

  const size_t out_len = sizeof(g_reference) / sizeof(g_reference[0]);
  CHECK(result.size() == out_len);

  std::cout << "librosa_filters_mel:result = \n";
  for (size_t i = 0; i < result.size(); i++) {
      std::cout << "[" << i << "] " << result[i] << "\n";
  }

  std::cout << "\nlibrosa_filters_mel:reference = \n";
  for (size_t i = 0; i < result.size(); i++) {
    std::cout << "[" << i << "] = " << g_reference[i] << "\n";
  }

  for (size_t i = 0; i < result.size(); i++) {
    if (g_reference[i] != Approx(result[i])) {
      std::cout << "librosa_filters_mel: diff at " << i << "\n";
    }
    CHECK(g_reference[i] == Approx(result[i]));
  }
}
