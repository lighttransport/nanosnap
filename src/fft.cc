#include "nanosnap/fft.h"
#include "nanosnap/nanosnap.h"

#include "pocketfft.h"

#include <cassert>
#include <cstring>
#include <memory>

#include <iostream>  // dbg

namespace nanosnap {

// 1D real fft.
bool rfft(const float *signal, const size_t fft_size, const size_t nframes,
          const size_t nrows, std::complex<float> *output,
          const bool normalize) {
  if (signal == nullptr) {
    return false;
  }

  if ((nframes < 1) || (nrows < 1)) {
    return false;
  }

  // TODO(LTE): Support fft_size != (nframes-1).
  assert(fft_size == (nframes - 1));

  // to double.
  std::vector<double> dsignal;
  dsignal.resize(nframes * nrows);

  for (size_t i = 0; i < nframes * nrows; i++) {
    dsignal[i] = double(signal[i]);
  }

  std::vector<double> dout;  // (real, img), (real, img), ...
  size_t output_nframes = (fft_size / 2) + 1;
  // 2x is for considering complex-type
  dout.resize(2 * output_nframes * nrows);
  memset(reinterpret_cast<char *>(dout.data()), 0,
         sizeof(double) * 2 * output_nframes * nrows);

  float norm_factor = 1.0f;
  if (normalize) {
    norm_factor = 1.0f / std::sqrt(float(nframes));
  }

  int n = rfft_forward_1d_array(dsignal.data(), int(fft_size), int(nframes),
                                int(nrows), norm_factor, dout.data());

  if (n != int(nrows)) {
    return false;
  }

  // from double.
  for (size_t j = 0; j < nrows; j++) {
    for (size_t i = 0; i < output_nframes; i++) {
      output[j * output_nframes + i] =
          std::complex<float>(float(dout[2 * (j * output_nframes + i) + 0]),
                              float(dout[2 * (j * output_nframes + i) + 1]));
    }
  }

  return true;
}

}  // namespace nanosnap
