#ifndef NANOSNAP_FFT_H_
#define NANOSNAP_FFT_H_

#include <complex>

namespace nanosnap {

///
/// Apply forward 1D FFT for real-typed array of 1D signal
/// Apply 1D FFT with length `fft_size` for each `cols` of `nframes` input signals.
///
/// @param[in] signal Input signal(real-value). The number of elements are `nframes * nrows`
/// @param[in] fft_size FFT length.
/// @param[in] nframes The number of frames.
/// @param[in] nrows The number of rows
/// @param[out] output Output signal(complex-value). The number of elements are `(fft_size/2 + 1) * nrows`.
/// @param[in] normalize Normalization mode(default is false = normalize by `1/n`, true = normalize by `1/sqrt(n`)
/// @return true upon success.
///
bool rfft(
  const float *signal,
  const size_t fft_size,
  const size_t nframes,
  const size_t nrows,
  std::complex<float> *output,
  const bool normalize = false);


} // namespace nanosnap


#endif // NANOSNAP_FFT_H_
