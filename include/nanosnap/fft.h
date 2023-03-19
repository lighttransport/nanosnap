///
/// @brief FFT module.
///

/*
The MIT License (MIT)

Copyright (c) 2019 - Present Light Transport Entertainment, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef NANOSNAP_FFT_H_
#define NANOSNAP_FFT_H_

#include <complex>
#include <vector>

namespace nanosnap {

///
/// @brief Apply forward 1D FFT for real-typed array of 1D signal
///
/// Apply 1D FFT with length `fft_size` for each `nrows` of `nframes` input
/// signals.
///
/// @code
/// Input
///
/// <--- nframes  ----->
/// +-------------------+
/// |    sample 0       |
/// +-------------------+
/// |    sample 1       |
/// +-------------------+
///         ...
/// +-------------------+
/// |  sample (nrows-1) |
/// +-------------------+
///
///
/// Output
/// <--- 2 * ((nfft_size / 2) + 1) ----->
/// +-----------------------------------+
/// | real, img, real, img, .. .        | sample 0
/// +-----------------------------------+
/// | real, img, real, img, ...         | sample 1
/// +-----------------------------------+
///                 ...
/// +-----------------------------------+
/// | real, img, real, img, ...         | sample (nrows-1)
/// +-----------------------------------+
///
/// @endcode
///
///
/// @param[in] signal Input signal(real-value). The number of elements are
/// `nframes * nrows`
/// @param[in] nframes The number of frames.
/// @param[in] nrows The number of rows
/// @param[in] fft_size FFT length.
/// @param[out] output Output signal(complex-value). The number of complex-value
/// elements are
/// `(fft_size/2 + 1) * nrows`.
/// @param[in] normalize Normalization mode(default is false. true = normalize
/// by `1/sqrt(n`)
/// @return true upon success.
///
bool rfft(const float *signal, const size_t nframes, const size_t nrows,
          const size_t fft_size, std::vector<std::complex<float>> *output,
          const bool normalize = false);

///
/// Double precision version.
///
bool rfft(const double *signal, const size_t nframes, const size_t nrows,
          const size_t fft_size, std::vector<std::complex<double>> *output,
          const bool normalize = false);

///
/// @brief Inverse FFT
///
/// @param[in] input Input array(complex-value). The number of elements are
/// `nframes * nrows`
/// @param[in] ncolumns The number of columns in `a`.
/// @param[in] nrows The number of rows in `a`.
/// @param[in] n Length of the transformed axis of the output.
/// @param[out] output Output signal(complex-value).
/// @return true upon success.
///
bool ifft(const std::complex<float> *input, const size_t ncolumns, const size_t nrows, const size_t n,
  std::vector<std::complex<float>> *output);

///
/// @brief Apply STFT for 1D signal
///
/// This is an implementation of librosa.stft
///
/// https://librosa.github.io/librosa/generated/librosa.core.stft.html
///
/// Assume window function is 'hann' and pad_mode is 'reflect'
///
/// TODO(LTE): Support other 'pad' mode.
/// TODO(LTE): Support other window function.
///
/// @param[in] signal Input signal(real-value).
/// @param[in] nsamples The number of samples in input signal.
/// @param[in] n_fft FFT length.
/// @param[in] hop_length Hop length. Default n_fft / 4.
/// @param[in] win_length Window length. Default n_fft.
/// @param[out] output Output signal(complex-value). The number of rows are
/// calculated by `output.size() / n_fft`.
/// @param[in] center Optional. Centerize input signal. Default true.
/// @return true upon success.
///
bool stft(const float *signal, const size_t nsamples, const size_t n_fft,
          const size_t hop_length, const size_t win_length,
          std::vector<std::complex<float>> *output, const bool center = true);

///
/// @brief Inverse STFT
///
/// This is an implementation of librosa.istft
///
/// https://librosa.github.io/librosa/generated/librosa.core.istft.html
///
/// Assume window function is 'hann'
///
/// TODO(LTE): Support other window function.
///
/// @param[in] stft Input STFT matrix as 1D array(complex-value).
/// @param[in] ncolumns The number of columns in stft matrix(`1 + n_fft /2`).
/// @param[in] nrows The number of rows in stft matrix(`t`).
/// @param[in] hop_length Hop length. Default n_fft / 4.
/// @param[in] win_length Window length. Default n_fft(=(ncolums-1)*2).
/// @param[out] output Output signal(real-value).
/// @param[in] center Optional. Centerize input signal. Default true.
/// @param[in] length Optional. When > 0, set output length to `length`(zero-padded or clipped)
/// @return true upon success.
///
bool istft(const std::complex<float> *stft, const size_t ncolumns,
           const size_t nrows, const size_t hop_length, const size_t win_length,
           std::vector<float> *output, const bool center = true, const int length = 0);

}  // namespace nanosnap

#endif  // NANOSNAP_FFT_H_
