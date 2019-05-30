///
/// @brief FFT module.
///
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
/// @brief Apply STFT for 1D input
///
/// This is an implementation of librosa.stft
///
/// https://librosa.github.io/librosa/generated/librosa.core.stft.html
///
/// Assume window function is 'hann' and pad_mode is 'reflect'
///
/// TODO(LTE): Support 'pad' mode.
/// TODO(LTE): Support window function.
///
/// @param[in] signal Input signal(real-value).
/// @param[in] nframes The number of frames of input signal.
/// @param[in] n_fft FFT length.
/// @param[in] hop_length Hop length. Default n_fft / 4.
/// @param[in] win_length Window length. Default n_fft.
/// @param[out] output Output signal(complex-value).
/// @param[in] center Optional. Centerize input signal. Default true.
/// @return true upon success.
///
bool stft(const float *signal, const size_t nframes, const size_t n_fft, const size_t hop_length, const size_t win_length, std::vector<std::complex<float>> *output, const bool center = true);

}  // namespace nanosnap

#endif  // NANOSNAP_FFT_H_
