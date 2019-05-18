///
/// @brief FFT module.
///
#ifndef NANOSNAP_FFT_H_
#define NANOSNAP_FFT_H_

#include <complex>


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
/// @param[in] fft_size FFT length.
/// @param[in] nframes The number of frames.
/// @param[in] nrows The number of rows
/// @param[out] output Output signal(complex-value). The number of complex-value
/// elements are
/// `(fft_size/2 + 1) * nrows`, thus the caller must pre-allocate enough buffer for this `output` at least `fft_size/2 + 1` elements.
/// @param[in] normalize Normalization mode(default is false. true = normalize by `1/sqrt(n`)
/// @return true upon success.
///
bool rfft(const float *signal, const size_t fft_size, const size_t nframes,
          const size_t nrows, std::complex<float> *output,
          const bool normalize = false);

///
/// Double precision version.
///
bool rfft(const double *signal, const size_t fft_size, const size_t nframes,
          const size_t nrows, std::complex<double> *output,
          const bool normalize = false);

}  // namespace nanosnap

#endif  // NANOSNAP_FFT_H_
