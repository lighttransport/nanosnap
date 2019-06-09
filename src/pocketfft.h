#ifndef NANOSNAP_POCKETFFT_H_
#define NANOSNAP_POCKETFFT_H_

#ifdef __cplusplus
extern "C" {
#endif

///
/// Apply forward 1D FFT for real-typed array of 1D signal
/// Repeat 1D FFT `nrows` times. Assume input signal has a shape of [nrows *
/// nsamples]
///
/// @param[in] input Input signal(real-value)
/// @param[in] fft_len FFT length
/// @param[in] nsamples The number of samples
/// @param[in] nrows The number of rows
/// @param[in] norm_factor Normalization factor(usually 1.0 or
/// `1/sqrt(nsamples)` )
/// @param[out] output Output signal(complex-value). The layout of elemens are:
/// `nrows * ((fft_len / 2) + 1) * 2(real, img)`
/// @return The number of FFTs processed(i.e. nrows = success). Zero or negative
/// value = error.
///
int rfft_forward_1d_array(const double *input, const int fft_len,
                          const int nsamples, const int nrows,
                          const float norm_factor, double *output);

///
/// Apply backword 1D FFT for complex-typed array
/// Repeat 1D IFFT `m` times. Assume input signal has a shape of [nrows *
/// nsamples]
///
/// @param[in] input Input signal(complex-value)
/// @param[in] nsamples The number of samples
/// @param[in] nrows The number of rows
/// @param[in] norm_factor Normalization factor(usually `1/nsamples` )
/// @param[out] output Output signal(complex-value). The shape of elements are
/// same as input.
/// @return The number of FFTs processed(i.e. nrows = success). Zero or negative
/// value = error.
///
int cfft_backward_1d_array(const double *input, const int nsamples,
                           const int nrows, const float norm_factor,
                           double *output);

#ifdef __cplusplus
}
#endif

#endif  // NANOSNAP_POCKETFFT_H_
