//
// NanoSNAP, nanoscale signal, noise and audio processing library in C++11.
//

/*
The MIT License (MIT)

Copyright (c) 2019 Light Transport Entertainment, Inc.

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
#ifndef NANOSNAP_H_
#define NANOSNAP_H_

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include <cmath>

#include "nanosnap/fft.h"

namespace nanosnap {

///
/// @brief NanoSNAP.
///

///
/// @brief pad 1D array with border mode `reflect`
///
/// numpy.pad()
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.pad.html
///
/// @param[in] input Input 1D array
/// @param[in] n The number of elements in `input`.
/// @param[in] pad_width_before Padding width to left side.
/// @param[in] pad_width_after Padding width to right side.
/// @param[out] output Padded output array
/// @return true upon success.
///

bool pad_reflect(const float *input, const size_t n,
                 const size_t pad_width_before, const size_t pad_width_after,
                 std::vector<float> *output);

///
/// @brief pad 1D array with border mode `constant`
///
/// numpy.pad()
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.pad.html
///
/// Supports sigle constant value for each border.
///
/// @param[in] input Input 1D array
/// @param[in] n The number of elements in `input`.
/// @param[in] pad_width_before Padding width to left side.
/// @param[in] pad_width_after Padding width to right side.
/// @param[out] output Padded output array
/// @param[in] pad_constant_value Padding constant value. Default 0.0
/// @return true upon success.
///

bool pad_constant(const float *input, const size_t n,
                  const size_t pad_width_before, const size_t pad_width_after,
                  std::vector<float> *output,
                  const float pad_constant_value = 0.0f);

///
/// @brief Create an array with the given shape and stride
///
/// numpy.lib.stride_tricks.as_strided()
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.lib.stride_tricks.as_strided.html
///
/// In contrast to numpy's `as_strided`. NanoSNAP implementation creates new
/// array, not a view. Input must be 1D and output shape must be up to 2D.
///
/// @param[in] x Input 1D array
/// @param[in] n The number of elements in `x`
/// @param[in] shape Shape information(2D).
/// @param[in] strides Strides(in bytes) information(2D).
/// @param[out] output Output array
///
/// @return true upon success.
///
bool reshape_with_strides(const float *x, const size_t n, const size_t shape[2],
                          const size_t strides[2], std::vector<float> *output);

///
/// @brief 1D Median filter
///
/// For k = 2m+1, y(i) is the median of x(i-m:i+m).
///
/// @param[in] x Input
/// @param[in] n The number of elements in `x`
/// @param[in] k Window size. Must be odd. Usually 3.
/// @param[out] y Output
/// @param[in] include_nan Optional. Omit NaN value when false.
/// @param[in] padding Optional. zero padding if true. false = truncate.
///
/// @return false when Window size is not odd number.
///
bool medfilt1(const float *x, const size_t n, const int k,
              std::vector<float> *y, bool include_nan = false,
              bool padding = true);

///
/// @brief Filter data along one-dimension with an IIR or FIR filter.
///
/// https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.signal.lfilter.html
///
/// Supports 1D or 2D array for `x`.
/// TODO(LTE): Support axis, zi
///
/// @param[in] b The numerator coefficient vector in a 1-D sequence.
/// @param[in] nb The length of `a`.
/// @param[in] a The denominator coefficient vector in a 1-D sequence. If a[0]
/// is not 1, then both a and b are normalized by a[0].
/// @param[in] na The length of `b`.
/// @param[in] x Input array
/// @param[in] nx The number of columns of `x`.
/// @param[in] mx The number of rows of `x`. 1 for 1D array.
///
/// @return true upon success.
///
bool lfilter(const float *b, const size_t nb, const float *a, const size_t na,
             const float *x, const size_t nx, const size_t mx);

///
/// @brief Returns the discrete, linear convolution of two 1D sequences.
///
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.convolve.html
///
/// If `v` is longer than `a`, the arrays are swapped before computation.
///
/// Note that `v` is accessed in reverse manner as described in numpy's
/// document.
///
/// @param[in] a First 1D input([n])
/// @param[in] v Second 1D input([m])
/// @param[in] mode full(0), same(1) or valid(2) (See numpy document for
/// details)
/// @param[out] output Discrete, linear convolution of `a` and `v`.
///
/// @return Discrete, linear convolution of `a` and `v`
///
bool convolve(const float *a, const size_t n, const float *v, const size_t m,
              std::vector<float> *output, const int mode);

///
/// @brief Get Hann windowing function
///
/// scipy.signal.windows.hann
/// https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.hann.html#scipy.signal.windows.hann
///
/// @param[in] m Number of sample points in output window
/// @param[in] symmetric Optional. true to create symmetric window. false to
/// create periodic window. Default true.
/// @return window function.
///
std::vector<float> window_hann(const size_t m, const bool symmetric = true);

///
/// @brief Get windowing function
///
/// scipy.signal.get_window
/// https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html
///
/// @param[in] window_type Name of Window function type. Currently only 'hann'
/// is implemented.
/// @param[in] nx The number of sampled in the window.
/// @param[out] output Samples of windowing function.
/// @param[in] periodic Optional. false to create symmertic window. Corresponds
/// to `fftbins` in `scipy.signal.get_window` parameter. Default true.
/// @return true upon succes. False when error(e.g. window_type is invalid or
/// unimplemented)
///
bool get_window(const std::string &window_type, const size_t nx,
                std::vector<float> *output, const bool periodic = true);

///
/// Generate sequence of uniformly random values of `n` elements.
/// Equivalent to `numpy.random.uniform`.
///
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.uniform.htm
///
/// This function actually generates deterministic values based on `seed` value.
///
/// Underling algorithm to generate random numbers is
/// MersenneTwister(std::mt19937).
///
/// @param[in] lowval Lower boundary of the output interval.
/// @param[in] highval Upper boundary of the output interval.
/// @param[in] n The number of RNG elements to generate.
/// @param[in] seed Seed value for RNG generator.
/// @return Generated random number values.
///
std::vector<float> random_uniform(const float lowval, const float highval,
                                  const size_t n, const uint32_t seed);

///
/// Draw random samples from a normal (Gaussian) distribution.
/// Equivalent to `numpy.random.normal`.
///
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.normal.htm
///
/// This function actually generates deterministic values based on `seed` value.
///
/// Underling algorithm to generate random numbers is
/// MersenneTwister(std::mt19937).
///
/// @param[in] mean Mean of the distribution. Default 0.0
/// @param[in] stddev Standard deviation of the distribution. Default 1.0
/// @param[in] n The number of RNG elements to generate.
/// @param[in] seed Seed value for RNG generator.
/// @return Drawn random samples.
///
std::vector<float> random_normal(const float mean, const float stddev,
                                 const size_t n, const uint32_t seed);

///
/// Shuffle array elements in random.
/// Equivalent to `numpy.random.shuffle`.
///
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.shuffle.htm
///
/// In contrast to `numpy.random.shuffle`, this function does not shuffle an
/// array in-place.
///
/// This function actually does deterministic shuffle(use deterministic RNGs)
/// based on `seed` value.
///
/// Underling algorithm to generate random numbers is
/// MersenneTwister(std::mt19937).
///
/// @param[in] x Input array.
/// @param[in] seed Seed value for RNG generator.
/// @return Shuffled array.
///
std::vector<float> random_shuffle(const float *x, const size_t n,
                                  const uint32_t seed);

///
/// Shuffle array elements randomly in-place.
/// Equivalent to `numpy.random.shuffle`.
///
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.shuffle.htm
///
/// This function actually does deterministic shuffle(use deterministic RNGs)
/// based on `seed` value.
///
/// Underling algorithm to generate random numbers is
/// MersenneTwister(std::mt19937).
///
/// @param[in,out] x Input/output array.
/// @param[in] seed Seed value for RNG generator.
///
void random_shuffle(float *x, const size_t n, const uint32_t seed);

///
/// @brief Read WAV file from a file.
///
/// API refers to `scipy.io.wavfile.read`
/// https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.wavfile.read.html
///
/// Data type is one of `float32`, `int32`, `int16`, `uint8` or `float64`.
/// After readig data, user must cast retured `data` depending on its dtype.
///
///
/// @param[in] filename WAV filename
/// @param[out] rate Sampling rate
/// @param[out] dtype Data type in string
/// @param[out] channels Audio channels(1 = mono, 2 = stereo)
/// @param[out] samples The number of audio samples.
/// @param[out] data WAV data(opaque binary data)
/// @param[out] err Error or warning message(optional)
///
/// @return true upon success. Return false when error or NanoSNAP is compiled
/// with `NANOSNAP_NO_STDIO`.
///
bool wav_read(const std::string &filename, uint32_t *rate, std::string *dtype,
              uint32_t *channels, uint64_t *samples, std::vector<uint8_t> *data,
              std::string *err);

///
/// @brief Write WAV file to a file as PCM format.
///
/// API refers to `scipy.io.wavfile.write`
/// https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.wavfile.read.html
///
/// Data type is one of `float32`, `int32`, `int16` or `uint8`.
/// After readig data, user must cast retured `data` depending on its dtype.
///
/// @param[in] filename WAV filename
/// @param[in] rate Sampling rate
/// @param[in] dtype Data type in string
/// @param[in] channels Audio channels(1 = mono, 2 = stereo)
/// @param[in] samples The number of audio samples.
/// @param[in] data WAV data(opaque binary data)
/// @param[out] err Error or warning message(optional)
///
/// @return true upon success. Return false when error or NanoSNAP is compiled
/// with `NANOSNAP_NO_STDIO`.
///
bool wav_write(const std::string &filename, const uint32_t rate,
               const std::string &dtype, const uint32_t channels,
               const uint64_t samples, const uint8_t *data, std::string *err);

// ----------------------------------------------
// From python_speech_features
// ----------------------------------------------

///
/// @brief Convert a value in Hertz to Mel.
///
template <typename T>
inline T hz2mel(const T hz) {
  // TODO(LTE): Use faster log10
  return static_cast<T>(2595) *
         std::log10(static_cast<T>(1) + hz / static_cast<T>(700));
}

///
/// @brief Convert a value in Mel to Hertz.
///
template <typename T>
inline T mel2hz(const T mel) {
  // TODO(LTE): Use faster pow
  return static_cast<T>(700) *
         (std::pow(static_cast<T>(10),
                   (mel / static_cast<T>(2595)) - static_cast<T>(1)));
}

///
/// @brief Default windowing function.
///
inline float fIdentityWindowingFunction(size_t n) {
  (void)n;
  return 1.0f;
}

///
/// @brief Apply a cepstral filter to the matrix of cepstra.
///
/// @param[in] cepstra The matrix of mel-cepstra with shape [nframes][ncoeffs]
/// @param[in] nframes The number of frames(columns) in `cepstra`
/// @param[in] ncoeffs The number of coefficients(rows) in `cepstra`
/// values.
/// @param[out] output Liftered output
/// @param[in] L (optional) the liftering coefficient to use. Default is 22. L
/// <= 0 disables lifter.
/// @return false when invalid input(e.g. cepstra.size() is not dividable by
/// `ncoeff`.
///
bool lifter(const float *cepstra, const size_t nframes, const size_t ncoeffs,
            std::vector<float> *output, const int L = 22);

///
/// Mel Frequency Cepstral Coefficient based on python_speech_features
///
/// @param[in] signal Input signal(1D).
/// @param[in] sig_len The length(samples) of signal.
/// @param[in] samplerate Sampling rate [Hz] Default 16000(16k Hz)
/// @param[in] winlen the length of the analysis window in seconds. default
/// 0.025 (25 [ms])
/// @param[in] winstep the step between successive windows in seconds. default
/// 0.01 (10 [ms])
/// @param[in] ncep The number of cepstrum to return. default 13.
/// @param[in] nfilt The number of filter banks. default 26.
/// @param[in] nfft The FFT size. Default 512.
/// @param[in] low_freq lowest band edge of mel filters. [Hz]. default 0.0.
/// @param[in] high_freq highest band edge of mel filters. [Hz]. Default
/// `samplerate/2`.
/// @param[in] preemph apply preemphasis filter with preemph as coefficient. 0
/// is no filter. default 0.97.
/// @param[in] cep_lifter apply a lifter to final cepstral coefficients. 0 is no
/// lifter. Default is 22.
/// @param[in] append_energy if this is true, the zeroth cepstral coefficient is
/// replaced with the log of the total frame energy.
/// @param[in] winfunc Windowing function. Use `fIdentityWindowingFunction` if
/// you don't need this feature.
/// @param[out] mfcc MFCC features.
///
/// @return true upon success. Return false when error.
bool mfcc(const float *signal, const size_t sig_len, const float samplerate,
          const float winlen, const float winstep, const size_t ncep,
          const size_t nfilt, const size_t nfft, const float low_freq,
          const float high_freq, const float preemph, const size_t cep_lifter,
          const bool append_energy, const std::function<float(int)> winfunc,
          std::vector<float> *mfcc);

///
/// \brief  Compute Mel-filterbank energy features from an audio signal.
///
/// Default values are from `python_speech_features`
///
/// @param[in] signal Input signal(1D).
/// @param[in] nframes The number of frames(samples) in input `signal`.
/// @param[in] samplerate Sampling rate [Hz] Default 16000(16k Hz)
/// @param[in] winlen the length of the analysis window in seconds. default
/// 0.025 (25 [ms])
/// @param[in] winstep the step between successive windows in seconds. default
/// 0.01 (10 [ms])
/// @param[in] winfunc Windowing function. Use `fIdentityWindowingFunction` if
/// you don't need this feature.
/// @param[in] nfilt The number of filter banks. default 26.
/// @param[in] nfft The FFT size. Default 512.
/// @param[in] lowfreq lowest band edge of mel filters. [Hz]. default 0.0.
/// @param[in] highfreq highest band edge of mel filters. [Hz]. Default
/// `samplerate/2`.
/// @param[in] preemph apply preemphasis filter with preemph as coefficient. 0
/// is no filter. default 0.97.
/// @param[out] features Output features.
/// @param[out] energies Optional. Output energies(can be nullptr).
///
/// @return The number of frames(positive value). Negative value means there was
/// an error.
///
int64_t fbank(const float *signal, const size_t nframes, const float samplerate,
              const float winlen, const float winstep,
              const std::function<float(int)> winfunc, const int nfilt,
              const int nfft, const float lowfreq, const float highfreq,
              const float preemph, std::vector<float> *features,
              std::vector<float> *energies);

///
/// \brief  Compute log Mel-filterbank energy features from an audio signal.
///
/// Default values are from `python_speech_features`
///
/// @param[in] signal Input signal(1D).
/// @param[in] nframes The number of frames(samples) in input `signal`.
/// @param[in] samplerate Sampling rate [Hz] Default 16000(16k Hz)
/// @param[in] winlen the length of the analysis window in seconds. default
/// 0.025 (25 [ms])
/// @param[in] winstep the step between successive windows in seconds. default
/// 0.01 (10 [ms])
/// @param[in] winfunc Windowing function. Use `fIdentityWindowingFunction` if
/// you don't need this feature.
/// @param[in] nfilt The number of filter banks. default 26.
/// @param[in] nfft The FFT size. Default 512.
/// @param[in] lowfreq lowest band edge of mel filters. [Hz]. default 0.0.
/// @param[in] highfreq highest band edge of mel filters. [Hz]. Default
/// `samplerate/2`.
/// @param[in] preemph apply preemphasis filter with preemph as coefficient. 0
/// is no filter. default 0.97.
/// @param[out] features Output features.
/// @param[out] energies Optional. Output energies(can be nullptr).
///
/// @return The number of frames(positive value). Negative value means there was
/// an error.
///
int64_t logfbank(const float *signal, const size_t nframes,
                 const float samplerate, const float winlen,
                 const float winstep, const std::function<float(int)> winfunc,
                 const int nfilt, const int nfft, const float lowfreq,
                 const float highfreq, const float preemph,
                 std::vector<float> *features, std::vector<float> *energies);

///
/// @brief Compute Spectral Subband Centroid features from an audio signal.
///
/// @param[in] signal Input signal(1D).
/// @param[in] sig_len The length(samples) of input `signal`.
/// @param[in] samplerate Sampling rate [Hz] Default 16000(16k Hz)
/// @param[in] winlen the length of the analysis window in seconds. default
/// 0.025 (25 [ms])
/// @param[in] winstep the step between successive windows in seconds. default
/// 0.01 (10 [ms])
/// @param[in] nfilt The number of filter banks. default 26.
/// @param[in] nfft The FFT size. Default 512.
/// @param[in] lowfreq lowest band edge of mel filters. [Hz]. default 0.0.
/// @param[in] highfreq highest band edge of mel filters. [Hz]. Default
/// `samplerate/2`.
/// @param[in] preemph apply preemphasis filter with preemph as coefficient. 0
/// is no filter. default 0.97.
/// @param[in] winfunc Windowing function. Use `fIdentityWindowingFunction` if
/// you don't need this feature.
/// @param[out] features Output SSC. The length = nframes(return value) * nfilt.
///
/// @return The number of frames(positive). Negative value when there was an
/// error.
///
int64_t ssc(const float *signal, const size_t sig_len, const int samplerate,
            const float winlen, const float winstep, const int nfilt,
            const int nfft, const int lowfreq, const int highfreq,
            const float preemph, std::function<float(int)> winfunc,
            std::vector<float> *features);

///
/// @brief Calculates the FFT size as a power of two greater than or equal to
/// the number of samples in a single window length.
///
///    Having an FFT less than the window length loses precision by dropping
///    many of the samples; a longer FFT than the window allows zero-padding
///    of the FFT buffer which is neutral in terms of frequency domain
///    conversion.
///
/// @param[in] sample_rate The sample rate of the signal we are working with, in
/// Hz.
/// @param[in] winlen The length of the analysis window in seconds.
///
inline int calculate_nfft(const float sample_rate, const float winlen) {
  const float window_length_samples = winlen * sample_rate;

  int nfft = 1;

  while (float(nfft) < window_length_samples) {
    nfft *= 2;
  }

  return nfft;
}

///
/// @brief Create a Filterbank matrix to combine FFT bins into Mel-frequency
/// bins
///
/// Implements `librosa.filters.mel`
///
/// @param[in] sample_rate Sampling rate of input signal.
/// @param[in] n_fft Number of FFT components.
/// @param[out] M Mel transform matrix(shape = [n_mels][1 + nfft/2]).
/// @param[in] n_mel Number of Mel bands to generate(default = 128).
/// @param[in] fmin lowest frequency([Hz]). default = 0.0
/// @param[in] fmax highest frequency(([Hz]). -1.0 = use `sample_rate / 2`.
/// @param[in] htk use HTK formula instead of Slaney(default false)`.
/// @param[in] norm if true, divide the triangular mel weights by the width of
/// the mel band (area normalization). Otherwise, leave all the triangles aiming
/// for a peak value of 1.0(default true)`.
///
/// @return true upon success. Return false when error.
///
bool mel_filter(const float sample_rate, const int n_fft, std::vector<float> *M,
                const int n_mel = 128, const float fmin = 0.0f,
                const float fmax = -1.0f, const bool htk = false,
                const bool norm = true);

///
/// @brief Load 1D or 2D array data in ASCII format.
/// TODO(LTE): Support delimiter character.
///
/// @param[in] filename File name
/// @param[out] values Array of loaded data.
/// @param[out] n The number of columns
/// @param[out] m The number of rows
/// @param[out] err Error message(if exists)
///
/// @return true upon success. false when error(and `err` will be filled)
///
bool loadtxt(const std::string &filename, std::vector<float> *values, int *n,
             int *m, std::string *err);

///
/// @brief Save 1D or 2D array data in ASCII format.
/// TODO(LTE): Support delimiter character.
///
/// @param[in] filename File name
/// @param[in] values Array of data.
/// @param[in] n The number of columns
/// @param[in] m The number of rows
/// @param[out] err Error message(if exists)
///
/// @return true upon success. false when error(and `err` will be filled)
///
bool savetxt(const std::string &filename, const float *values, const int n,
             const int m, std::string *err);

///
/// @brief resize image.
///
/// Resize 32bit float image with bilinear interpolation.
/// Pixel data in memory are LLLL...(gray scale), RGBRGBRGB...(RGB) or
/// RGBARGBARGBA...(RGBA) (This means HWC format in ML field).
///
/// Behavior of interpolation matches with cv2.resize.
/// (TensorFlow resize_image v2 or PyTorch 0.4+)
///
/// Assume image is in linear color space.
/// Up to 2GB of image data.
///
/// @param[in] src Input source image.
/// @param[in] src_width Width of source image in pixels.
/// @param[in] src_width_stride Width stride of source image in pixels.
/// `src_width_stride` must be greater or equal to `width`. Map to `src_width`
/// when set to `0`
/// @param[in] src_height Height of source image in pixels.
/// @param[in] channels The number of pixel channels. 1(grayscale), 3(RGB) or
/// @param[out] dst Output resized image.
/// @param[out] dst_width Width of source image in pixels.
/// @param[out] dst_width_stride Width stride of source image in pixels.
/// `src_width_stride` must be greater or equal to `width`. Map to `src_width`
/// when set to `0`
/// @param[out] dst_height Height of source image in pixels.
/// 4(RGBA) are supported.
/// @return true upon success.
///
bool resize_bilinear(const float *src, const int32_t src_width,
                     const int32_t src_width_stride, const int32_t src_height,
                     const int32_t channels, const int32_t dst_width,
                     const int32_t dst_width_stride, const int32_t dst_height,
                     std::vector<float> *dst);

///
/// @brief Load LDR image from a file.
///
/// Load LDR(e.g. JPG/PNG) image and convert to 32bit float image.
/// [0, 255] are mapped to [0.0, 1.0]
///
/// Image data are stored in memory with HWC format(e.g. RGBRGBRGB...).
///
/// sRGB to Linear conversion is applied.
///
/// @param[in] filename Filename.
/// @param[out] image Image data.
/// @param[out] width Width of loaded image.
/// @param[out] height Height of loaded image.
/// @param[out] channels The number of pixel channels in loaded image(e.g. 3 for
/// RGB)
/// @param[in] srgb_to_linear Optional. Apply sRGB to Linear conversion. Default
/// true. false to load image as is.
///
bool imread(const std::string &filename, std::vector<float> *image,
            int32_t *width, int32_t *height, int32_t *channels,
            const bool srgb_to_linear = true);

///
/// @brief Save image to a file with LDR format.
///
/// Assume input image is in linear space.
/// Pixel value [0.0, 1.0] are mapped to [0, 255].
/// Linear to sRGB conversion is applied.
///
/// Image data are stored in memory with HWC format(e.g. RGBRGBRGB...).
///
/// TODO(LTE): Support width stride.
///
/// @param[in] filename Filename.
/// @param[in] image Image data.
/// @param[in] width Width of loaded image.
/// @param[in] height Height of loaded image.
/// @param[in] channels The number of pixel channels in loaded image(e.g. 3 for
/// RGB)
/// @param[in] linear_to_srgb Optional. Apply Linear to sRGB conversion. Default
/// true. false to load image in linear space.
///
bool imsave(const std::string &filename, std::vector<float> &image,
            const int32_t width, const int32_t height, const int32_t channels,
            const bool linear_to_srgb = true);

}  // namespace nanosnap

#endif  // NANOSNAP_H_
