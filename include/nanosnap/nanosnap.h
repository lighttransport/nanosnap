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
/// @brief 1D Median filter
///
/// For k = 2m+1, y(i) is the median of x(i-m:i+m).
///
/// @param[in] n The number of elements in `x`
/// @param[in] x Input
/// @param[in] k Window size. Must be odd. Usually 3.
/// @param[in] y Output
/// @param[in] include_nan Optional. Omit NaN value when false.
/// @param[in] padding Optional. zero padding if true. false = truncate.
///
/// @return false when Window size is not odd number.
///
bool medfilt1(const size_t n, const float *x, const int k, float *y,
              bool include_nan = false, bool padding = true);

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
                                  const size_t n, const size_t seed);

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
                                 const size_t n, const size_t seed);

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
                                  const size_t seed);

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
void random_shuffle(float *x, const size_t n, const size_t seed);

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
/// @param[in] cep_lifter apply a lifter to final cepstral coefficients. 0 is no lifter. Default is 22.
/// @param[in] append_energy if this is true, the zeroth cepstral coefficient is replaced with the log of the total frame energy.
/// @param[in] winfunc Windowing function. Use `fIdentityWindowingFunction` if
/// you don't need this feature.
/// @param[out] mfcc MFCC features.
///
/// @return true upon success. Return false when error.
bool mfcc(const float *signal, const size_t sig_len, const float samplerate,
          const float winlen, const float winstep, const size_t ncep,
          const size_t nfilt, const size_t nfft, const float low_freq,
          const float high_freq, const float preemph, const size_t cep_lifter,
          const bool append_energy,
          const std::function<float(int)> winfunc,
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
/// @return The number of frames(positive value). Negative value means there was an error.
///
ssize_t fbank(const float *signal, const size_t nframes, const float samplerate,
           const float winlen, const float winstep,
           const std::function<float(int)> winfunc,
           const int nfilt, const int nfft, const float lowfreq,
           const float highfreq, const float preemph,
           std::vector<float> *features, std::vector<float> *energies);

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
/// @return The number of frames(positive value). Negative value means there was an error.
///
ssize_t logfbank(const float *signal, const size_t nframes, const float samplerate,
              const float winlen, const float winstep,
              const std::function<float(int)> winfunc,
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
/// @return The number of frames(positive). Negative value when there was an error.
///
ssize_t ssc(const float *signal, const size_t sig_len, const int samplerate,
         const float winlen, const float winstep,
         const int nfilt,
         const int nfft, const int lowfreq, const int highfreq,
         const float preemph,
        std::function<float(int)> winfunc,
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
/// @param[in] samplerate The sample rate of the signal we are working with, in
/// Hz.
/// @param[in] winlen The length of the analysis window in seconds.
///
inline int calculate_nfft(const float samplerate, const float winlen) {
  const float window_length_samples = winlen * samplerate;

  int nfft = 1;

  while (float(nfft) < window_length_samples) {
    nfft *= 2;
  }

  return nfft;
}

}  // namespace nanosnap

#endif  // NANOSNAP_H_
