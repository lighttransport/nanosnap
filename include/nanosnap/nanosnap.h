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
#include <string>
#include <vector>

#include <cmath>

#include "nanosnap/fft.h"

namespace nanosnap {

///
/// 1D Median filter
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
/// Median filter
///
void medfilt(float x, float *y);

///
/// Read WAV file from a file.
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
/// Write WAV file to a file as PCM format.
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
/// Convert a value in Hertz to Mel.
///
template<typename T>
inline T hz2mel(const T hz)
{
    // TODO(LTE): Use faster log10
    return static_cast<T>(2595) * std::log10(static_cast<T>(1)+hz/static_cast<T>(700));
}

///
/// Convert a value in Mel to Hertz.
///
template<typename T>
inline T mel2hz(const T mel)
{
    // TODO(LTE): Use faster pow
    return static_cast<T>(700)*(std::pow(static_cast<T>(10), (mel/static_cast<T>(2595))-static_cast<T>(1)));
}

///
/// Apply a cepstral filter to the matrix of cepstra.
/// @param[in] cepstra The matrix of mel-cepstra with shape [nframes][ncoeffs]
/// @param[in] nframes The number of frames(columns) in `cepstra`
/// @param[in] ncoeffs The number of coefficients(rows) in `cepstra`
/// @param[out] output Pointer where result will be written. It has same shape of `cepstral` and should have enough memory to store `nframes * ncoeffs` values.
/// @param[in] L (optional) the liftering coefficient to use. Default is 22. L
/// <= 0 disables lifter.
/// @return false when invalid input(e.g. cepstra.size() is not dividable by
/// `ncoeff`.
///
bool lifter(const float *cepstra, const size_t nframes, const size_t ncoeffs,
            float *output, const int L = 22);

///
/// Mel Frequency Cepstral Coefficient based on python_speech_features
///
/// @return true upon success. Return false when error.
bool mfcc(float a);

}  // namespace nanosnap

#endif  // NANOSNAP_H_
