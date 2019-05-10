#ifndef NANOSNP_H_
#define NANOSNP_H_

#include <cstdint>
#include <string>
#include <vector>

namespace nanosnp {

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
/// Data type is one of `float32`, `int32`, `int16` or `uint8`.
/// After readig data, user must cast retured `data` depending on its dtype.
///
///
/// @param[in] filename WAV filename
/// @param[out] rate Sampling rate
/// @param[out] dtype Data type in string
/// @param[out] channels Audio channels(1 = mono, 2 = stereo)
/// @param[out] data WAV data(opaque binary data)
/// @param[out] err Error or warning message(optional)
///
/// @return true upon success. Return false when error or NanoSNP is compiled with `NANOSNP_NO_STDIO`.
///
bool wav_read(const std::string &filename, int *rate, std::string *dtype, int *channels, std::vector<uint8_t> *data, std::string *err);

///
/// Write WAV file to a file.
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
/// @param[in] data WAV data(opaque binary data)
/// @param[out] err Error or warning message(optional)
///
/// @return true upon success. Return false when error or NanoSNP is compiled with `NANOSNP_NO_STDIO`.
///
bool wav_write(const std::string &filename, const int rate, const std::string &dtype, const int channels, std::vector<uint8_t> *data, std::string *err);

} // namespace nanosnp

#endif // NANOSNP_H_
