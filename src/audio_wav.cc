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
#include "nanosnap/nanosnap.h"

#ifdef NANOSNAP_NO_STDIO
#ifndef DR_WAV_NO_STDIO
#define DR_WAV_NO_STDIO
#endif
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wmissing-declarations"
#endif

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

namespace nanosnap {

bool wav_read(const std::string &filename, uint32_t *rate, std::string *dtype,
              uint32_t *channels, uint64_t *samples, std::vector<uint8_t> *data,
              std::string *err) {
#ifdef NANOSNAP_NO_STDIO
  (void)filename;
  (void)rate;
  (void)dtype;
  (void)channels;
  (void)data;

  if (err) {
    (*err) +=
        "IO is disabled in this build. Use `wav_read_from_buffer` instead.\n";
  }

  return false;
#else

  if (!rate) {
    if (err) {
      (*err) += "nullptr is set to `rate` parameter.\n";
    }
    return false;
  }

  if (!channels) {
    if (err) {
      (*err) += "nullptr is set to `channels` parameter.\n";
    }
    return false;
  }

  if (!dtype) {
    if (err) {
      (*err) += "nullptr is set to `dtype` parameter.\n";
    }
    return false;
  }

  if (!data) {
    if (err) {
      (*err) += "nullptr is set to `data` parameter.\n";
    }
    return false;
  }

  drwav *pwav = drwav_open_file(filename.c_str());

  if (!pwav) {
    if (err) {
      (*err) += "Error opening WAV file. File not found or not a WAV file? : " +
                filename + "\n";
    }
    return false;
  }

  if (pwav->bitsPerSample == 8) {
    (*dtype) = "uint8";
  } else if (pwav->bitsPerSample == 16) {
    (*dtype) = "int16";
  } else if (pwav->bitsPerSample == 32) {
    if (pwav->fmt.formatTag == DR_WAVE_FORMAT_IEEE_FLOAT) {
      (*dtype) = "float32";
    } else {
      (*dtype) = "float32";
    }
  } else if ((pwav->bitsPerSample == 64) &&
             (pwav->fmt.formatTag == DR_WAVE_FORMAT_IEEE_FLOAT)) {
    (*dtype) = "float64";
  } else {
    if (err) {
      (*err) = "Unsupported format type. bitsPerSample = " +
               std::to_string(int(pwav->bitsPerSample)) +
               ", format = " + std::to_string(int(pwav->fmt.formatTag)) + "\n";
    }

    drwav_close(pwav);

    return false;
  }

  size_t data_len = size_t(pwav->totalPCMFrameCount * pwav->channels *
                           (pwav->bitsPerSample / 8));
  data->resize(data_len);
  drwav_uint64 frame_bytes_read = drwav_read_pcm_frames(
      pwav, pwav->totalPCMFrameCount, reinterpret_cast<void *>(data->data()));

  if (frame_bytes_read != pwav->totalPCMFrameCount) {
    if (err) {
      (*err) = "The number of frames read(" + std::to_string(frame_bytes_read) +
               ") does not match PCM frames(" +
               std::to_string(pwav->totalPCMFrameCount) + ")\n";
    }

    drwav_close(pwav);

    return false;
  }

  (*channels) = pwav->channels;
  (*rate) = pwav->sampleRate;
  (*samples) = pwav->totalPCMFrameCount;

  drwav_close(pwav);

  return true;

#endif  // NANOSNAP_NO_STDIO
}

bool wav_write(const std::string &filename, const uint32_t rate,
               const std::string &dtype, const uint32_t channels,
               const uint64_t samples, const uint8_t *data, std::string *err) {
#ifdef NANOSNAP_NO_STDIO
  (void)filename;
  (void)rate;
  (void)dtype;
  (void)channels;
  (void)data;

  if (err) {
    (*err) +=
        "IO is disabled in this build. Use `wav_write_to_buffer` instead.\n";
  }

  return false;
#else

  drwav_uint32 bps = 0;
  if (dtype.compare("float32") == 0) {
    bps = 32;
  } else if (dtype.compare("int32") == 0) {
    bps = 32;
  } else if (dtype.compare("int16") == 0) {
    bps = 16;
  } else if (dtype.compare("uint8") == 0) {
    bps = 8;
  }

  assert(bps > 0);

  drwav_data_format data_format;
  data_format.container = drwav_container_riff;
  data_format.format = DR_WAVE_FORMAT_PCM;
  data_format.channels = channels;
  data_format.sampleRate = rate;
  data_format.bitsPerSample = bps;

  drwav *pwav = drwav_open_file_write(filename.c_str(), &data_format);
  if (!pwav) {
    if (err) {
      (*err) +=
          "Failed to open WAV file to write. WAV format(channels/rate) is "
          "invalid or cannot write to disk : " +
          filename + "\n";
    }
    return false;
  }

  drwav_uint64 sz = drwav_write_pcm_frames(pwav, samples, reinterpret_cast<const void *>(data));

  if (sz == 0) {
    if (err)  {
      (*err) += "Failed to write PCM frames.\n";
    }
    return false;
  }

  drwav_close(pwav);

  return true;

#endif  // !NANOSNAP_NO_STDIO
}

}  // namespace nanosnap
