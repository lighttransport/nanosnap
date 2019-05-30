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
    (*err) +=
        "Failed to open WAV file to write. WAV format(channels/rate) is "
        "invalid or cannot write to disk : " +
        filename + "\n";
    return false;
  }

  // TODO(LTE): Implement

  (void)data;
  (void)samples;

  drwav_close(pwav);

  return true;

#endif  // !NANOSNAP_NO_STDIO
}

}  // namespace nanosnap
