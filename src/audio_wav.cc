#include "nanosnp.h"

#ifdef NANOSNP_NO_STDIO
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

namespace nanosnp {

bool wav_read(const std::string &filename, int *rate, std::string *dtype,
              int *channels, std::vector<uint8_t> *data, std::string *err) {
#ifdef NANOSNP_NO_STDIO
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
  }

  if (!channels) {
    if (err) {
      (*err) += "nullptr is set to `channels` parameter.\n";
    }
  }

  if (!dtype) {
    if (err) {
      (*err) += "nullptr is set to `dtype` parameter.\n";
    }
  }

  if (!data) {
    if (err) {
      (*err) += "nullptr is set to `data` parameter.\n";
    }
  }

  drwav wav;
  if (!drwav_init_file(&wav, filename.c_str())) {
    if (err) {
      (*err) += "Error opening WAV file. File not found or not a WAV file? : " + filename + "\n";
    }
    return false;
  }

  return true;

#endif  // NANOSNP_NO_STDIO
}

bool wav_write(const std::string &filename, const int rate, const std::string &dtype,
               const int channels, const std::vector<uint8_t> &data, std::string *err) {
#ifdef NANOSNP_NO_STDIO
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
  }

  if (!channels) {
    if (err) {
      (*err) += "nullptr is set to `channels` parameter.\n";
    }
  }

  if (!dtype) {
    if (err) {
      (*err) += "nullptr is set to `dtype` parameter.\n";
    }
  }

  if (!data) {
    if (err) {
      (*err) += "nullptr is set to `data` parameter.\n";
    }
  }

  drwav wav;
  if (!drwav_init_file(&wav, filename.c_str())) {
    if (err) {
      (*err) += "Error opening WAV file. File not found or not a WAV file? : " + filename + "\n";
    }
    return false;
  }

  return true;

#endif  // NANOSNP_NO_STDIO
}

}  // namespace nanosnp
