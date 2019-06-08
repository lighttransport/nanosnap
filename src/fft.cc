#include "nanosnap/fft.h"
#include "nanosnap/nanosnap.h"

#include "pocketfft.h"

#include <cassert>
#include <cstring>
#include <deque>
#include <memory>

#include <iostream>  // dbg

namespace nanosnap {

namespace {

static bool pad_center(const std::vector<float> &data, const size_t n, std::vector<float> *output)
{
  if (data.size() > n) {
    return false;
  }

  const size_t pad_width = (n - data.size()) / 2;

  bool ret = pad_constant(data.data(), data.size(), pad_width, pad_width, output, /* pad constant value */0.0f);

  return ret;
}

// librosa.utils.frame
// 1D signal input only
// output will have shape [out_rows][frame_length]
static bool frame_sig(const std::vector<float> &y,
                      std::vector<float> *output,
                      size_t *out_rows,
                      const size_t frame_length = 2048,
                      const size_t hop_length = 512) {
  if (y.size() < frame_length) {
    return false;
  }

  if (hop_length < 1) {
    return false;
  }

  if (hop_length > frame_length) {
    return false;
  }

  // Compute the number of frames that will fit. The end may get truncated.
  const size_t n_frames = 1 + ((y.size() - frame_length) / hop_length);

  // Vertical stride is one sample
  // Horizontal stride is `hop_length` samples
  // y_frames = as_strided(y, shape=(frame_length, n_frames),
  //                       strides=(y.itemsize, hop_length * y.itemsize))

  // y_frames[j][i] == y[i * hop_length + j]

  std::cout << "y.size = " << y.size() << std::endl;
  std::cout << "frame_sig.n_frames = " << n_frames << std::endl;
  std::cout << "frame_length = " << frame_length << std::endl;
  std::cout << "hop_length = " << hop_length << std::endl;

  output->resize(n_frames * frame_length);
  for (size_t j = 0; j < frame_length; j++) {
    for (size_t i = 0; i < n_frames; i++) {
      (*output)[j * n_frames + i] = y[i * hop_length + j];
    }
  }

  (*out_rows) = n_frames;

  return true;
}

}  // namespace

// 1D real fft.
bool rfft(const float *signal, const size_t nframes, const size_t nrows,
          const size_t fft_size, std::vector<std::complex<float>> *output,
          const bool normalize) {
  if (signal == nullptr) {
    return false;
  }

  if ((nframes < 1) || (nrows < 1)) {
    return false;
  }

  // TODO(LTE): Support fft_size != (nframes-1).
  assert(fft_size == (nframes - 1));

  // to double.
  std::vector<double> dsignal;
  dsignal.resize(nframes * nrows);

  for (size_t i = 0; i < nframes * nrows; i++) {
    dsignal[i] = double(signal[i]);
  }

  std::vector<double> dout;  // (real, img), (real, img), ...
  size_t output_nframes = (fft_size / 2) + 1;
  // 2x is for considering complex-type
  dout.resize(2 * output_nframes * nrows);
  memset(reinterpret_cast<char *>(dout.data()), 0,
         sizeof(double) * 2 * output_nframes * nrows);

  float norm_factor = 1.0f;
  if (normalize) {
    norm_factor = 1.0f / std::sqrt(float(nframes));
  }

  int n = rfft_forward_1d_array(dsignal.data(), int(fft_size), int(nframes),
                                int(nrows), norm_factor, dout.data());

  if (n != int(nrows)) {
    return false;
  }

  // Cast from double to float.
  output->resize(output_nframes * nrows);
  for (size_t j = 0; j < nrows; j++) {
    for (size_t i = 0; i < output_nframes; i++) {
      (*output)[j * output_nframes + i] =
          std::complex<float>(float(dout[2 * (j * output_nframes + i) + 0]),
                              float(dout[2 * (j * output_nframes + i) + 1]));
    }
  }

  return true;
}

// librosa.stft(1D input only)
bool stft(const float *signal, const size_t sig_len, const size_t n_fft,
          const size_t hop_length, const size_t win_length,
          std::vector<std::complex<float>> *output, const bool center) {
  // TODO(LTE): support pad mode.
  // Assume pad_mode is 'reflect'

  std::vector<float> y;

  if (center) {
    // Pad the time series
    // y = np.pad(y, int(n_fft // 2), mode=pad_mode)
    bool ret = pad_reflect(signal, sig_len, n_fft / 2, n_fft / 2, &y);
    if (!ret) {
      return false;
    }
  } else {
    y.resize(sig_len);
    memcpy(y.data(), signal, sig_len * sizeof(float));
  }

  std::vector<float> _window;
  {
    bool ret = get_window("hann", win_length, &_window, /* periodic */true);
    if (!ret) {
      return false;
    }
  }

  // Pad the window out to n_fft size
  std::vector<float> window;
  {
    bool ret= pad_center(_window, n_fft, &window);
    if (!ret) {
      return false;
    }
  }

  for (size_t i = 0; i < window.size(); i++) {
    std::cout << "window[" << i << "] = " << window[i] << "\n";
  }

  std::vector<float> y_frames;
  size_t n_rows = 0;
  {
    bool ret = frame_sig(y, &y_frames, &n_rows, n_fft, hop_length);
    if (!ret) {
      return false;
    }
  }

  // Flip dimension to match numpy(librosa)'s result.
  // TODO(LTE): Rewrite frame_sig so that we don't flip dimension to avoid confusion.

  size_t frame_length = n_rows;
  size_t nframes = n_fft;

  std::cout << "window.size = [" << window.size() << "\n";
  std::cout << "nframes = " << nframes << "\n";
  std::cout << "frame_length = [[" << frame_length << "\n";

  for (size_t i = 0; i < y_frames.size(); i++) {
    std::cout << "y_frames[" << i << "] = " << y_frames[i] << "\n";
  }

  std::vector<float> y_modulated(y_frames.size());

  // Apply window function.
  // window.shape = [frame_length][1]
  // (window * y).shape = [frame_length][nframes] in C++
  for (size_t j = 0; j < nframes; j++) {
    for (size_t i = 0; i < frame_length; i++) {
      //y_modulated[j * frame_length + i] *= window[j];
      // transpose
      y_modulated[j * frame_length + i] = y_frames[i * nframes + j] * window[j];
    }
  }

  for (size_t i = 0; i < y_modulated.size(); i++) {
    std::cout << "modulated y_frames[" << i << "] = " << y_modulated[i] << "\n";
  }

  // HACK
  bool ret = rfft(y_modulated.data(), /* frame len */nframes, /* num frames */frame_length,  /* fft size */nframes, output);

  for (size_t i = 0; i < output->size(); i++) {
    std::cout << "output[" << i << "] = " << (*output)[i] << std::endl;
  }

  return ret;
}

}  // namespace nanosnap
