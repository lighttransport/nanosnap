#include "nanosnap/fft.h"
#include "nanosnap/nanosnap.h"

#include "pocketfft.h"

#include <cassert>
#include <cstring>
#include <memory>
#include <deque>

#include <iostream>  // dbg

namespace nanosnap {

#if 0
namespace {

// Implementation of np.pad() with 'reflect' pad mode.
bool pad_reflect(const float *input, const size_t n, const size_t pad_len_left, const size_t pad_len_right, std::vector<float> *output) {

  //
  //  in [1, 2, 3, 4, 5, 6]
  //  pad = 3, 'reflect' gives,
  //
  //  out [4, 3, 2, 1, 2, 3, 4, 5, 6, 5, 4, 3]
  //

  if (n < 1) {
    return false;
  }

  size_t output_length = n + pad_len_left + pad_len_right;

  if ((pad_len_left == 0) && (pad_len_right == 0)) {
    // just return input

    output->resize(output_length);

    // broadcast input[0]
    for (size_t i = 0; i < output_length; i++) {
      (*output)[i] = input[i];
    }

    return true;
  }

  if (n == 1) {
    output->resize(output_length);

    // broadcast input[0]
    for (size_t i = 0; i < output_length; i++) {
      (*output)[i] = input[0];
    }

    return true;
  }

  // TODO(LTE): Optimize code.
  std::deque<float> buf;

  // left
  if (pad_len_left > 0) {
    int dir = 1; // 1 or -1
    size_t idx = 1;
    for (size_t i = 1; i <= pad_len_left; i++) {
      if ((i % (n - 1)) == 0) {
        // reflect direction
        dir *= -1;
      }

      buf.push_front(input[idx]);

      if (dir > 0) {
        idx++;
      } else {
        idx--;
      }
    }
  }

  // input
  for (size_t i = 0; i < n; i++) {
    buf.push_back(input[i]);
  }

  // right
  if (pad_len_right > 0) {
    int dir = -1; // 1 or -1
    size_t idx = n - 2;
    for (size_t i = 1; i <= pad_len_right; i++) {
      if ((i % (n - 1)) == 0) {
        // reflect direction
        dir *= -1;
      }

      buf.push_back(input[idx]);

      if (dir > 0) {
        idx++;
      } else {
        idx--;
      }
    }
  }

  output->resize(output_length);
  for (size_t i = 0; i < output_length; i++) {
    (*output)[i] = buf[i];
  }

  return true;
}


} // namespace
#endif

// 1D real fft.
bool rfft(const float *signal, const size_t nframes,
          const size_t nrows, const size_t fft_size, std::vector<std::complex<float>> *output,
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
  output->resize(dout.size());
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
bool stft(const float *signal, const size_t sig_len, const size_t n_fft, const size_t hop_length, const size_t win_length, std::vector<std::complex<float>> *output, const bool center)
{
  // TODO(LTE): support pad mode.
  // Assume pad_mode is 'reflect'

  std::vector<float> input;

  if (center) {
    // Pad the time series
    // y = np.pad(y, int(n_fft // 2), mode=pad_mode)
    bool ret = pad_reflect(signal, sig_len, n_fft / 2, n_fft / 2, &input);
    if (!ret) {
      return false;
    }
  } else {
    input.resize(sig_len);
    memcpy(input.data(), signal, sig_len * sizeof(float));
  }

  (void)hop_length;
  (void)win_length;
  (void)output;

#if 0
    # Window the time series.
    y_frames = util.frame(y, frame_length=n_fft, hop_length=hop_length)

    # Pre-allocate the STFT matrix
    stft_matrix = np.empty((int(1 + n_fft // 2), y_frames.shape[1]),
                           dtype=dtype,
                           order='F')

    # how many columns can we fit within MAX_MEM_BLOCK?
    n_columns = int(util.MAX_MEM_BLOCK / (stft_matrix.shape[0] *
                                          stft_matrix.itemsize))

    for bl_s in range(0, stft_matrix.shape[1], n_columns):
        bl_t = min(bl_s + n_columns, stft_matrix.shape[1])

        stft_matrix[:, bl_s:bl_t] = fft.rfft(fft_window *
                                             y_frames[:, bl_s:bl_t],
                                             axis=0)
    return stft_matrix

    std::vector<std::complex<float>> output;
    bool ret = rfft(y, y_nframes, /* 1D array */1,  /* fft size */y_nframes, &output);

  if (ret) {
    return false;
  }
bool rfft(const float *signal, const size_t nframes,
          const size_t nrows, const size_t fft_size, std::vector<std::complex<float>> *output,
          const bool normalize) {
#endif

  return false;
}

}  // namespace nanosnap
