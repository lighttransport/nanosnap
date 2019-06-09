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
#include "nanosnap/fft.h"
#include "nanosnap/nanosnap.h"

#include "pocketfft.h"

#include <cassert>
#include <cstring>
#include <deque>
#include <memory>
#include <limits>
#include <algorithm>

//#include <iostream>  // dbg

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

  output->resize(n_frames * frame_length);
  for (size_t j = 0; j < frame_length; j++) {
    for (size_t i = 0; i < n_frames; i++) {
      (*output)[j * n_frames + i] = y[i * hop_length + j];
    }
  }

  (*out_rows) = n_frames;

  return true;
}

// librosa.window_sumsquare
// 'hann' window only. norm = 'None'
bool hann_window_sumsquare(
  const size_t n_frames, const size_t hop_length, const size_t win_length, const size_t n_fft,
  std::vector<float> *output)
{
  const size_t n = n_fft + hop_length * (n_frames - 1);
  std::vector<float> x(n, 0.0f);

  //  Compute the squared window at the desired length
  std::vector<float> win_sq;
  {
    std::vector<float> _win_sq;
    bool ret = get_window("hann", win_length, &_win_sq);
    if (!ret) {
      return false;
    }

    // abs(x)^2
    for (size_t i = 0; i < _win_sq.size(); i++) {
      _win_sq[i] = _win_sq[i] * _win_sq[i];
    }

    // pad_center
    ret = pad_center(_win_sq, n_fft, &win_sq);
    if (!ret) {
      return false;
    }
  }

  // Fill the envelope
  // __window_ss_fill(x, win_sq, n_frames, hop_length)
  {
    for (size_t i = 0; i < n_frames; i++) {
      size_t sample = i * hop_length;
      size_t end_idx = (std::min)(x.size(), sample + n_fft);

      // x[sample:min(n, sample + n_fft)] += win_sq[:max(0, min(n_fft, n - sample))]
      for (size_t k = sample; k < end_idx; k++) {
        // clamp for safety
        size_t win_idx = (std::min)(n_fft - 1, k - sample);
        x[k] += win_sq[win_idx];
      }
    }
  }

  (*output) = x; // copy

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

// 1D ifft
bool ifft(const std::complex<float> *input, const size_t ncolumns, const size_t nrows,
          const size_t n, std::vector<std::complex<float>> *output) {
  if (input == nullptr) {
    return false;
  }

  if ((ncolumns < 1) || (nrows < 1)) {
    return false;
  }

  size_t output_n = n;

  // to pair of doubles.
  std::vector<double> dinput;
  dinput.resize(2 * output_n * nrows);

  // clear with zeros.
  memset(dinput.data(), 0, sizeof(double) * dinput.size());

  // Create an array of pair of doubles.
  // Also resize columns
  for (size_t j = 0; j < nrows; j++) {
    for (size_t i = 0; i < output_n; i++) {
      const size_t src_idx = j * ncolumns + i;
      dinput[2 * (j * output_n + i) + 0] = double(input[src_idx].real());
      dinput[2 * (j * output_n + i) + 1] = double(input[src_idx].imag());
    }
  }

  std::vector<double> dout;  // (real, img), (real, img), ...
  // 2x is for considering complex-type
  dout.resize(2 * output_n * nrows);
  memset(reinterpret_cast<char *>(dout.data()), 0,
         sizeof(double) * 2 * output_n * nrows);

  float norm_factor = 1.0f / float(n);

  // TODO(LTE): Implement norm = "ortho"
  // if (normalize) {
  //   norm_factor = 1.0f / sqrt(n);
  // }

  (void)norm_factor;

  int m = cfft_backward_1d_array(dinput.data(), int(output_n), int(nrows), norm_factor, dout.data());

  if (m != int(nrows)) {
    return false;
  }

  // Cast from double to float.
  output->resize(output_n * nrows);
  for (size_t j = 0; j < nrows; j++) {
    for (size_t i = 0; i < output_n; i++) {
      (*output)[j * output_n + i] =
          std::complex<float>(float(dout[2 * (j * output_n + i) + 0]),
                              float(dout[2 * (j * output_n + i) + 1]));
    }
  }

  return true;
}


// librosa.stft(1D input only)
bool stft(const float *signal, const size_t sig_len, const size_t n_fft,
          const size_t hop_length, const size_t win_length,
          std::vector<std::complex<float>> *output, const bool center) {

  // NOTE(LTE): librosa's stft implementation uses Fortran order of array.
  // Need to consider flipping row and columns in appropriate place.
  //
  // TODO(LTE): Rewrite stft without transposing shape to avoid confusion.

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

  std::vector<float> y_frames;
  size_t n_rows = 0;
  {
    bool ret = frame_sig(y, &y_frames, &n_rows, n_fft, hop_length);
    if (!ret) {
      return false;
    }
  }

  // Flip dimension to match numpy(librosa)'s Fortran order.

  size_t frame_length = n_rows;
  size_t nframes = n_fft;

  //std::cout << "window.size = [" << window.size() << "\n";
  //std::cout << "nframes = " << nframes << "\n";
  //std::cout << "frame_length = [[" << frame_length << "\n";
  //for (size_t i = 0; i < y_frames.size(); i++) {
  //  std::cout << "y_frames[" << i << "] = " << y_frames[i] << "\n";
  //}

  std::vector<float> y_modulated(y_frames.size());

  // Apply window function.
  for (size_t j = 0; j < nframes; j++) {
    for (size_t i = 0; i < frame_length; i++) {
      //y_modulated[j * frame_length + i] = y_frames[j * frame_length + i] * window[j];

      // transpose y_modulated
      y_modulated[i * nframes + j] = y_frames[j * frame_length + i] * window[j];
    }
  }

  //for (size_t i = 0; i < y_modulated.size(); i++) {
  //  std::cout << "modulated y_frames[" << i << "] = " << y_modulated[i] << "\n";
  //}

  //std::cout << "nframes = " << nframes << std::endl;
  //std::cout << "frame_len = " << frame_length << std::endl;

  // transpose frame_len and num_frames
  bool ret = rfft(y_modulated.data(), /* frame len */nframes, /* num frames */frame_length,  /* fft size */nframes, output);

  //for (size_t i = 0; i < output->size(); i++) {
  //  std::cout << "output[" << i << "] = " << (*output)[i] << std::endl;
  //}

  return ret;
}

bool istft(const std::complex<float> *stft, const size_t ncolumns,
           const size_t nrows, const size_t hop_length, const size_t win_length,
           std::vector<float> *output, const bool center)
{
  size_t n_fft = 2 * (ncolumns -1);

  //std::cout << "win_length = " << win_length << std::endl;
  //std::cout << "n_fft = " << n_fft << std::endl;

  std::vector<float> _ifft_window;
  {
    bool ret = get_window("hann", win_length, &_ifft_window, /* periodic */true);
    if (!ret) {
      return false;
    }
  }

  // Pad out to match n_fft
  std::vector<float> ifft_window;
  {
    bool ret = pad_center(_ifft_window, n_fft, &ifft_window);
    if (!ret) {
      return false;
    }
  }

  // TODO(LTE): Implement

  size_t n_frames = nrows;
  size_t expected_signal_len = n_fft + hop_length * (n_frames - 1);

  // y = np.zeros(expected_signal_len, dtype=dtype)
  std::vector<float> y(expected_signal_len, 0.0f);

  /*
    for i in range(n_frames):
        sample = i * hop_length
        spec = stft_matrix[:, i].flatten()
        spec = np.concatenate((spec, spec[-2:0:-1].conj()), 0)
        ytmp = ifft_window * fft.ifft(spec).real

        y[sample:(sample + n_fft)] = y[sample:(sample + n_fft)] + ytmp
  */
  for (size_t i = 0; i < n_frames; i++) {
    size_t sample = i * hop_length;

    // spec = stft_matrix[:, i].flatten()
    // => extract each row

    // spec = np.concatenate((spec, spec[-2:0:-1].conj()), 0)
    // => reflect with conjugate.
    // e.g. spec [a, b, c, d, e]
    //
    // np.concatenate((spec, spec[-2:0:-1].conj()), 0)
    //
    // => [a, b, c, d, e, d.conj(), c.conj(), b.conj()]

    std::vector<std::complex<float>> spec(n_fft);

    // spec = stft_matrix[:, i].flatten()
    for (size_t k = 0; k < ncolumns; k++) {
      spec[k] = stft[i * ncolumns + k];
    }

    // concat: spec[-2:0:-1].conj
    for (size_t k = 0; k < ncolumns - 2; k++) {
      spec[ncolumns + k] = std::conj(spec[ncolumns - 2 - k]);
    }

    // ytmp = ifft_window * fft.ifft(spec).real

    std::vector<float> ytemp(n_fft);
    std::vector<std::complex<float>> ifft_spec;

    bool ret = ifft(spec.data(), n_fft, /* 1D */1, /* n */n_fft, &ifft_spec);

    if (!ret) {
      return false;
    }


    for (size_t k = 0; k < n_fft; k++) {
      ytemp[k] = ifft_window[k] * ifft_spec[k].real();
    }

    // y[sample:(sample + n_fft)] = y[sample:(sample + n_fft)] + ytmp
    for (size_t k = 0; k < n_fft; k++) {
      y[sample + k] += ytemp[k];
    }
  }

  // Normalize by sum of squared window
  std::vector<float> ifft_window_sum;
  {
    bool ret = hann_window_sumsquare(n_frames, hop_length, win_length, n_fft, &ifft_window_sum);
    if (!ret) {
      return false;
    }
  }

  // approx_nonzero_indices = ifft_window_sum > util.tiny(ifft_window_sum)
  // y[approx_nonzero_indices] /= ifft_window_sum[approx_nonzero_indices]
  assert(y.size() == ifft_window_sum.size());

  for (size_t i = 0; i < y.size(); i++) {
    if (y[i] > std::numeric_limits<float>::min()) {
      y[i] /= ifft_window_sum[i];
    }
  }

  //  if length is None:
  //      # If we don't need to control length, just do the usual center trimming
  //      # to eliminate padded data
  //      if center:
  //          y = y[int(n_fft // 2):-int(n_fft // 2)]

  // TODO(LTE): Support `length` parameter.
  if (center) {
    // eliminate padded data.
    const size_t s_idx = n_fft / 2;
    const size_t e_idx = size_t(int64_t(y.size()) - (int64_t(n_fft) / 2));

    output->resize(e_idx - s_idx);
    for (size_t i = s_idx; i < e_idx; i++) {
      (*output)[i - s_idx] = y[i];
    }
  } else {
    (*output) = y;
  }

  return true;
}

}  // namespace nanosnap
