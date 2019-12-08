#include <algorithm>
#include <limits>

// #include <iostream> // DBG

#include "nanosnap/nanosnap.h"

namespace nanosnap {

namespace {

template <typename T>
inline T safe_div(const T a, const T b) {
  if (std::fabs(b) < std::numeric_limits<T>::epsilon()) {
    return static_cast<T>(0.0);
  } else {
    return a / b;
  }
}

// numpy.linspace in C++
static std::vector<float> linspace(const float start, const float stop,
                                   const size_t num = 50,
                                   const bool end_point = true) {
  std::vector<float> w(num);

  int denom = end_point ? (int(num) - 1) : int(num);

  float step = safe_div(stop - start, float(denom));

  for (size_t i = 0; i < num; i++) {
    w[i] = i * step;
  }

  return w;
}

std::vector<float> fft_frequencies(const float sr, const int n_fft) {
  return linspace(/* start */ 0.0f, /* end */ sr / 2.0f,
                  /* num */ size_t(1 + n_fft / 2), /* endpoint */ true);
}

float hz_to_mel(const float freq, bool htk = false) {
  if (htk) {
    return hz2mel(freq);
  }

  // Fill in the linear part
  const float f_min = 0.0f;
  const float f_sp = 200.0f / 3.0f;

  float mels = (freq - f_min) / f_sp;

  // Fill in the log-scale part

  const float min_log_hz = 1000.0f;  // beginning of log region (Hz)
  const float min_log_mel = (min_log_hz - f_min) / f_sp;  // same (Mels)
  const float logstep = std::log(6.4f) / 27.0f;  // step size for log region

  if (freq >= min_log_hz) {
    mels = min_log_mel + std::log(freq / min_log_hz) / logstep;
  }

  return mels;
}

float mel_to_hz(const float mel, bool htk = false) {
  if (htk) {
    return mel2hz(mel);
  }

  // Fill in the linear scale
  const float f_min = 0.0f;
  const float f_sp = 200.0f / 3.0f;
  float freqs = f_min + f_sp * mel;

  // And now the nonlinear scale
  const float min_log_hz = 1000.0f;  // beginning of log region (Hz)
  const float min_log_mel = (min_log_hz - f_min) / f_sp;  // same (Mels)
  const float logstep = std::log(6.4f) / 27.0f;  // step size for log region

  if (mel >= min_log_mel) {
    freqs = min_log_hz * std::exp(logstep * (mel - min_log_mel));
  }

  return freqs;
}

std::vector<float> mel_frequencies(const int n_mels = 128,
                                   const float fmin = 0.0f,
                                   const float fmax = 11025.0f,
                                   const bool htk = false) {
  // 'Center freqs' of mel bands - uniformly spaced between limits
  const float min_mel = hz_to_mel(fmin, htk);
  const float max_mel = hz_to_mel(fmax, htk);

  std::vector<float> mels = linspace(min_mel, max_mel, size_t(n_mels));

  std::transform(mels.begin(), mels.end(), mels.begin(),
                 [htk](float mel) { return mel_to_hz(mel, htk); });

  return mels;
}

}  // namespace

bool mel_filter(const float sample_rate, const int n_fft, std::vector<float> *M,
                const int n_mels, const float fmin, const float _fmax,
                const bool htk, const bool norm) {
  float fmax = _fmax;
  if (fmax < 0.0f) {
    fmax = sample_rate / 2.0f;
  }

  if (n_mels < 1) {
    return false;
  }

  size_t height = size_t(n_mels);
  size_t width = size_t(1 + n_fft / 2);

  if (width < 2) {
    // TOO short
    return false;
  }

  std::vector<float> weights(width * height);

  //  # Center freqs of each FFT bin. 1D
  const std::vector<float> fftfreqs = fft_frequencies(sample_rate, n_fft);

  //  # 'Center freqs' of mel bands - uniformly spaced between limits. 1D
  const std::vector<float> mel_f = mel_frequencies(n_mels + 2, fmin, fmax, htk);

  // fdif = np.diff(mel_f). 1D
  std::vector<float> fdiff(mel_f.size() - 1);

  for (size_t x = 0; x < fdiff.size(); x++) {
    fdiff[x] = mel_f[(x + 1)] - mel_f[x];
  }

  //  ramps = np.subtract.outer(mel_f, fftfreqs)
  std::vector<float> ramps(mel_f.size() * fftfreqs.size());  // 2D

  //std::cout << "fftfreqs..size" << fftfreqs.size() << ", width = " << width << "\n";
  for (size_t i = 0; i < mel_f.size(); i++) {
    for (size_t j = 0; j < fftfreqs.size(); j++) {
      ramps[i * fftfreqs.size() + j] = mel_f[i] - fftfreqs[j];
    }
  }

  for (size_t y = 0; y < height; y++) {
    // lowwer and upper slopes for all bins.
    for (size_t x = 0; x < width; x++) {
      const float lower = -ramps[y * width + x] / fdiff[y];
      const float upper = ramps[(y + 2) * width + x] /
                          fdiff[y + 1];  // +2 is safe since we create `mel_f`
                                         // with `n_mels + 2` size

      weights[y * width + x] = std::max(0.0f, std::min(lower, upper));
    }
  }

  if (norm) {
    // Slaney-style mel is scaled to be approx constant energy per channel
    // enorm = 2.0 / (mel_f[2:n_mels+2] - mel_f[:n_mels])
    for (size_t i = 0; i < size_t(n_mels); i++) {
      float enorm = 2.0f / (mel_f[i + 2] - mel_f[i]);
      // std::cout << "[" << i << "] enorm = " << enorm << "\n";

      for (size_t x = 0; x < width; x++) {
        weights[i * width + x] *= enorm;
      }
    }
  }

  // only check weights if f_mel[0] is positive
  {
    // if not np.all((mel_f[:-2] == 0) | (weights.max(axis=1) > 0))
    for (size_t y = 0; y < height; y++) {
      float maxval = weights[y * width];
      for (size_t x = 1; x < width; x++) {
        maxval = std::max(weights[y * width + x], maxval);
      }

      if ((maxval > 0.0f) ||
          (std::fabs(mel_f[y]) < std::numeric_limits<float>::epsilon())) {
        // ok
      } else {
        // TODO(LTE): print warn
        // warnings.warn('Empty filters detected in mel frequency basis. '
        //               'Some channels will produce empty responses. '
        //               'Try increasing your sampling rate (and fmax) or '
        //               'reducing n_mels.')
        break;
      }
    }
  }

  (*M) = std::move(weights);

  return true;
}

}  // namespace nanosnap
