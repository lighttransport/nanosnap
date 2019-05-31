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

#include "stack_vector.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <algorithm>

#include <iostream>  // dbg

namespace nanosnap {

namespace {

constexpr float kPI = 3.141592f;

#if 0
template <typename T>
inline T safe_div(const T a, const T b) {
  if (std::fabs(b) < std::numeric_limits<T>::epsilon()) {
    return static_cast<T>(0.0);
  } else {
    return a / b;
  }
}

// numpy.linspace in C++
static void linspace(const float start, const float stop,
                     std::vector<float> *out, const size_t num = 50,
                     const bool end_point = true) {
  out->resize(size_t(num));

  int denom = end_point ? (int(num) - 1) : int(num);

  float step = safe_div(stop - start, float(denom));

  for (size_t i = 0; i < num; i++) {
    (*out)[i] = i * step;
  }
}

// samplerate : [Hz]
// winlen : [s]
size_t calcualte_nfft(float samplerate, const float winlen) {
  (void)samplerate;
  (void)winlen;
  return 0;
}
#endif

static std::vector<float> preemphasis(const float *signal, size_t sig_len,
                                      float coeff) {
  std::vector<float> preemph;
  preemph.resize(sig_len);

  for (size_t i = sig_len - 1; i >= 1; i--) {
    preemph[i] = signal[i] - signal[i - 1] * coeff;
  }
  preemph[0] = signal[0];

  return preemph;
}

static std::vector<float> framesig(
    const float *signal, const int sig_len, const int frame_len,
    const int padded_frame_len, const int frame_step,
    std::function<float(int)> winfunc = fIdentityWindowingFunction) {
  // Based on c_speech_features
  int frame_width = (std::max)(padded_frame_len, frame_len);

  int nframes;
  if (sig_len > frame_len) {
    nframes = 1 + int(std::ceil((sig_len - frame_len) / float(frame_step)));
  } else {
    nframes = 1;
  }

  std::vector<int> indices;
  indices.resize(size_t(nframes * frame_len));

  for (int i = 0, idx = 0; i < nframes; i++) {
    int base = i * frame_step;
    for (int j = 0; j < frame_len; j++, idx++) {
      indices[size_t(idx)] = base + j;
    }
  }

  std::vector<float> frames;
  frames.resize(size_t(nframes * frame_width));

  for (int i = 0, idx = 0, iidx = 0; i < nframes; i++) {
    for (int j = 0; j < frame_len; j++, idx++, iidx++) {
      int index = indices[size_t(iidx)];
      frames[size_t(idx)] = index < sig_len ? signal[size_t(index)] : 0.0f;

      // apply windowing function
      frames[size_t(idx)] *= winfunc(j);
    }
    for (int j = frame_len; j < padded_frame_len; j++, idx++) {
      frames[size_t(idx)] = 0.0;
    }
  }

  return frames;
}

// Input is a 2D array. frames[nrows][nframes]
bool magspec(const float *frames, const size_t nframes, const size_t nrows,
             const size_t NFFT, std::vector<float> *output) {
  /*
    """Compute the magnitude spectrum of each frame in frames. If frames is an
    NxD matrix, output will be Nx(NFFT/2+1).

    :param frames: the array of frames. Each row is a frame.
    :param NFFT: the FFT length to use. If NFFT > frame_len, the frames are
    zero-padded. :returns: If frames is an NxD matrix, output will be
    Nx(NFFT/2+1). Each row will be the magnitude spectrum of the corresponding
    frame.
    """
    if numpy.shape(frames)[1] > NFFT:
        logging.warn(
            'frame length (%d) is greater than FFT size (%d), frame will be
    truncated. Increase NFFT to avoid.', numpy.shape(frames)[1], NFFT)
    complex_spec = numpy.fft.rfft(frames, NFFT)
    return numpy.absolute(complex_spec)
  */

  if (nframes > NFFT) {
    // TODO(LTE): Use zero-adding and continue
    return false;
  }

  std::vector<std::complex<float>> complex_output;

  bool ret = rfft(frames, NFFT, nframes, nrows, &complex_output);
  if (!ret) {
    return false;
  }

  output->resize(complex_output.size());

  // Take an absolute of complex value.
  for (size_t i = 0; i < complex_output.size(); i++) {
    (*output)[i] = std::abs(complex_output[i]);
  }

  return true;
}

bool powspec(const float *frames, const size_t nframes, const size_t nrows,
             const size_t NFFT, std::vector<float> *output) {
  /*
  """Compute the power spectrum of each frame in frames. If frames is an NxD
  matrix, output will be Nx(NFFT/2+1).

  :param frames: the array of frames. Each row is a frame.
  :param NFFT: the FFT length to use. If NFFT > frame_len, the frames are
  zero-padded. :returns: If frames is an NxD matrix, output will be
  Nx(NFFT/2+1). Each row will be the power spectrum of the corresponding frame.
  """
  return 1.0 / NFFT * numpy.square(magspec(frames, NFFT))
  */

  bool ret = magspec(frames, nframes, nrows, NFFT, output);

  if (!ret) {
    return false;
  }

  // Scale by `1 / NFFT`
  const float k = 1.0f / float(NFFT);
  for (size_t i = 0; i < output->size(); i++) {
    (*output)[i] *= k;
  }

  return true;
}

static std::vector<float> get_filterbanks(const int nfilt, const int nfft,
                                          const float samplerate,
                                          const float lowfreq,
                                          const float highfreq) {
  // Based on c_speech_features.
  int feat_width = nfft / 2 + 1;
  float lowmel = hz2mel(lowfreq);
  float highmel = hz2mel((highfreq <= lowfreq) ? samplerate / 2.0f : highfreq);
  StackVector<int, 32> bin;
  bin->resize(size_t(nfilt + 2));

  std::vector<float> fbank(size_t(nfilt * feat_width));

  for (size_t i = 0; i < size_t(nfilt + 2); i++) {
    const float melpoint = ((highmel - lowmel) / float(nfilt + 1) * i) + lowmel;
    bin[i] = int(std::floor(((nfft + 1) * mel2hz(melpoint) / samplerate)));
  }

  for (size_t i = 0, idx = 0; i < size_t(nfilt);
       i++, idx += size_t(feat_width)) {
    int start = std::min(bin[i], bin[i + 1]);
    int end = std::max(bin[i], bin[i + 1]);
    for (int j = start; j < end; j++) {
      fbank[idx + size_t(j)] = (j - bin[i]) / float(bin[i + 1] - bin[i]);
    }
    start = std::min(bin[i + 1], bin[i + 2]);
    end = std::max(bin[i + 1], bin[i + 2]);
    for (int j = start; j < end; j++) {
      fbank[idx + size_t(j)] =
          (bin[i + 2] - j) / float(bin[i + 2] - bin[i + 1]);
    }
  }

  return fbank;
}

}  // namespace

bool lifter(const float *cepstra, const size_t nframes, const size_t ncoeff,
            std::vector<float> *output, const int L) {
  // """Apply a cepstral lifter the the matrix of cepstra. This has the effect
  // of increasing the magnitude of the high frequency DCT coeffs.

  // :param cepstra: the matrix of mel-cepstra, will be numframes * numcep in
  // size. :param L: the liftering coefficient to use. Default is 22. L <= 0
  // disables lifter.
  // """
  // if L > 0:
  //     nframes,ncoeff = numpy.shape(cepstra)
  //     n = numpy.arange(ncoeff)
  //     lift = 1 + (L/2.)*numpy.sin(numpy.pi*n/L)
  //     return lift*cepstra
  // else:
  //     # values of L <= 0, do nothing
  //     return cepstra

  if (ncoeff == 0) {
    return false;
  }

  if (nframes == 0) {
    return false;
  }

  output->resize(nframes * ncoeff);

  if (L <= 0) {
    // Do nothing. Simply copy `cepstra`
    memcpy(output->data(), cepstra, sizeof(float) * nframes * ncoeff);
    return true;
  }

  StackVector<float, 128> lift;
  lift->resize(ncoeff);

  for (size_t i = 0; i < ncoeff; i++) {
    float n = float(i);
    lift[i] = 1.0f + (float(L) / 2.0f) * std::sin(kPI * n / float(L));
  }

  for (size_t f = 0; f < nframes; f++) {
    for (size_t i = 0; i < ncoeff; i++) {
      (*output)[f * ncoeff + i] = lift[i] * cepstra[f * ncoeff + i];
    }
  }

  return true;
}

int64_t fbank(const float *_signal, const size_t sig_len,
              const float samplerate, const float winlen, const float winstep,
              const std::function<float(int)> winfunc, const int nfilt,
              const int nfft, const float lowfreq, const float highfreq,
              const float preemph, std::vector<float> *features,
              std::vector<float> *energies) {
  /*
    highfreq= highfreq or samplerate/2
    signal = sigproc.preemphasis(signal,preemph)
    frames = sigproc.framesig(signal, winlen*samplerate, winstep*samplerate,
    winfunc) pspec = sigproc.powspec(frames,nfft) energy = numpy.sum(pspec,1) #
    this stores the total energy in each frame energy = numpy.where(energy ==
    0,numpy.finfo(float).eps,energy) # if energy is zero, we get problems with
    log

    fb = get_filterbanks(nfilt,nfft,samplerate,lowfreq,highfreq)
    feat = numpy.dot(pspec,fb.T) # compute the filterbank energies
    feat = numpy.where(feat == 0,numpy.finfo(float).eps,feat) # if feat is zero,
    we get problems with log

    return feat,energy
  */
  if (samplerate <= 0.0f) {
    // invalid parameter
    return -1;
  }

  if (!features) {
    return -2;
  }

  // input signal is 1D with length `nframes`.
  // signal = sigproc.preemphasis(signal,preemph)
  std::vector<float> signal;
  {
    // y[0] = x[0]
    // y[n] = x[n] - coeff * x[n-1] for n > 0
    signal.resize(sig_len);

    signal[0] = _signal[0];
    for (size_t i = 1; i < sig_len; i++) {
      signal[i] = _signal[i] - preemph * signal[i - 1];
    }
  }

  int frame_len = int(std::round(winlen * samplerate));
  int frame_step = int(std::round(winstep * samplerate));

  //  frames = sigproc.framesig(signal, winlen*samplerate, winstep*samplerate,
  //  winfunc)
  std::vector<float> frames =
      framesig(signal.data(), int(sig_len), frame_len,
               /* padded framelen */ nfft, frame_step, winfunc);

  size_t nframes = frames.size();

  // pspec = sigproc.powspec(frames,nfft)
  std::vector<float> pspec;
  {
    bool ret =
        powspec(frames.data(), nframes, /* nrows */ 1, size_t(nfft), &pspec);
    if (!ret) {
      return false;
    }
  }

  const size_t feature_width = size_t((nfft / 2) + 1);

  //  energy = numpy.sum(pspec,1) # this stores the total energy in each frame
  if (energies) {
    energies->resize(nframes);
    {
      for (size_t j = 0; j < nframes; j++) {
        float energy = 0.0f;
        for (size_t i = 0; i < feature_width; i++) {
          energy += pspec[j * feature_width + i];
        }

        // energy = numpy.where(energy == 0,numpy.finfo(float).eps,energy)
        (*energies)[j] = std::max(energy, std::numeric_limits<float>::min());
      }
    }
  }

  //  fb = get_filterbanks(nfilt,nfft,samplerate,lowfreq,highfreq)
  std::vector<float> fb =
      get_filterbanks(nfilt, nfft, samplerate, lowfreq, highfreq);

  //  feat = numpy.dot(pspec,fb.T) # compute the filterbank energies
  {
    for (size_t i = 0, idx = 0, pidx = 0; i < size_t(nframes);
         i++, idx += size_t(nfilt), pidx += feature_width) {
      for (size_t j = 0, fidx = 0; j < size_t(nfilt); j++) {
        float feat = 0.0f;
        for (size_t k = 0; k < feature_width; k++, fidx++) {
          feat += pspec[pidx + k] * fb[fidx];
        }

        // avoid zero
        (*features)[idx + j] =
            std::max(feat, std::numeric_limits<float>::min());
      }
    }
  }

  return int64_t(nframes);
}

bool mfcc(const float *signal, const size_t sig_len, const float samplerate,
          const float winlen, const float winstep, const size_t ncep,
          const size_t nfilt, const size_t nfft, const float low_freq,
          const float high_freq, const float preemph, const size_t cep_lifter,
          const bool append_energy, const std::function<float(int)> winfunc,
          std::vector<float> *out_mfcc) {
  //
  // nfft = nfft or calculate_nfft(samplerate, winlen)
  // feat,energy =
  // fbank(signal,samplerate,winlen,winstep,nfilt,nfft,lowfreq,highfreq,preemph,winfunc)
  // feat = numpy.log(feat)
  // feat = dct(feat, type=2, axis=1, norm='ortho')[:,:numcep]
  // feat = lifter(feat,ceplifter)
  // if appendEnergy: feat[:,0] = numpy.log(energy) # replace first cepstral
  // coefficient with log of frame energy return feat

  std::vector<float> feat;
  std::vector<float> energy;

  int64_t nframes = fbank(signal, sig_len, samplerate, winlen, winstep, winfunc,
                          int(nfilt), int(nfft), low_freq, high_freq, preemph,
                          &feat, append_energy ? &energy : nullptr);

  if (nframes <= 0) {
    return false;
  }

  // DCT is based on c_speech_features code.

  // Allocate an array so we can calculate the inner loop multipliers
  // in the DCT-II just one time.
  std::vector<float> dct2f;
  dct2f.resize(nfilt * ncep);

  // Perform DCT-II
  float sf1 = std::sqrt(1 / (4 * float(nfilt)));
  float sf2 = std::sqrt(1 / (2 * float(nfilt)));

  std::vector<float> mfcc;
  mfcc.resize(size_t(nframes) * ncep);

  for (size_t i = 0, idx = 0, fidx = 0; i < size_t(nframes);
       i++, idx += ncep, fidx += nfilt) {
    for (size_t j = 0, didx = 0; j < ncep; j++) {
      float sum = 0.0;
      for (size_t k = 0; k < nfilt; k++, didx++) {
        if (i == 0) {
          dct2f[didx] = std::cos(kPI * j * (2 * k + 1) / float(2 * nfilt));
        }
        sum += float(feat[fidx + k]) * dct2f[didx];
      }
      mfcc[idx + j] = float(sum * 2.0f * ((i == 0 && j == 0) ? sf1 : sf2));
    }
  }

  // Apply a cepstral lifter
  if (cep_lifter > 0) {
    bool ret = lifter(mfcc.data(), size_t(nframes), size_t(ncep), out_mfcc,
                      int(cep_lifter));
    if (!ret) {
      return false;
    }
  } else {
    (*out_mfcc) = mfcc;
  }

  // Append energies
  if (append_energy) {
    for (size_t i = 0, idx = 0; i < size_t(nframes); i++, idx += ncep) {
      (*out_mfcc)[idx] = std::log(energy[i]);
    }
  }

  return true;
}

int64_t ssc(const float *signal, const size_t sig_len, const int samplerate,
            const float winlen, const float winstep, const int nfilt,
            const int nfft, const int low_freq, const int high_freq,
            const float preemph_coeff, std::function<float(int)> winfunc,
            std::vector<float> *features) {
  // based on c_speech_features
  std::vector<float> preemph = preemphasis(signal, sig_len, preemph_coeff);

  int frame_len = int(std::round(winlen * samplerate));
  int frame_step = int(std::round(winstep * samplerate));
  int feat_width = nfft / 2 + 1;

  // Frame the signal into overlapping frames
  std::vector<float> frames = framesig(preemph.data(), int(sig_len), frame_len,
                                       nfft, frame_step, winfunc);

  size_t nframes = frames.size();

  // Compute the power spectrum of the frames
  std::vector<float> pspec;
  {
    if (!powspec(frames.data(), nframes, /* nrows */ 1, size_t(nfft), &pspec)) {
      return -1;
    }
  }

  // Make sure there are no zeroes in the power spectrum
  for (size_t i = 0, idx = 0; i < nframes; i++) {
    for (size_t j = 0; j < size_t(feat_width); j++, idx++) {
      pspec[idx] = std::max(pspec[idx], std::numeric_limits<float>::min());
    }
  }

  // Compute the filter-bank energies
  std::vector<float> fbank =
      get_filterbanks(nfilt, nfft, float(samplerate), float(low_freq), float(high_freq));
  std::vector<float> feat;
  feat.resize(nframes * size_t(nfilt));

  for (size_t i = 0, idx = 0, pidx = 0; i < nframes;
       i++, idx += size_t(nfilt), pidx += size_t(feat_width)) {
    for (size_t j = 0, fidx = 0; j < size_t(nfilt); j++) {
      for (size_t k = 0; k < size_t(feat_width); k++, fidx++) {
        feat[idx + j] += pspec[pidx + k] * fbank[fidx];
      }
    }
  }

  // Calculate Spectral Sub-band Centroid features
  std::vector<float> ssc;
  ssc.resize(nframes * size_t(nfilt));

  const float r = ((samplerate / 2) - 1) / float(feat_width - 1);
  for (size_t i = 0, idx = 0, pidx = 0; i < nframes;
       i++, idx += size_t(nfilt), pidx += size_t(feat_width)) {
    for (size_t j = 0, fidx = 0; j < size_t(nfilt); j++) {
      float R = 1;
      for (size_t k = 0; k < size_t(feat_width); k++, fidx++) {
        ssc[idx + j] += pspec[pidx + k] * R * fbank[fidx];
        R += r;
      }
      ssc[idx + j] /= feat[idx + j];
    }
  }

  (*features) = std::move(ssc);

  return int64_t(nframes);
}

}  // namespace nanosnap
