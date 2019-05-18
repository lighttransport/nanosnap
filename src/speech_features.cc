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

#include <iostream>  // dbg

// pocketfft

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

//
// Frame a signal into overlapping frames
//
// sig: the audio signal to frame.
// sig_len: length of `sig`..
// _frame_len: length of each frame measured in samples
// _frame_step: length of each frame measured in samples
// winfunc: Window function(default =
//
static std::vector<float> framesig(const float *sig, const size_t sig_len, const float _frame_len,
              const float _frame_step,
              std::function<std::vector<float>(int)> winfunc = fIdentityWindowingFunction,
              const bool stride_trick = false) {
  // TODO(LTE): Implement;
  (void)winfunc;
  (void)stride_trick;

  // Fortunately, C++11 has `round` function for round half up.
  int frame_len = int(std::round(_frame_len));
  int frame_step = int(std::round(_frame_step));

  int num_frames;
  if (int(sig_len) <= frame_len) {
    num_frames = 1;
  } else {
    num_frames = 1 + int(std::ceil((float(sig_len) - float(frame_len)) / float(frame_step)));
  }

  int pad_len = (num_frames - 1) * frame_step + frame_len;

  std::vector<float> padded_signal;
  padded_signal.resize(sig_len);
  memcpy(padded_signal.data(), sig, sizeof(float) * sig_len);

  // apped zeros.
  for (int i = 0; i < pad_len; i++) {
    padded_signal.emplace_back(0.0f);
  }

  std::vector<float> frames; // length =
  std::vector<float> win; // length = frame_len
  if (stride_trick) {
    // evaluate windowing function.
    win = winfunc(frame_len);
  } else {
  }

  // frames * win
  assert(frames.size() == win.size());
  for (size_t i = 0; i < win.size(); i++) {
    frames[i] *= win[i];
  }

  return frames;
}

// Input is a 2D array. frames[D][N]
bool magspec(const size_t nframes, const size_t nrows, const float *frames,
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
  size_t output_num_elements = nrows * (NFFT / 2 + 1);
  complex_output.resize(output_num_elements);

  bool ret = rfft(frames, NFFT, nframes, nrows, complex_output.data());
  if (!ret) {
    return false;
  }

  output->resize(output_num_elements);

  // Take an absolute of complex value.
  for (size_t i = 0; i < output_num_elements; i++) {
    (*output)[i] = std::fabs(complex_output[i]);
  }

  return true;
}

bool powspec(const size_t nframes, const size_t nrows, const float *frames,
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

  bool ret = magspec(nframes, nrows, frames, NFFT, output);

  if (!ret) {
    return false;
  }

  // 1 / NFFT
  const float k = 1.0f / float(NFFT);
  for (size_t i = 0; i < output->size(); i++) {
    (*output)[i] *= k;
  }

  return true;
}

static std::vector<float> get_filterbanks(const int nfilt = 20, const int nfft = 512,
                                   const float samplerate = 16000, const float lowfreq = 0,
                                   const float highfreq = -1.0f) {
  // TODO(LTE): Implenment
  std::vector<float> fb;

  (void)nfilt;
  (void)nfft;
  (void)samplerate;
  (void)lowfreq;
  (void)highfreq;

  return fb;
}


}  // namespace

bool lifter(const float *cepstra, const size_t nframes, const size_t ncoeff,
            float *output, const int L) {
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

  if (L <= 0) {
    // Do nothing. Simply copy `cepstra`
    memcpy(output, cepstra, sizeof(float) * nframes * ncoeff);
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
      output[f * ncoeff + i] = lift[i] * cepstra[f * ncoeff + i];
    }
  }

  return true;
}

bool fbank(const float *_signal, const size_t nframes, const float samplerate,
           const float winlen, const float winstep,
           const std::function<std::vector<float>(int)> winfunc,
           const int nfilt, const int nfft, const float lowfreq,
           const float highfreq, const float preemph) {
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
    return false;
  }

  // input signal is 1D with length `nframes`.
  // signal = sigproc.preemphasis(signal,preemph)
  std::vector<float> signal;
  {
    // y[0] = x[0]
    // y[n] = x[n] - coeff * x[n-1] for n > 0
    signal.resize(nframes);

    signal[0] = _signal[0];
    for (size_t i = 1; i < nframes; i++) {
      signal[i] = _signal[i] - preemph * signal[i - 1];
    }
  }

  //  frames = sigproc.framesig(signal, winlen*samplerate, winstep*samplerate,
  //  winfunc)
  std::vector<float> frames = framesig(signal.data(), nframes, winlen * samplerate, winstep * samplerate, winfunc);

  size_t nrows = 1;  // 1D array

  // pspec = sigproc.powspec(frames,nfft)
  std::vector<float> pspec;
  {
    bool ret = powspec(nframes, nrows, frames.data(), size_t(nfft), &pspec);
    if (!ret) {
      return false;
    }
  }

  std::vector<float> feat;

  //  energy = numpy.sum(pspec,1) # this stores the total energy in each frame
  std::vector<float> energy;
  energy.resize(nrows);
  {
    for (size_t j = 0; j < nrows; j++) {
      energy[j] = 0.0f;
      for (size_t i = 0; i < nframes; i++) {
        energy[j] += pspec[j * nframes + i];
      }

      // energy = numpy.where(energy == 0,numpy.finfo(float).eps,energy)
      energy[j] = std::max(energy[j], std::numeric_limits<float>::epsilon());
    }
  }

  //  fb = get_filterbanks(nfilt,nfft,samplerate,lowfreq,highfreq)
  std::vector<float> fb =
      get_filterbanks(nfilt, nfft, samplerate, lowfreq, highfreq);

  //  feat = numpy.dot(pspec,fb.T) # compute the filterbank energies
  {
    for (size_t j = 0; j < nrows; j++) {
      for (size_t i = 0; i < nframes; i++) {
        feat[j * nframes + i] = 0.0f;
        for (size_t k = 0; k < nframes; k++) {
          // we don't need to compute transpose of `fb`
          feat[j * nframes + i] += pspec[j * nframes + k] * fb[j * nframes + k];
        }
      }
    }

    // avoid zero
    for (size_t j = 0; j < feat.size(); j++) {
      feat[j] = std::max(feat[j], std::numeric_limits<float>::epsilon());
    }
  }

  return true;
}

#if 0
bool mfcc(const std::vector<float> &signal, std::vector<float> *output, const float winlen, const size_t nfft,
          const bool calculate_nfft, const bool append_energy) {
  // Based on python_speech_features
  // nfft = nfft or calculate_nfft(samplerate, winlen)
  // feat,energy =
  // fbank(signal,samplerate,winlen,winstep,nfilt,nfft,lowfreq,highfreq,preemph,winfunc)
  // feat = numpy.log(feat)
  // feat = dct(feat, type=2, axis=1, norm='ortho')[:,:numcep]
  // feat = lifter(feat,ceplifter)
  // if appendEnergy: feat[:,0] = numpy.log(energy) # replace first cepstral
  // coefficient with log of frame energy return feat

  size_t nfft;
  if (calculate_nfft) {
    nfft = CalculateNfft(samplerate, winlen)
  } else {
    nfft = _nfft;
  }

#if 0
  std::vector<float> feat;
  std::vector<float> energy;

  feat, energy = fbank(signal, samplerate, winlen, winstep, nfilt, nfft,
                       lowfreq, highfreq, preemph, winfunc)

        // Take a log.
        for (size_t i = 0; i < feat.size(); t++) {
    feat[i] = std::log(feat[i]);
  }

  // feat = dct(feat, type=2, axis=1, norm='ortho')[:,:numcep]
  // feat = lifter(feat,ceplifter)

  if (appendEnergy) {
    // replace first cepstral coefficient with log of frame energy
    // feat[:,0] = numpy.log(energy) # replace first cepstral coefficient with
    // log of frame energy
  }

  return feat
#endif
}
#endif

}  // namespace nanosnap
