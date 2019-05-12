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
#include <memory>

#include <iostream> // dbg

namespace nanosnap {

namespace {

// numpy.linspace in C++
static void linspace(const float start, const float stop, std::vector<float> *out, const size_t num = 50, const bool end_point = true)
{
  out->resize(size_t(num));

  int denom = end_point ? (num - 1) : num;

  float step = (stop - start) / float(denom);

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

//
// Frame a signal into overlapping frames
//
// sig: the audio signal to frame.
// _frame_len: length of each frame measured in samples
// _frame_step: length of each frame measured in samples
//
void framesig(const float *sig, const float _frame_len, const float _frame_step)
{
  // TODO(LTE): Implement;
  (void)sig;

  // Fortunately, C++11 has `round` function for round half up.
  float iframe_len = std::round(_frame_len);
  float iframe_step = std::round(_frame_step);
  (void)iframe_len;
  (void)iframe_step;
}

#if 0
// Input is a 2D array. frames[D][N]
bool magspec(const size_t N, const size_t D, const double *frames, const size_t NFF, std::vector<double> *output)
{
  /*
    """Compute the magnitude spectrum of each frame in frames. If frames is an NxD matrix, output will be Nx(NFFT/2+1).

    :param frames: the array of frames. Each row is a frame.
    :param NFFT: the FFT length to use. If NFFT > frame_len, the frames are zero-padded.
    :returns: If frames is an NxD matrix, output will be Nx(NFFT/2+1). Each row will be the magnitude spectrum of the corresponding frame.
    """
    if numpy.shape(frames)[1] > NFFT:
        logging.warn(
            'frame length (%d) is greater than FFT size (%d), frame will be truncated. Increase NFFT to avoid.',
            numpy.shape(frames)[1], NFFT)
    complex_spec = numpy.fft.rfft(frames, NFFT)
    return numpy.absolute(complex_spec)
  */

  if (D > NFFT) {
    // TODO(LTE): Use zero-adding and continue
    return false;
  }

  output->resize(N * (NFFT/2+1));

  // numpy.absolute
  for (size_t i < 0; i < data_length; i++) {
    (*output)[i] = std::fabs((*output)[i]);
  }
}


void powspec(const size_t frame_len, size_t const float *frames, const size_t NFFT) {
    /*
    """Compute the power spectrum of each frame in frames. If frames is an NxD matrix, output will be Nx(NFFT/2+1).

    :param frames: the array of frames. Each row is a frame.
    :param NFFT: the FFT length to use. If NFFT > frame_len, the frames are zero-padded.
    :returns: If frames is an NxD matrix, output will be Nx(NFFT/2+1). Each row will be the power spectrum of the corresponding frame.
    """
    return 1.0 / NFFT * numpy.square(magspec(frames, NFFT))
    */

}
#endif



#define kPI (3.141592f)

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

//void get_filterbanks(const int nfilt=20, const int nfft=512,samplerate=16000,lowfreq=0,highfreq=None)
//{
//
//}

#if 0
bool fbank(const size_t n, const float *_signal, const float samplerate, const float winlen, const float winstep, const float nfilt, const float nfft, const float lowfreq, const float _highfreq, const float preemph)
{
  /*
    highfreq= highfreq or samplerate/2
    signal = sigproc.preemphasis(signal,preemph)
    frames = sigproc.framesig(signal, winlen*samplerate, winstep*samplerate, winfunc)
    pspec = sigproc.powspec(frames,nfft)
    energy = numpy.sum(pspec,1) # this stores the total energy in each frame
    energy = numpy.where(energy == 0,numpy.finfo(float).eps,energy) # if energy is zero, we get problems with log

    fb = get_filterbanks(nfilt,nfft,samplerate,lowfreq,highfreq)
    feat = numpy.dot(pspec,fb.T) # compute the filterbank energies
    feat = numpy.where(feat == 0,numpy.finfo(float).eps,feat) # if feat is zero, we get problems with log

    return feat,energy
  */
  if (samplerate <= 0.0f) {
    // invalid parameter
    return false;
  }

  float highfreq = _highfreq;
  if (highfreq < 0.0f) {
    // Use samplerate /2
    highfreq = samplerate / 2.0f;
  }

  // preemphasis. No emphasis when preemph == 0
  // y[0] = x[0]
  // y[n] = x[n] - coeff * x[n-1] for n > 0
  StackVector<float, 512> signal;
  signal->resize(n);

  signal[0] = _signal[0];
  for (size_t i = 1; i < n; i++) {
    signal[i] = _signal[i] - preemph * signal[i-1];
  }

  StackVector<float, 512> frames;
  FrameSig();

  StackVector<float, 512> pspec;

  StackVector<float, 512> feat;


  //feat = numpy.dot(pspec,fb.T) # compute the filterbank energies

  return true;


}
#endif

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
