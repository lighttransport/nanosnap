#!/usr/bin/env python

import numpy
import librosa

from common_util import print_c_array

def gen():
    nframes = 128
    numpy.random.seed(42)
    y = numpy.random.rand(nframes).astype('float32')
    n_fft = 64
    win_length = n_fft  # by default win_length = n_fft
    hop_length = win_length // 4 # default hop length is win_length / 4
    # window function = 'hann'
    D = librosa.stft(y, n_fft=n_fft, hop_length=hop_length, win_length=win_length)
    #print(D.shape)

    # NOTE: librosa uses Fortran order for stft matrix.
    if numpy.isfortran(D):
        D = numpy.transpose(D)

    print('const float g_input[] = {' + print_c_array(y) + '};')
    print('const float g_reference[] = {' + print_c_array(D) + '};')
    print('const size_t k_n_fft = {};'.format(n_fft))
    print('const size_t k_nframes = {};'.format(nframes))
    print('const size_t k_out_dim0 = {};'.format(D.shape[0]))
    print('const size_t k_out_dim1 = {};'.format(D.shape[1]))
    #print('const size_t k_nrows = {};'.format(1))
    print('const size_t k_hop_length = {};'.format(hop_length))
    print('const size_t k_win_length = {};'.format(win_length))

gen()
