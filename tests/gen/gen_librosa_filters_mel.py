#!/usr/bin/env python

import numpy
import librosa

from common_util import print_c_array

def gen():
    sr = 22050
    n_fft = 2048
    M = librosa.filters.mel(sr, n_fft)

    if numpy.isfortran(M):
        M = numpy.transpose(M)

    print('const float g_reference[] = {' + print_c_array(M) + '};')
    print('const size_t k_sr = {};'.format(sr))
    print('const size_t k_nfft = {};'.format(n_fft))

gen()
