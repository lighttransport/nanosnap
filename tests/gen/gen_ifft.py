#!/usr/bin/env python

import numpy
from python_speech_features import *

from common_util import print_c_array

def gen():
    nframes = 16
    nrows = 4

    numpy.random.seed(42)
    volume = numpy.random.rand(nrows, nframes).astype(numpy.float32)

    fft_len = nframes-1
    f = numpy.fft.rfft(volume, fft_len)
    #print(len(y))
    #print(numpy.absolute(y))

    #a = numpy.absolute(y)
    #print(a.shape)

    #print('sum = ', numpy.sum(a, 1))

    n = 6
    y = numpy.fft.ifft(f, n)

    print('const float g_input[] = {' + print_c_array(f) + '};')
    print('const float g_reference[] = {' + print_c_array(y) + '};')
    print('const size_t k_n = {};'.format(n))
    print('const size_t k_ncolumns = {};'.format(f.shape[1]))
    print('const size_t k_nrows = {};'.format(f.shape[0]))

gen()
