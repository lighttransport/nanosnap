#!/usr/bin/env python

import scipy.signal
import numpy

from common_util import print_c_array

def gen():
    n = 32
    ksize = 3
    numpy.random.seed(42)
    volume = numpy.random.rand(n).astype(numpy.float32)
    out = scipy.signal.medfilt(volume, kernel_size=ksize)

    print('const float g_input[] = {' + print_c_array(volume) + '};')
    print('const float g_reference[] = {' + print_c_array(out) + '};')
    print('const size_t k_input_n = {};'.format(n))
    print('const size_t k_window_size = {};'.format(ksize))

gen()
