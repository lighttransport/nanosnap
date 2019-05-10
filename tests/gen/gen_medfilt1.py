#!/usr/bin/env python

import scipy.signal
import numpy

# arr : 1D array
def print_c_array(arr):
    c_arr = []
    for i in arr:
        c_arr.append(str(float(i)))

    c_str = ', '.join(c_arr)

    return c_str

def gen():
    n = 32
    ksize = 3
    volume = numpy.random.rand(n).astype(numpy.float32)
    out = scipy.signal.medfilt(volume, kernel_size=ksize)

    print('const float input[] = {' + print_c_array(volume) + '};')
    print('const float reference[] = {' + print_c_array(out) + '};')
    print('const size_t input_n = {};'.format(n))
    print('const size_t window_size = {};'.format(ksize))

gen()
