#!/usr/bin/env python

import scipy.signal
import numpy

# arr : 1D array
def print_c_array(arr):
    c_arr = []
    for i in arr:
        c_arr.append(str(i))

    c_str = ', '.join(c_arr)

    return c_str

def gen():
    volume = numpy.random.rand(32)
    out = scipy.signal.medfilt(volume, kernel_size=3)

    print('const float input[] = {' + print_c_array(volume) + '};')
    print('const float reference[] = {' + print_c_array(out) + '};')

gen()
