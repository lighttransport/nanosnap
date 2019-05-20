#!/usr/bin/env python

import numpy

from common_util import print_c_array

def gen():
    n = 5
    m = 3

    seed = 42
    numpy.random.seed(seed)
    a = numpy.random.rand(n).astype(numpy.float32)
    v = numpy.random.rand(m).astype(numpy.float32)
    c = numpy.convolve(a, v, mode='valid')

    assert len(c) == (n - m + 1)

    print('const float g_a[] = {' + print_c_array(a) + '};')
    print('const float g_v[] = {' + print_c_array(v) + '};')
    print('const float g_reference[] = {' + print_c_array(c) + '};')
    print('const size_t k_n = {};'.format(n))
    print('const size_t k_m = {};'.format(m))
    print('const size_t k_outnum = {};'.format(len(c)))

gen()
