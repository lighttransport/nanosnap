#!/usr/bin/env python

import numpy

from common_util import print_c_array

def gen():
    n = 32
    seed = 42
    numpy.random.seed(seed)
    volume = numpy.random.rand(n).astype(numpy.float32)

    print('const float g_reference[] = {' + print_c_array(volume) + '};')
    print('const size_t k_n = {};'.format(n))
    print('const size_t k_seed = {};'.format(seed))

gen()
