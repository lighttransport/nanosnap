#!/usr/bin/env python

import numpy
from python_speech_features import *

from common_util import print_c_array

def gen():
    nframes = 16
    ncoeffs = 5

    numpy.random.seed(42)
    volume = numpy.random.rand(nframes, ncoeffs).astype(numpy.float32)

    y = lifter(volume)

    print('const float g_input[] = {' + print_c_array(volume) + '};')
    print('const float g_reference[] = {' + print_c_array(y) + '};')
    print('const size_t k_nframes = {};'.format(nframes))
    print('const size_t k_ncoeffs = {};'.format(ncoeffs))

gen()
