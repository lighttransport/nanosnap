#!/usr/bin/env python

import numpy
import scipy.signal

from common_util import print_c_array

def gen():

    win_length = 64
    y = scipy.signal.get_window('hann', win_length) # fftbins=True = peridoc

    print('const float g_reference[] = {' + print_c_array(y) + '};')
    print('const size_t k_win_length = {};'.format(win_length))

gen()
