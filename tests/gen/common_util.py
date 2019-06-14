import numpy

# NOTE(LTE): We cannot use numpy.iscomplex() whether input value is complex-value,
# since it will return False for complex-value whose imaginary part is zero.
#
# arr : numpy 1D or 2D array
def print_c_array(arr):
    c_arr = []

    if numpy.isfortran(arr):
        arr = numpy.transpose(arr)

    if len(arr.shape) == 1:
        for i in arr:
            if isinstance(i, complex) or isinstance(i, numpy.complex) or isinstance(i, numpy.complex64):
                c_arr.append(str(float(i.real)) + 'f')
                c_arr.append(str(float(i.imag)) + 'f')
            else:
                c_arr.append(str(float(i)) + 'f')

        c_str = ', '.join(c_arr)
    elif len(arr.shape) == 2:
        for row in arr:
            for i in row:
                if isinstance(i, complex) or isinstance(i, numpy.complex) or isinstance(i, numpy.complex64):
                    c_arr.append(str(float(i.real)) + 'f')
                    c_arr.append(str(float(i.imag)) + 'f')
                else:
                    c_arr.append(str(float(i)) + 'f')

        c_str = ', '.join(c_arr)
    else:
        # Unsupported.
        raise

    return c_str

