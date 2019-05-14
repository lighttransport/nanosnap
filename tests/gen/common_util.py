# arr : numpy 1D or 2D array
def print_c_array(arr):
    c_arr = []

    if len(arr.shape) == 1:
        for i in arr:
            if isinstance(i, complex):
                c_arr.append(str(float(i.real)) + 'f')
                c_arr.append(str(float(i.imag)) + 'f')
            else:
                c_arr.append(str(float(i)) + 'f')

        c_str = ', '.join(c_arr)
    elif len(arr.shape) == 2:
        for row in arr:
            for i in row:
                if isinstance(i, complex):
                    c_arr.append(str(float(i.real)) + 'f')
                    c_arr.append(str(float(i.imag)) + 'f')
                else:
                    c_arr.append(str(float(i)) + 'f')

        c_str = ', '.join(c_arr)
    else:
        # Unsupported.
        raise

    return c_str

