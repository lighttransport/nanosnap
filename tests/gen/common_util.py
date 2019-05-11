# arr : numpy 1D or 2D array
def print_c_array(arr):
    c_arr = []

    if len(arr.shape) == 1:
        for i in arr:
            c_arr.append(str(float(i)))

        c_str = ', '.join(c_arr)
    elif len(arr.shape) == 2:
        for row in arr:
            for i in row:
                c_arr.append(str(float(i)))

        c_str = ', '.join(c_arr)
    else:
        # Unsupported.
        raise

    return c_str

