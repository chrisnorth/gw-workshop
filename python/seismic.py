# Seismic isolation code

import scipy.io

def read_data(file):

    mat = scipy.io.loadmat(file)
    return mat